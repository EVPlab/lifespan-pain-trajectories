# ===================================================================
# Lifespan Pain Trajectories - IPD Model with Analytical Statistics
# FINAL CORRECTED VERSION 3 (with Deployment Exports)
# Updated: extend top age range to 100 (was 90)
# Updated: Adds Step 5e for Public Benchmarking Package
# Updated: Adds case_definition covariate (with per-site override)
#          Prediction standardized to case_definition = "general"
# Updated: Adds "Both" sex-standardized outputs (observed sex ratio)
#          - Appends Sex = "Both" to <site>_Global_PastMonth_Results.csv
#          - Adds Group = "Both" to <site>_Global_PastMonth_Slopes.csv
#          - Appends Sex = "Both" to Deployment_RefTable.csv
# ===================================================================

# --- Step 0: Load Libraries and Configure Environment ---
library(lme4)
library(dplyr)
library(data.table)
library(splines)
library(ggplot2)
library(marginaleffects)
library(foreach)
library(doParallel)
library(clubSandwich)   # CR2 robust vcov

set.seed(123)
options(contrasts = c("contr.sum", "contr.poly"))

# --- Step 1: CONFIGURATION ---
pain_sites_to_process <- c(
  "Pain","JointPain",
  "Headache", "FacialPain", "NeckShoulderPain",
  "ChestPain", "StomachAbdominalPain",
  "BackPain", "HipPain", "KneePain",
  "HandPain", "ElbowPain", "FootPain"
)

# # For quick testing (Uncomment to test single site):
# pain_sites_to_process <- c(
#   "Pain"
#   )
# case_definition settings
CASE_LEVELS <- c("general", "chronic_explicit")
PRED_CASEDEF_DEFAULT <- "general"

# Configuration for derived metrics
REP_AGES <- c(5, 18, 35, 50, 65, 80, 100)
BANDS <- list(
  `5-17`   = c(5, 17),
  `18-35`  = c(18, 35),
  `35-50`  = c(35, 50),
  `50-65`  = c(50, 65),
  `65-80`  = c(65, 80),
  `80-100` = c(80, 100)
)

# --- Step 2: Load and Prepare Data ---
data_path <- "/Users/Patty/Desktop/R/LifespanPain/Data/PainData_CS.csv"

required_columns <- c(
  "Study","Sex","Age","Global_Region", "time_window_category",
  "case_definition",
  paste0("time_window_category_", pain_sites_to_process),
  paste0("case_definition_", pain_sites_to_process),
  pain_sites_to_process
)

# Robust selection (does nothing different if all columns exist)
header_names <- names(fread(data_path, nrows = 0))
missing_cols <- setdiff(required_columns, header_names)
if (length(missing_cols) > 0) {
  message("NOTE: These requested columns were not found in the CSV and will be skipped:\n  ",
          paste(missing_cols, collapse = ", "))
}
select_cols <- intersect(required_columns, header_names)

data <- fread(
  data_path,
  select = select_cols
)

# Ensure global case_definition column exists (site-specific columns can still override)
if (!"case_definition" %in% names(data)) {
  data[, case_definition := NA_character_]
}

# Prepare factor variables
data[, `:=`(
  Study = as.factor(Study),
  Sex = factor(as.integer(Sex), levels = c(1, 2), labels = c("Male", "Female")),
  Global_Region = as.factor(trimws(as.character(Global_Region))),
  time_window_category = as.factor(trimws(as.character(time_window_category))),
  case_definition = trimws(as.character(case_definition))
)]

if (any(nchar(levels(data$time_window_category)) == 0)) {
  data$time_window_category <- droplevels(data$time_window_category)
}

# --- Step 3: Define Global Parameters and Helper Functions ---
age_all <- data$Age[is.finite(data$Age)]

BKN <- c(5, 100)

IK3 <- as.numeric(quantile(age_all, c(.33, .67), na.rm = TRUE))
IK3 <- pmin(pmax(IK3, BKN[1] + 0.01), BKN[2] - 0.01)

message(sprintf(
  "Global knots: Boundary [%.2f, %.2f]; Internal df=3: %s",
  BKN[1], BKN[2], paste(round(IK3, 2), collapse = ", ")
))

# Helper for Area Under Curve (AUC) calculation
trapz_between <- function(x, y, L, U) {
  xs <- c(L, x[x > L & x < U], U)
  ys <- approx(x, y, xout = xs, rule = 2)$y
  sum((head(ys, -1) + tail(ys, 1)) * diff(xs) / 2)
}

# Helper: choose a prediction level that exists in the fitted model
pick_level <- function(mod, var, desired) {
  lv <- levels(model.frame(mod)[[var]])
  if (length(lv) == 0) return(desired)
  if (desired %in% lv) return(desired)
  message(sprintf("NOTE: '%s' not in levels(%s). Using '%s' for prediction.",
                  desired, var, lv[1]))
  lv[1]
}

# --- Step 4: Parallel Backend Setup ---
total_cores <- parallel::detectCores()
num_cores <- max(1, min(total_cores - 1, 8))
cl <- parallel::makeCluster(num_cores)

parallel::clusterEvalQ(cl, {
  library(lme4)
  library(data.table)
  library(splines)
  library(marginaleffects)
  library(clubSandwich)
})

doParallel::registerDoParallel(cl)

clusterExport(
  cl,
  varlist = c(
    "IK3","BKN","REP_AGES","BANDS","trapz_between","data","pain_sites_to_process",
    "CASE_LEVELS","PRED_CASEDEF_DEFAULT","pick_level"
  ),
  envir = environment()
)

# --- Step 5: The Master Function for Modeling and Statistics ---
model_and_stats_ipd_replicated <- function(
    pain_site, data,
    csv_path = "/Users/Patty/Desktop/LifespanPain/Results/Global_trajectories/Stats/"
) {
  
  # Diagnostic check
  if (!pain_site %in% names(data)) {
    stop(sprintf(
      "CRITICAL ERROR: The column '%s' was not found in the loaded data. Please restart your R session and run the entire script from top to bottom.",
      pain_site
    ))
  }
  
  dir.create(csv_path, showWarnings = FALSE, recursive = TRUE)
  
  # --- 5a. Data Preparation ---
  analysis_data <- data[
    !is.na(get(pain_site)) & !is.na(Age) & !is.na(Sex) &
      !is.na(Study) & !is.na(Global_Region) & !is.na(time_window_category)
  ]
  
  analysis_data[[pain_site]] <- as.numeric(as.character(analysis_data[[pain_site]]))
  
  analysis_data <- analysis_data[, `:=`(
    Sex = droplevels(Sex),
    Global_Region = droplevels(Global_Region),
    time_window_category = droplevels(time_window_category),
    Study = droplevels(Study)
  )]
  
  analysis_data[, Sex := relevel(factor(Sex), ref = "Female")]
  
  # Per-site override for time_window_category (existing behavior)
  win_col <- paste0("time_window_category_", pain_site)
  if (win_col %in% names(analysis_data)) {
    analysis_data[, time_window_category :=
                    ifelse(!is.na(get(win_col)),
                           as.character(get(win_col)),
                           as.character(time_window_category))
    ]
  }
  
  analysis_data[, time_window_category := factor(trimws(as.character(time_window_category)))]
  analysis_data$time_window_category <- droplevels(analysis_data$time_window_category)
  if ("past_month" %in% levels(analysis_data$time_window_category)) {
    analysis_data[, time_window_category := relevel(time_window_category, ref = "past_month")]
  }
  
  # Per-site override for case_definition (new behavior, mirrors time_window_category logic)
  if (!"case_definition" %in% names(analysis_data)) {
    analysis_data[, case_definition := NA_character_]
  }
  case_col <- paste0("case_definition_", pain_site)
  analysis_data[, case_definition := trimws(as.character(case_definition))]
  
  if (case_col %in% names(analysis_data)) {
    analysis_data[, case_definition :=
                    ifelse(!is.na(get(case_col)),
                           trimws(as.character(get(case_col))),
                           trimws(as.character(case_definition)))
    ]
  }
  
  analysis_data[case_definition == "", case_definition := NA_character_]
  analysis_data[!(case_definition %in% CASE_LEVELS), case_definition := NA_character_]
  analysis_data <- analysis_data[!is.na(case_definition)]
  analysis_data[, case_definition := factor(case_definition)]
  if ("general" %in% levels(analysis_data$case_definition)) {
    analysis_data[, case_definition := relevel(case_definition, ref = "general")]
  }
  
  gA <- "Female"; gB <- "Male"
  
  # Sex weights for "Both" (observed in this pain_site analytic sample)
  sex_wts <- prop.table(table(analysis_data$Sex))
  sex_wts <- sex_wts[c(gA, gB)]
  sex_wts[is.na(sex_wts)] <- 0
  if (sum(sex_wts) == 0) sex_wts[] <- 0.5
  sex_wts <- as.numeric(sex_wts / sum(sex_wts))
  names(sex_wts) <- c(gA, gB)
  
  # --- Age weights for "overall prevalence" (age-standardized to this analytic sample) ---
  # Use 1-year bins aligned to integer ages 5:100
  age_int <- pmin(pmax(round(analysis_data$Age), BKN[1]), BKN[2])
  age_levels <- BKN[1]:BKN[2]
  
  age_tab <- table(factor(age_int, levels = age_levels))
  age_wts <- as.numeric(age_tab) / sum(age_tab)
  names(age_wts) <- as.character(age_levels)
  
  # --- 5b. Fit the IPD Model (ONCE) ---
  message(sprintf("[%s] Fitting IPD GLMM...", pain_site))
  
  form_full <- as.formula(paste0(
    pain_site, " ~ ns(Age, df = 3, Boundary.knots = BKN) * Sex + ",
    "time_window_category + case_definition + (1 | Study) + (1 | Global_Region)"
  ))
  
  fit <- glmer(
    formula = form_full,
    data = analysis_data,
    family = binomial(link = "logit"),
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
  )
  
  # CR2 vcov function for marginaleffects
  vc_cr2 <- function(mod) {
    clus <- model.frame(mod)$Study
    V <- clubSandwich::vcovCR(mod, cluster = clus, type = "CR2")
    as.matrix(V)
  }
  
  # --- 5c. Define Prediction Grid ---
  age_grid <- seq(BKN[1], BKN[2], length.out = 101)
  
  TW_PRED <- pick_level(fit, "time_window_category", "past_month")
  CASE_PRED <- pick_level(fit, "case_definition", PRED_CASEDEF_DEFAULT)
  
  # --- 5d. Generate All Statistics and Save to Replicated CSVs ---
  
  # File 1: <site>_Global_PastMonth_Results.csv
  pred_grid <- datagrid(
    model = fit,
    Age = age_grid,
    Sex = c(gA, gB),
    time_window_category = TW_PRED,
    case_definition = CASE_PRED
  )
  
  preds <- avg_predictions(
    fit,
    newdata = pred_grid,
    by = c("Age","Sex"),
    type = "response",
    re.form = NA,
    vcov = vc_cr2
  )
  
  preds_dt <- as.data.table(preds)[, .(
    Age, Sex,
    Predicted_Prob = estimate,
    SE = std.error,
    Lower_CI = conf.low,
    Upper_CI = conf.high
  )]
  
  # Add "Both" (sex-standardized over observed sex ratio)
  pred_grid_both <- pred_grid
  pred_grid_both$wts <- ifelse(pred_grid_both$Sex == gA, sex_wts[gA], sex_wts[gB])
  
  preds_both <- avg_predictions(
    fit,
    newdata = pred_grid_both,
    by = "Age",
    type = "response",
    re.form = NA,
    vcov = vc_cr2,
    wts = "wts"
  )
  
  preds_both_dt <- as.data.table(preds_both)[, .(
    Age,
    Sex = "Both",
    Predicted_Prob = estimate,
    SE = std.error,
    Lower_CI = conf.low,
    Upper_CI = conf.high
  )]
  
  preds_dt <- rbindlist(list(preds_dt, preds_both_dt), use.names = TRUE, fill = TRUE)
  setorder(preds_dt, Sex, Age)
  
  preds_dt[, Rate_Of_Change_Annual := c(NA, diff(Predicted_Prob)) / c(NA, diff(Age)), by = Sex]
  preds_dt[, Rate_Of_Change_Annual := nafill(Rate_Of_Change_Annual, type = "nocb"), by = Sex]
  preds_dt[, Rate_Of_Change_Decade := Rate_Of_Change_Annual * 10]
  
  fwrite(preds_dt, file = file.path(csv_path, paste0(pain_site, "_Global_PastMonth_Results.csv")))
  
  # ------------------------------------------------------------------
  # NEW FILE: <site>_Global_PastMonth_OverallPrevalence.csv
  # Overall prevalence (age-standardized) with robust CR2 delta-method CI
  # ------------------------------------------------------------------
  message(sprintf("[%s] Calculating: Overall prevalence (age-standardized)", pain_site))
  
  age_grid_int <- BKN[1]:BKN[2]
  
  # Female + Male overall (age-weighted within sex)
  nd_overall_sex <- datagrid(
    model = fit,
    Age = age_grid_int,
    Sex = c(gA, gB),
    time_window_category = TW_PRED,
    case_definition = CASE_PRED
  )
  nd_overall_sex$wts <- age_wts[as.character(nd_overall_sex$Age)]
  
  overall_sex <- avg_predictions(
    fit,
    newdata = nd_overall_sex,
    by = "Sex",
    type = "response",
    re.form = NA,
    vcov = vc_cr2,
    wts = "wts"
  )
  
  overall_sex_dt <- as.data.table(overall_sex)[, .(
    Sex,
    Predicted_Prob = estimate,
    SE = std.error,
    Lower_CI = conf.low,
    Upper_CI = conf.high
  )]
  
  # "Both" overall (age weights + observed sex weights)
  nd_overall_both <- nd_overall_sex
  nd_overall_both$wts <- nd_overall_both$wts * ifelse(
    nd_overall_both$Sex == gA, sex_wts[gA], sex_wts[gB]
  )
  
  overall_both <- avg_predictions(
    fit,
    newdata = nd_overall_both,
    type = "response",
    re.form = NA,
    vcov = vc_cr2,
    wts = "wts"
  )
  
  overall_both_dt <- data.table(
    Sex = "Both",
    Predicted_Prob = overall_both$estimate,
    SE = overall_both$std.error,
    Lower_CI = overall_both$conf.low,
    Upper_CI = overall_both$conf.high
  )
  
  overall_prev_dt <- rbindlist(list(overall_sex_dt, overall_both_dt), use.names = TRUE, fill = TRUE)
  
  # Add metadata so it is explicit what this number is
  overall_prev_dt[, `:=`(
    Age_Standardization = "AnalyticSample_5to100_1y",
    time_window_category = TW_PRED,
    case_definition = CASE_PRED
  )]
  
  fwrite(
    overall_prev_dt,
    file = file.path(csv_path, paste0(pain_site, "_Global_PastMonth_OverallPrevalence.csv"))
  )
  
  # File 2: <site>_Global_PastMonth_RD_Bands.csv
  message(sprintf("[%s] Calculating: RD Bands", pain_site))
  
  rd_curve <- comparisons(
    fit,
    newdata = datagrid(
      model = fit,
      Age = age_grid,
      time_window_category = TW_PRED,
      case_definition = CASE_PRED
    ),
    variables = "Sex",
    type = "response",
    comparison = "difference",
    re.form = NA,
    vcov = vc_cr2
  )
  
  rd_bands_dt <- data.table(
    Age = rd_curve$Age,
    RD = rd_curve$estimate,
    RD_L = rd_curve$conf.low,
    RD_U = rd_curve$conf.high
  )
  
  fwrite(rd_bands_dt, file = file.path(csv_path, paste0(pain_site, "_Global_PastMonth_RD_Bands.csv")))
  
  # File 3: <site>_Global_PastMonth_IntegratedDifference.csv
  message(sprintf("[%s] Calculating: Integrated Difference (AUC)", pain_site))
  
  auc_pt_prop <- trapz_between(rd_curve$Age, rd_curve$estimate, BKN[1], BKN[2])
  
  auc_dt <- data.table(
    Age_Min = BKN[1], Age_Max = BKN[2],
    AUC_RD_pp_years = auc_pt_prop * 100,
    AUC_RD_pp_years_L = NA_real_, AUC_RD_pp_years_U = NA_real_,
    Mean_RD_pp = (auc_pt_prop * 100) / (BKN[2] - BKN[1]),
    Mean_RD_pp_L = NA_real_, Mean_RD_pp_U = NA_real_,
    GroupA = gA, GroupB = gB
  )
  
  fwrite(auc_dt, file = file.path(csv_path, paste0(pain_site, "_Global_PastMonth_IntegratedDifference.csv")))
  
  # File 4: <site>_Global_PastMonth_Differences.csv
  message(sprintf("[%s] Calculating: Differences by Age Band", pain_site))
  
  diffs_list <- list()
  
  rd_overall <- avg_comparisons(
    fit,
    variables = "Sex",
    type = "response",
    re.form = NA,
    comparison = "difference",
    newdata = datagrid(
      model = fit,
      Age = age_grid,
      time_window_category = TW_PRED,
      case_definition = CASE_PRED
    ),
    vcov = vc_cr2
  )
  
  rr_overall <- avg_comparisons(
    fit,
    variables = "Sex",
    type = "response",
    re.form = NA,
    comparison = "ratio",
    newdata = datagrid(
      model = fit,
      Age = age_grid,
      time_window_category = TW_PRED,
      case_definition = CASE_PRED
    ),
    vcov = vc_cr2
  )
  
  diffs_list[["overall"]] <- data.table(
    age_band = "overall",
    RD_pp = rd_overall$estimate * 100,
    RD_pp_L = rd_overall$conf.low * 100,
    RD_pp_U = rd_overall$conf.high * 100,
    RR = rr_overall$estimate,
    RR_L = rr_overall$conf.low,
    RR_U = rr_overall$conf.high
  )
  
  for (bn in names(BANDS)) {
    br <- BANDS[[bn]]
    
    rd_band <- avg_comparisons(
      fit,
      variables = "Sex",
      type = "response",
      re.form = NA,
      comparison = "difference",
      newdata = datagrid(
        model = fit,
        Age = br,
        time_window_category = TW_PRED,
        case_definition = CASE_PRED
      ),
      vcov = vc_cr2
    )
    
    rr_band <- avg_comparisons(
      fit,
      variables = "Sex",
      type = "response",
      re.form = NA,
      comparison = "ratio",
      newdata = datagrid(
        model = fit,
        Age = br,
        time_window_category = TW_PRED,
        case_definition = CASE_PRED
      ),
      vcov = vc_cr2
    )
    
    diffs_list[[bn]] <- data.table(
      age_band = bn,
      RD_pp = rd_band$estimate * 100,
      RD_pp_L = rd_band$conf.low * 100,
      RD_pp_U = rd_band$conf.high * 100,
      RR = rr_band$estimate,
      RR_L = rr_band$conf.low,
      RR_U = rr_band$conf.high
    )
  }
  
  diffs_dt <- rbindlist(diffs_list)
  fwrite(diffs_dt, file = file.path(csv_path, paste0(pain_site, "_Global_PastMonth_Differences.csv")))
  
  # File 5: <site>_Global_PastMonth_RepAges.csv
  message(sprintf("[%s] Calculating: Differences at Representative Ages", pain_site))
  
  rd_ages <- comparisons(
    fit,
    newdata = datagrid(
      model = fit,
      Age = REP_AGES,
      time_window_category = TW_PRED,
      case_definition = CASE_PRED
    ),
    variables = "Sex",
    type = "response",
    re.form = NA,
    comparison = "difference",
    vcov = vc_cr2
  )
  
  rr_ages <- comparisons(
    fit,
    newdata = datagrid(
      model = fit,
      Age = REP_AGES,
      time_window_category = TW_PRED,
      case_definition = CASE_PRED
    ),
    variables = "Sex",
    type = "response",
    re.form = NA,
    comparison = "ratio",
    vcov = vc_cr2
  )
  
  rep_dt <- data.table(
    Age = rd_ages$Age,
    RD_pp = rd_ages$estimate * 100,
    RD_pp_L = rd_ages$conf.low * 100,
    RD_pp_U = rd_ages$conf.high * 100,
    RR = rr_ages$estimate,
    RR_L = rr_ages$conf.low,
    RR_U = rr_ages$conf.high,
    GroupA = gA, GroupB = gB
  )
  
  fwrite(rep_dt, file = file.path(csv_path, paste0(pain_site, "_Global_PastMonth_RepAges.csv")))
  
  # File 6: <site>_Global_PastMonth_Slopes.csv
  message(sprintf("[%s] Calculating: Slopes by Age Band", pain_site))
  
  slope_rows <- list()
  
  for (nm in names(BANDS)) {
    win <- BANDS[[nm]]
    
    # Explicit Sex grid so Female, Male, and Both are well-defined
    nd_slp <- datagrid(
      model = fit,
      Age = win,
      Sex = c(gA, gB),
      time_window_category = TW_PRED,
      case_definition = CASE_PRED
    )
    nd_slp$wts <- ifelse(nd_slp$Sex == gA, sex_wts[gA], sex_wts[gB])
    
    slp <- avg_slopes(
      fit,
      newdata = nd_slp,
      by = "Sex",
      type = "response",
      re.form = NA,
      variable = "Age",
      vcov = vc_cr2
    )
    
    slp_both <- avg_slopes(
      fit,
      newdata = nd_slp,
      type = "response",
      re.form = NA,
      variable = "Age",
      vcov = vc_cr2,
      wts = "wts"
    )
    
    # Female row
    slope_rows[[length(slope_rows) + 1]] <- data.table(
      Window = nm, Age_L = win[1], Age_H = win[2],
      Group = gA,
      Slope_pp_per_decade = slp[slp$Sex == gA, ]$estimate * 1000,
      L = slp[slp$Sex == gA, ]$conf.low * 1000,
      U = slp[slp$Sex == gA, ]$conf.high * 1000
    )
    
    # Male row
    slope_rows[[length(slope_rows) + 1]] <- data.table(
      Window = nm, Age_L = win[1], Age_H = win[2],
      Group = gB,
      Slope_pp_per_decade = slp[slp$Sex == gB, ]$estimate * 1000,
      L = slp[slp$Sex == gB, ]$conf.low * 1000,
      U = slp[slp$Sex == gB, ]$conf.high * 1000
    )
    
    # Both row (sex-standardized)
    slope_rows[[length(slope_rows) + 1]] <- data.table(
      Window = nm, Age_L = win[1], Age_H = win[2],
      Group = "Both",
      Slope_pp_per_decade = slp_both$estimate * 1000,
      L = slp_both$conf.low * 1000,
      U = slp_both$conf.high * 1000
    )
    
    # Female - Male row (kept for backward compatibility)
    # Note: CI here is computed by the model-based method implicit in avg_slopes only for each sex.
    # If you need a formally correct CI for the difference-in-slopes, we can compute it via a dedicated contrast step.
    f_est <- slp[slp$Sex == gA, ]$estimate
    m_est <- slp[slp$Sex == gB, ]$estimate
    
    slope_rows[[length(slope_rows) + 1]] <- data.table(
      Window = nm, Age_L = win[1], Age_H = win[2],
      Group = paste0(gA, "-", gB),
      Slope_pp_per_decade = (f_est - m_est) * 1000,
      L = NA_real_,
      U = NA_real_
    )
  }
  
  slopes_dt <- rbindlist(slope_rows)
  fwrite(slopes_dt, file = file.path(csv_path, paste0(pain_site, "_Global_PastMonth_Slopes.csv")))
  
  # File 7: <site>_Global_PastMonth_Peaks_And_Interaction.csv
  message(sprintf("[%s] Calculating: Peak Difference and Interaction", pain_site))
  
  peak_info <- rd_bands_dt[which.max(abs(RD))]
  
  excl0 <- (rd_bands_dt$RD_L > 0 & rd_bands_dt$RD_U > 0) |
    (rd_bands_dt$RD_L < 0 & rd_bands_dt$RD_U < 0)
  
  r <- rle(excl0)
  clear_start <- NA_real_
  clear_end <- NA_real_
  
  if (any(r$values)) {
    idx_run <- which.max(ifelse(r$values, r$lengths, 0))
    end_pos <- cumsum(r$lengths)[idx_run]
    start_pos <- end_pos - r$lengths[idx_run] + 1
    clear_start <- age_grid[start_pos]
    clear_end <- age_grid[end_pos]
  }
  
  # Interaction LRT: remove the spline-by-sex interaction term
  form_add <- update(form_full, . ~ . - ns(Age, df = 3, Boundary.knots = BKN):Sex)
  fit_add <- update(fit, formula. = form_add)
  
  lrt <- anova(fit_add, fit, test = "LRT")
  p_interaction <- as.numeric(lrt$`Pr(>Chisq)`[2])
  
  peaks_dt <- data.table(
    Peak_RD_pp = peak_info$RD * 100,
    Peak_RD_pp_L = peak_info$RD_L * 100,
    Peak_RD_pp_U = peak_info$RD_U * 100,
    Peak_Age = peak_info$Age,
    Peak_Age_L = NA_real_,
    Peak_Age_U = NA_real_,
    ClearIntervalStart = clear_start,
    ClearIntervalEnd = clear_end,
    Interaction_p = p_interaction,
    Boots_Successful = NA_integer_
  )
  
  fwrite(peaks_dt, file = file.path(csv_path, paste0(pain_site, "_Global_PastMonth_Peaks_And_Interaction.csv")))
  
  # --- 5e. EXPORT DEPLOYMENT PACKAGE (UNIVERSAL VERSION) ---
  
  deploy_dir <- file.path(csv_path, "Deployment_Package")
  dir.create(deploy_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1. Calculate Taus (Standard and Robust)
  re <- ranef(fit)
  grp_names <- names(re)
  study_grp_name <- grp_names[grepl("Study", grp_names)][1]
  u_vals <- re[[study_grp_name]][, 1]
  
  vc <- VarCorr(fit)
  sd_region_std <- if ("Global_Region" %in% names(vc)) attr(vc[["Global_Region"]], "stddev") else 0
  sd_study_std  <- attr(vc[[study_grp_name]], "stddev")
  tau_standard  <- sqrt(sd_region_std^2 + sd_study_std^2)
  
  mad_region <- 0
  if ("Global_Region" %in% names(re)) {
    u_reg <- re[["Global_Region"]][, 1]
    mad_region <- mad(u_reg, center = median(u_reg), constant = 1.4826)
  }
  mad_study <- mad(u_vals, center = median(u_vals), constant = 1.4826)
  tau_robust <- sqrt(mad_region^2 + mad_study^2)
  
  # 2. Generate Diagnostic Plot
  png(filename = file.path(deploy_dir, paste0(pain_site, "_Heterogeneity_Diagnostic.png")),
      width = 800, height = 600)
  
  par(mar = c(5, 4, 4, 2) + 0.1)
  hist(u_vals, breaks = 30, probability = TRUE, col = "grey90", border = "grey60",
       main = paste0("Heterogeneity Distribution: ", pain_site),
       xlab = "Study Deviation (Logits)",
       sub = "Bars = Study BLUPs | Blue = Robust (MAD) | Red = Standard (SD)")
  
  lines(density(u_vals), lwd = 2, col = "black")
  abline(v = c(-tau_robust, tau_robust), col = "blue", lwd = 3, lty = 1)
  abline(v = c(-tau_standard, tau_standard), col = "red", lwd = 2, lty = 2)
  
  legend("topright",
         legend = c(sprintf("Robust Tau (MAD): %.2f", tau_robust),
                    sprintf("Standard Tau (SD): %.2f", tau_standard)),
         col = c("blue", "red"), lty = c(1, 2), lwd = c(3, 2))
  
  dev.off()
  
  # 3. Save Params
  param_dt <- data.table(
    pain_site = pain_site,
    tau_standard = tau_standard,
    tau_robust = tau_robust,
    n_studies = length(u_vals),
    sex_weight_female = sex_wts[gA],
    sex_weight_male = sex_wts[gB],
    note = "Reference Table contains ALL combinations of Window and CaseDef. Sex includes Female, Male, Both."
  )
  fwrite(param_dt, file = file.path(deploy_dir, paste0(pain_site, "_Deployment_Params.csv")))
  
  # 4. Export UNIVERSAL Reference Table
  windows_in_model <- levels(model.frame(fit)$time_window_category)
  cases_in_model   <- levels(model.frame(fit)$case_definition)
  
  deploy_grid <- expand.grid(
    Age = seq(BKN[1], BKN[2], by = 1),
    Sex = c("Female", "Male"),
    time_window_category = windows_in_model,
    case_definition = cases_in_model,
    stringsAsFactors = FALSE
  )
  
  preds_deploy <- avg_predictions(
    fit,
    newdata = deploy_grid,
    by = c("Age", "Sex", "time_window_category", "case_definition"),
    type = "response",
    re.form = NA
  )
  
  ref_table <- as.data.table(preds_deploy)[, .(
    Age, Sex, time_window_category, case_definition,
    Expected_Prob = estimate
  )]
  
  # Add Sex = Both (sex-standardized)
  ref_both <- ref_table[Sex %in% c(gA, gB)][, .(
    Expected_Prob = sum(Expected_Prob * ifelse(Sex == gA, sex_wts[gA], sex_wts[gB]))
  ), by = .(Age, time_window_category, case_definition)]
  ref_both[, Sex := "Both"]
  
  ref_table <- rbindlist(list(ref_table, ref_both), use.names = TRUE, fill = TRUE)
  
  fwrite(ref_table, file = file.path(deploy_dir, paste0(pain_site, "_Deployment_RefTable.csv")))
  
  message(sprintf("[%s] Universal Deployment package saved (Robust Tau=%.3f).",
                  pain_site, tau_robust))
  
  message(sprintf("SUCCESS: Completed IPD analysis for %s. All files saved.", pain_site))
  invisible(NULL)
}

clusterExport(cl, varlist = c("model_and_stats_ipd_replicated"), envir = environment())

# --- Step 6: Execute Parallel Processing ---
res <- foreach(
  site = pain_sites_to_process,
  .packages = character(0),
  .inorder = FALSE
) %dopar% {
  model_and_stats_ipd_replicated(site, data)
}

# --- Step 7: Clean Up ---
stopCluster(cl)
gc()
cat("Global IPD modeling and analytical statistics completed for all pain sites.\n")

# ---------------- Setup ----------------
library(data.table)
library(ggplot2)
library(scales)

# Paths
res_dir  <- "/Users/Patty/Desktop/LifespanPain/Results/LOSO/"
fig_dir  <- "/Users/Patty/Desktop/LifespanPain/Figures/Extended/Sensitivity_Analysis/LOSO/"

# ---------------- 1. Load Data ----------------
file_path <- file.path(res_dir, "LOSO_SUMMARY_ALL_SITES.csv")
if (!file.exists(file_path)) stop("File not found: ", file_path)

all_tab <- fread(file_path)
all_tab <- all_tab[!is.na(MAE)]

# ---------------- 2. Settings ----------------
pub_width <- 7.5
base_font_size <- 10

# ---------------- 3. Define Sort Orders (Highest to Lowest) ----------------

# A. Order Pain Sites: Highest Average MAE -> Lowest Average MAE
# The minus sign (-) in reorder sorts Descending
all_tab[, Pain_Site := reorder(factor(Pain_Site), -MAE, FUN = mean)]

# B. Order Studies: Highest Average MAE (across all sites) -> Lowest
# We calculate this manually to help split the pages correctly later
study_stats <- all_tab[, .(MeanMAE = mean(MAE, na.rm = TRUE)), by = Excluded_Study]
setorder(study_stats, -MeanMAE) # Descending order (Most influential first)
ordered_studies <- study_stats$Excluded_Study

# Apply this order to the main table factor levels
# Note: ggplot plots Y-axis from bottom to top.
# To have the "Highest" study at the TOP, we reverse the levels here.
all_tab[, Excluded_Study := factor(Excluded_Study, levels = rev(ordered_studies))]

# ---------------- 4. Summary Boxplot ----------------

p_summary <- ggplot(all_tab, aes(x = Pain_Site, y = MAE)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "grey92", width = 0.6) +
  geom_jitter(aes(color = Pain_Site), width = 0.15, alpha = 0.7, size = 1.2) +
  
  scale_y_continuous(labels = percent_format(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.05))) +
  
  labs(
    title = "Sensitivity Analysis: Stability of Global Trajectories",
    subtitle = "Ordered by average sensitivity (Highest to Lowest)",
    y = "Mean Absolute Difference vs. Full Model",
    x = NULL
  ) +
  
  theme_classic(base_size = base_font_size) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
  )

ggsave(file.path(fig_dir, "LOSO_Summary_Boxplot.pdf"),
       p_summary, width = pub_width, height = 5, units = "in")


# ---------------- 5. Summary Heatmap (Split into 2 Pages) ----------------

# Split the ORDERED list of studies
n_studies <- length(ordered_studies)
mid_point <- ceiling(n_studies / 2)

# Page 1: Top half (Most Influential)
studies_p1 <- ordered_studies[1:mid_point]
# Page 2: Bottom half (Least Influential)
studies_p2 <- ordered_studies[(mid_point + 1):n_studies]

study_batches <- list("Page1" = studies_p1, "Page2" = studies_p2)

# Global color limits for consistency
mae_range <- range(all_tab$MAE, na.rm = TRUE)

for (page_name in names(study_batches)) {
  
  batch_studies <- study_batches[[page_name]]
  
  # Filter data
  plot_dat <- all_tab[Excluded_Study %in% batch_studies]
  
  # Title info
  page_desc <- ifelse(page_name == "Page1",
                      "(Most Influential Studies)",
                      "(Least Influential Studies)")
  
  p_heat <- ggplot(plot_dat, aes(x = Pain_Site, y = Excluded_Study, fill = MAE)) +
    geom_tile(color = "white", linewidth = 0.1) +
    
    scale_fill_viridis_c(option = "mako", direction = -1,
                         labels = percent_format(accuracy = 0.1),
                         name = "MAE",
                         limits = mae_range) +
    
    labs(
      title = paste("Study Influence Matrix:", page_desc),
      subtitle = "Ordered by average influence (Highest to Lowest)",
      x = NULL,
      y = "Excluded Study"
    ) +
    
    theme_minimal(base_size = base_font_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 8, color = "black"),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "grey30")
    )
  
  out_name <- sprintf("LOSO_Summary_Heatmap_%s.pdf", page_name)
  ggsave(file.path(fig_dir, out_name),
         p_heat, width = pub_width, height = 8.5, units = "in")
}

message("Done. Plots ordered Highest -> Lowest MAE.")
