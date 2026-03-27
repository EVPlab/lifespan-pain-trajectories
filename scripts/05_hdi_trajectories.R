# HDI non-linear #
# --- Step 0: Libraries & Setup ---
library(lme4)
library(dplyr)
library(data.table)
library(splines)
library(ggplot2)
library(scales)
library(forcats)
library(patchwork)
library(marginaleffects)
library(clubSandwich)  # for CR2 cluster-robust vcov

set.seed(123)
options(contrasts = c("contr.sum", "contr.poly"))

# --- Step 1: Configuration ---
pain_sites_to_process <- c("BackPain", "JointPain", "Pain", "Headache")
# pain_sites_to_process <- c("BackPain")

pain_sites_all <- pain_sites_to_process

data_path <- "/Users/Patty/Desktop/R/LifespanPain/Data/PainData_CS.csv"
results_output_path <- "/Users/Patty/Desktop/LifespanPain/Results/HDI_UN_Threshold_Analysis/"

# case_definition settings (match your global script)
CASE_LEVELS <- c("general", "chronic_explicit")
PRED_CASEDEF_DEFAULT <- "general"

# Representative HDI points for prediction (labels are for display)
hdi_points_to_predict <- c(0.55, 0.70, 0.80, 0.90)
hdi_labels <- c("Low (HDI=0.55)", "Medium (HDI=0.70)", "High (HDI=0.80)", "Very High (HDI=0.90)")
hdi_colors <- c(
  "Low (HDI=0.55)"       = "#1F4E79",
  "Medium (HDI=0.70)"    = "#3A86C8",
  "High (HDI=0.80)"      = "#72B2D7",
  "Very High (HDI=0.90)" = "#A6CEE3"
)

# --- Step 2: Load & Prepare Data ---
message("--- Loading and preparing data for HDI analysis... ---")

time_window_cols <- paste0("time_window_category_", pain_sites_all)
case_def_cols    <- paste0("case_definition_", pain_sites_all)

required_columns <- unique(c(
  "Study","Sex","Age","HDI","Global_Region",
  "time_window_category",
  "case_definition",
  time_window_cols,
  case_def_cols,
  pain_sites_all
))

if (!file.exists(data_path)) stop("Data file not found at path: ", data_path)

# Robust selection (won't error if some requested columns are missing)
header_names <- names(fread(data_path, nrows = 0))
missing_cols <- setdiff(required_columns, header_names)
if (length(missing_cols) > 0) {
  message("NOTE: These requested columns were not found in the CSV and will be skipped:\n  ",
          paste(missing_cols, collapse = ", "))
}
select_cols <- intersect(required_columns, header_names)

data <- fread(data_path, select = select_cols)
# data <- data %>% sample_frac(.75)

# Ensure global case_definition column exists even if missing in CSV
if (!"case_definition" %in% names(data)) {
  data[, case_definition := NA_character_]
}

# Factors
data[, Sex := factor(Sex, levels = c(1, 2), labels = c("Male", "Female"))]
data[, `:=`(
  Study = as.factor(Study),
  Global_Region = as.factor(trimws(as.character(Global_Region))),
  time_window_category = as.factor(trimws(as.character(time_window_category))),
  case_definition = trimws(as.character(case_definition))
)]

if (any(nchar(levels(data$time_window_category)) == 0)) {
  data$time_window_category <- forcats::fct_drop(data$time_window_category)
}

dir.create(results_output_path, showWarnings = FALSE, recursive = TRUE)

# Global spline knots
age_all <- data$Age[is.finite(data$Age)]
BKN <- range(age_all, na.rm = TRUE)
IK2 <- as.numeric(quantile(age_all, 0.50, na.rm = TRUE))
message(sprintf("Global spline knots set. Boundary: [%.2f, %.2f]; Internal(df=2): %.2f",
                BKN[1], BKN[2], IK2))

# Helper: choose a prediction level that exists in the fitted model
pick_level <- function(mod, var, desired) {
  lv <- levels(model.frame(mod)[[var]])
  if (length(lv) == 0) return(desired)
  if (desired %in% lv) return(desired)
  message(sprintf("NOTE: '%s' not in levels(%s). Using '%s' for prediction.", desired, var, lv[1]))
  lv[1]
}

# --- Step 3: Analysis Function (non-sex-stratified; spline×spline HDI; RD & RR) ---
generate_hdi_prevalence_results <- function(pain_site, full_data, results_path) {
  message(paste0("\n--- Processing Pain Site by HDI (non-sex-stratified, spline×spline): ", pain_site, " ---"))
  tryCatch({
    # 1) DATA
    analysis_data <- full_data[
      !is.na(get(pain_site)) & !is.na(Age) & !is.na(Sex) &
        !is.na(HDI) & !is.na(Study) & !is.na(Global_Region)
    ]
    if (nrow(analysis_data) == 0L) {
      warning("No data after filtering.")
      return(NULL)
    }
    
    analysis_data[, `:=`(
      Study = droplevels(Study),
      Global_Region = droplevels(Global_Region)
    )]
    
    # Per-site override for time_window_category
    win_col <- paste0("time_window_category_", pain_site)
    if (win_col %in% names(analysis_data)) {
      analysis_data[, time_window_category := ifelse(
        !is.na(get(win_col)),
        as.character(get(win_col)),
        as.character(time_window_category)
      )]
    }
    analysis_data[, time_window_category := factor(
      trimws(time_window_category),
      levels = c("past_week","past_month","past_3_plus_months","chronic_explicit")
    )]
    
    # Per-site override for case_definition
    if (!"case_definition" %in% names(analysis_data)) {
      analysis_data[, case_definition := NA_character_]
    }
    case_col <- paste0("case_definition_", pain_site)
    analysis_data[, case_definition := trimws(as.character(case_definition))]
    
    if (case_col %in% names(analysis_data)) {
      analysis_data[, case_definition := ifelse(
        !is.na(get(case_col)),
        trimws(as.character(get(case_col))),
        trimws(as.character(case_definition))
      )]
    }
    
    analysis_data[case_definition == "", case_definition := NA_character_]
    analysis_data[!(case_definition %in% CASE_LEVELS), case_definition := NA_character_]
    analysis_data <- analysis_data[!is.na(case_definition)]
    analysis_data[, case_definition := factor(case_definition)]
    if ("general" %in% levels(analysis_data$case_definition)) {
      analysis_data[, case_definition := relevel(case_definition, ref = "general")]
    }
    
    # Ensure numeric outcome
    if (!is.numeric(analysis_data[[pain_site]])) {
      analysis_data[[pain_site]] <- as.numeric(as.character(analysis_data[[pain_site]]))
    }
    
    # Report HDI range
    message(sprintf("HDI range in data: [%.3f, %.3f], N = %d",
                    min(analysis_data$HDI, na.rm = TRUE),
                    max(analysis_data$HDI, na.rm = TRUE),
                    nrow(analysis_data)))
    
    # 2) MODEL — ns(Age) × ns(HDI) for flexible surface
    message("Fitting mixed-effects model: ns(Age, df=2) * ns(HDI, df=2) + ns(Age, df=2) * Sex ...")
    model_formula <- as.formula(paste0(
      pain_site, " ~ ",
      "ns(Age, df=2) * ns(HDI, df=2) + ",
      "ns(Age, df = 2) *  Sex +",  
      "time_window_category + case_definition + ",
      "(1|Study) + (1|Global_Region:Study)"
    ))
    model <- glmer(
      model_formula,
      data = analysis_data,
      family = binomial(link = "logit"),
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
    )
    if (isSingular(model, tol = 1e-4)) {
      message(">> Warning: model is near-singular. Consider simplifying if needed.")
    }
    
    # CR2 vcov function (Study-clustered)
    vc_cr2 <- function(mod) {
      cl <- model.frame(mod)$Study
      as.matrix(clubSandwich::vcovCR(mod, cluster = cl, type = "CR2"))
    }
    
    # 3) PREDICTION GRID
    min_age <- max(min(analysis_data$Age, na.rm = TRUE), BKN[1])
    max_age <- min(max(analysis_data$Age, na.rm = TRUE), BKN[2])
    if (!is.finite(min_age) || !is.finite(max_age) || min_age >= max_age) {
      warning("Age range invalid after boundaries; skipping.")
      return(NULL)
    }
    age_grid <- seq(min_age, max_age, length.out = 100)
    
    sex_levels <- levels(analysis_data$Sex)
    sex_w <- setNames(rep(1 / length(sex_levels), length(sex_levels)), sex_levels)
    
    win_levels <- levels(analysis_data$time_window_category)
    target_window <- if ("past_month" %in% win_levels) {
      "past_month"
    } else {
      analysis_data[, .N, by = time_window_category][order(-N), as.character(time_window_category[1])]
    }
    
    CASE_PRED <- pick_level(model, "case_definition", PRED_CASEDEF_DEFAULT)
    
    # Non-sex-stratified prevalence (CR2-robust) — using raw HDI values
    grid_resp <- CJ(
      Age = age_grid,
      Sex = sex_levels,
      HDI = hdi_points_to_predict,
      time_window_category = factor(target_window, levels = win_levels),
      case_definition = factor(CASE_PRED, levels = levels(model.frame(model)$case_definition)),
      sorted = FALSE
    )
    grid_resp[, w := sex_w[Sex]]
    
    ap <- marginaleffects::avg_predictions(
      model,
      newdata = grid_resp,
      by      = c("Age", "HDI"),
      type    = "response",
      re.form = NA,
      wts     = grid_resp$w,
      vcov    = vc_cr2
    )
    res <- as.data.table(ap)
    res[, HDI_Level := factor(
      sprintf("%.2f", HDI),
      levels = sprintf("%.2f", hdi_points_to_predict),
      labels = hdi_labels
    )]
    
    prevalence_out <- res[, .(
      Pain_Site = pain_site,
      HDI_Level,
      HDI_Value = round(HDI, 3),
      Age = round(Age, 2),
      Predicted_Probability = estimate,
      Lower_CI = conf.low,
      Upper_CI = conf.high
    )]
    
    # RD & RR (Very High vs Low), CR2-robust
    grid_contrast <- CJ(
      Age = age_grid,
      Sex = sex_levels,
      HDI = c(hdi_points_to_predict[1], hdi_points_to_predict[4]),
      time_window_category = factor(target_window, levels = win_levels),
      case_definition = factor(CASE_PRED, levels = levels(model.frame(model)$case_definition)),
      sorted = FALSE
    )
    grid_contrast[, w := sex_w[Sex]]
    
    rd_cmp <- marginaleffects::comparisons(
      model,
      variables  = list(HDI = c(hdi_points_to_predict[1], hdi_points_to_predict[4])),
      newdata    = grid_contrast,
      by         = "Age",
      type       = "response",
      comparison = "difference",
      re.form    = NA,
      wts        = grid_contrast$w,
      vcov       = vc_cr2
    )
    rd_df <- as.data.table(rd_cmp)[, .(
      Age = as.numeric(Age),
      RD  = estimate,
      RD_LCL = conf.low,
      RD_UCL = conf.high
    )]
    
    rr_cmp <- marginaleffects::comparisons(
      model,
      variables  = list(HDI = c(hdi_points_to_predict[1], hdi_points_to_predict[4])),
      newdata    = grid_contrast,
      by         = "Age",
      type       = "response",
      comparison = "ratio",
      re.form    = NA,
      wts        = grid_contrast$w,
      vcov       = vc_cr2
    )
    rr_df <- as.data.table(rr_cmp)[, .(
      Age = as.numeric(Age),
      RR  = estimate,
      RR_LCL = conf.low,
      RR_UCL = conf.high
    )]
    
    # 4) SAVE
    fwrite(prevalence_out,
           file.path(results_path, paste0("Predicted_Prevalence_HDI_nonSex_", pain_site, ".csv")))
    fwrite(rd_df,
           file.path(results_path, paste0("RD_VHigh_vs_Low_", pain_site, ".csv")))
    fwrite(rr_df,
           file.path(results_path, paste0("RR_VHigh_vs_Low_", pain_site, ".csv")))
    
    # 5) PLOTS
    min_age_plot <- round(min_age/10)*10
    max_age_plot <- round(max_age/10)*10
    
    p_prev <- ggplot(prevalence_out,
                     aes(x = Age, y = Predicted_Probability,
                         color = HDI_Level, fill = HDI_Level)) +
      geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.10, linetype = 0) +
      geom_line(linewidth = 1.2) +
      scale_y_continuous(name = "Predicted Prevalence",
                         labels = percent_format(accuracy = 1),
                         expand = expansion(mult = c(0.05, 0.05))) +
      scale_x_continuous(name = NULL, breaks = seq(min_age_plot, max_age_plot, by = 10)) +
      scale_color_manual(values = hdi_colors, name = "HDI Level") +
      scale_fill_manual(values = hdi_colors, name = "HDI Level", guide = "none") +
      labs(
        title = paste("Predicted Lifespan Prevalence of", pain_site, "by HDI"),
        subtitle = "Non-sex-stratified; ns(Age,2)×ns(HDI,2); standardized to past-month recall & general case definition"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
    
    p_rd <- ggplot(rd_df, aes(x = Age, y = RD)) +
      geom_ribbon(aes(ymin = RD_LCL, ymax = RD_UCL), alpha = 0.2) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      scale_y_continuous(name = "Risk Difference (V.High − Low)",
                         labels = percent_format(accuracy = 1)) +
      scale_x_continuous(name = "Age (Years)", breaks = seq(min_age_plot, max_age_plot, by = 10)) +
      theme_minimal(base_size = 12)
    
    p_rr <- ggplot(rr_df, aes(x = Age, y = RR)) +
      geom_ribbon(aes(ymin = RR_LCL, ymax = RR_UCL), alpha = 0.2) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      scale_y_continuous(name = "Risk Ratio (V.High / Low)") +
      scale_x_continuous(name = "Age (Years)", breaks = seq(min_age_plot, max_age_plot, by = 10)) +
      theme_minimal(base_size = 12)
    
    combined <- p_prev / (p_rd | p_rr) + plot_layout(heights = c(3, 1.2))
    
    ggsave(file.path(results_path, paste0("HDI_nonSex_", pain_site, "_prev_RD_RR.pdf")),
           combined, width = 9, height = 6.5)
    
    list(plot = combined,
         prevalence = prevalence_out,
         RD = rd_df,
         RR = rr_df)
    
  }, error = function(e) {
    message(paste("<<<< ERROR processing", pain_site, "by HDI (non-sex) >>>>"))
    message(conditionMessage(e))
    NULL
  })
}

# --- Step 4: Run Analysis ---
message("\n--- Starting analysis for all specified pain sites (non-sex-stratified, spline×spline HDI) ---")
all_results <- vector("list", length(pain_sites_to_process))
names(all_results) <- pain_sites_to_process

for (site in pain_sites_to_process) {
  all_results[[site]] <- generate_hdi_prevalence_results(
    pain_site = site,
    full_data = data,
    results_path = results_output_path
  )
}

# --- Step 5: Save & Display ---
message("\n--- All HDI analyses complete. Processing and displaying results... ---")
successful <- vapply(all_results, function(x) !is.null(x), logical(1))
ok_names   <- names(all_results)[successful]

# Helper to tag and bind
bind_with_tag <- function(lst, slot) {
  rbindlist(lapply(ok_names, function(nm) {
    dt <- copy(lst[[nm]][[slot]])
    dt[, Pain_Site := nm]
    dt
  }), use.names = TRUE, fill = TRUE)
}

if (any(successful)) {
  prev_all <- bind_with_tag(all_results, "prevalence")
  rd_all   <- bind_with_tag(all_results, "RD")
  rr_all   <- bind_with_tag(all_results, "RR")
  
  fwrite(prev_all, file.path(results_output_path, "Predicted_Prevalence_ALL_SITES_HDI_nonSex.csv"))
  fwrite(rd_all,   file.path(results_output_path, "RD_VHigh_vs_Low_ALL_SITES.csv"))
  fwrite(rr_all,   file.path(results_output_path, "RR_VHigh_vs_Low_ALL_SITES.csv"))
  
  for (site_name in ok_names) {
    message(paste("\nDisplaying HDI plot (non-sex, spline×spline) for:", site_name))
    print(all_results[[site_name]]$plot)
  }
} else {
  message("No results were generated successfully.")
}

message("\n--- Script finished. ---")