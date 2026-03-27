# -------------------------------
# Age-adjusted RR for BMI, Smoking, Income (level vs reference)
# Separate GLMM per exposure × pain site; one combined CSV per exposure
# -------------------------------

library(lme4)
library(dplyr)
library(data.table)
library(splines)
library(foreach)
library(doParallel)
library(tidyr)
library(MASS)

set.seed(123)
options(contrasts = c("contr.sum", "contr.poly"))

# --- Config ---
tw_levels <- c("past_week","past_month","past_3_plus_months","chronic_explicit")
AGE_MIN <- 5
AGE_MAX <- 100
NSIMS   <- 1000

pain_sites <- c('Headache','FacialPain','NeckShoulderPain','ChestPain',
                'StomachAbdominalPain','BackPain','HipPain','KneePain',
                'HandPain', 'FootPain', 'ElbowPain','Pain')

# --- Define Exposure Variables ---
# NOTE: Smoking now uses SmokingStatus2 (collapsed to Non-Smoker vs Smoker),
# but file_suffix and output path remain unchanged.
exposure_analyses <- list(
  BMI = list(
    var_name   = "BMI_Category",
    ref_level  = "Normal",
    levels     = c("Normal", "Overweight", "Obese"),
    file_suffix= "BMI"
  ),
  Smoking = list(
    var_name   = "SmokingStatus2",                # <— changed to collapsed variable
    ref_level  = "Non-Smoker",
    levels     = c("Non-Smoker", "Smoker"),
    file_suffix= "Smoking"                        # <— unchanged (keeps same filename)
  ),
  Income = list(
    var_name   = "Income_class_quintile",
    ref_level  = "Upper",
    levels     = c("Lower", "Lower Middle", "Middle", "Upper Middle", "Upper"),
    file_suffix= "Income"
  )
)

# --- Data ---
time_window_cols <- paste0("time_window_category_", pain_sites)
required_columns <- c(
  "Study","Global_Region","Age","Sex",
  "BMI_Category","Smoking_category","SmokingStatus","Income_class_quintile",  # <— add SmokingStatus
  "time_window_category", time_window_cols, pain_sites
)

data <- fread("/Users/Patty/Desktop/R/LifespanPain/Data/PainData_CS.csv",
              select = required_columns)

# Factor setup (do not collapse levels)
data[, `:=`(
  Study         = factor(Study),
  Global_Region = factor(trimws(as.character(Global_Region))),
  Sex           = factor(Sex)
)]

# --- Harmonize smoking once, then define 2-level (never vs current); exclude ex-smokers ---
data <- data %>%
  mutate(
    Smoking_category = trimws(as.character(Smoking_category)),
    SmokingStatus    = trimws(as.character(SmokingStatus)),
    SmokingStatus = case_when(
      !is.na(SmokingStatus) ~ SmokingStatus,                # keep 3-level when present
      Smoking_category %in% c("Smoker","2","1") ~ "Smoker", # common encodings; adjust if needed
      Smoking_category %in% c("Non-Smoker","0") ~ "Non-Smoker",
      TRUE ~ NA_character_
    ),
    SmokingStatus = factor(SmokingStatus, levels = c("Non-Smoker","Ex-Smoker","Smoker")),
    # NEW: 2-level variable for this RR analysis: Ex-Smoker -> NA (drop)
    SmokingStatus2 = case_when(
      SmokingStatus == "Smoker" ~ "Smoker",
      SmokingStatus == "Non-Smoker" ~ "Non-Smoker",
      SmokingStatus == "Ex-Smoker" ~ NA_character_
    ),
    SmokingStatus2 = factor(SmokingStatus2, levels = c("Non-Smoker","Smoker"))
  )


# --- Parallel ---
total_cores <- parallel::detectCores()
num_cores   <- min(total_cores - 1, 4)
cl          <- makeCluster(num_cores)
registerDoParallel(cl)

# --- Core function (exposure-agnostic) ---
compute_rr_for_exposure <- function(pain_site, exposure, dat,
                                    out_dir = "/Users/Patty/Desktop/LifespanPain/Results/Exposure_AgeAdjusted/") {
  library(lme4); library(data.table); library(dplyr); library(splines); library(tidyr); library(MASS)
  
  exp_var   <- exposure$var_name
  exp_ref   <- exposure$ref_level
  exp_lvls  <- exposure$levels
  
  # Site-specific time window column (fallback to global)
  tw_col <- paste0("time_window_category_", pain_site)
  if (!tw_col %in% names(dat)) tw_col <- "time_window_category"
  
  # Filter usable rows & enforce exposure levels/order
  df <- as.data.table(dat)[
    !is.na(get(pain_site)) & is.finite(Age) &
      !is.na(Study) & !is.na(get(exp_var)) &
      !is.na(Sex) & !is.na(Global_Region) & !is.na(get(tw_col))
  ]
  
  df <- df[get(exp_var) %in% exp_lvls]
  df[, (exp_var) := factor(get(exp_var), levels = exp_lvls)]
  
  # Ensure binary numeric outcome
  if (!is.numeric(df[[pain_site]])) {
    df[[pain_site]] <- as.numeric(as.character(df[[pain_site]]))
  }
  
  # TW factor with present levels
  tw_fit <- tw_levels[tw_levels %in% unique(as.character(df[[tw_col]]))]
  df[, tw := factor(get(tw_col), levels = tw_fit)]
  
  # GLMM (no Age×Exposure interaction)
  form <- as.formula(paste0(
    "`", pain_site, "` ~ ns(Age, df = 2) + `", exp_var, "` + ",
    "Sex + tw + (1 | Global_Region) + (1 | Study)"
  ))
  
  fit <- tryCatch(
    glmer(form, data = df, family = binomial(link = "logit"),
          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))),
    error = function(e) { message("glmer error (", pain_site, "/", exp_var, "): ", e$message); NULL }
  )
  if (is.null(fit)) return(NULL)
  
  # Prediction grid (uniform age, equal sex, fixed TW)
  age_seq    <- seq(AGE_MIN, AGE_MAX, length.out = 150)
  sex_levels <- levels(df$Sex)
  exp_levels <- levels(df[[exp_var]])
  target_tw  <- if ("past_month" %in% levels(df$tw)) "past_month" else df[, .N, by = tw][order(-N)]$tw[1]
  
  grid <- tidyr::crossing(
    Age = age_seq,
    !!exp_var := factor(exp_levels, levels = exp_levels),
    Sex = factor(sex_levels, levels = sex_levels),
    tw  = factor(target_tw, levels = levels(df$tw))
  ) |>
    group_by(!!as.symbol(exp_var)) |>
    mutate(w = 1 / n()) |>
    ungroup() |>
    as.data.table()
  
  # Fixed-effects design for newdata
  tt   <- terms(fit, fixed.only = TRUE)
  ffix <- delete.response(tt)
  Xnew <- model.matrix(ffix, grid)
  
  # Point predictions (population-average)
  beta_hat <- lme4::fixef(fit)
  eta_hat  <- as.numeric(Xnew %*% beta_hat)
  p_hat    <- plogis(eta_hat)
  
  # Average within exposure level
  grp <- grid[[exp_var]]
  avg_point <- rowsum(p_hat * grid$w, group = grp) / rowsum(grid$w, group = grp)
  pbar_hat  <- as.numeric(avg_point[,1])
  names(pbar_hat) <- rownames(avg_point)
  
  if (!(exp_ref %in% names(pbar_hat))) return(NULL)
  
  # RR (level vs reference)
  comp_levels <- setdiff(names(pbar_hat), exp_ref)
  rr_point <- pbar_hat[comp_levels] / pbar_hat[exp_ref]
  
  # Parametric simulation over fixed effects
  Vbeta <- as.matrix(vcov(fit))
  sims  <- MASS::mvrnorm(n = NSIMS, mu = beta_hat, Sigma = Vbeta)
  
  eta_sim <- sims %*% t(Xnew)   # (NSIMS x N)
  p_sim   <- plogis(eta_sim)
  
  # Average by exposure level for each draw
  idx_list <- split(seq_len(nrow(grid)), grp)
  w_vec    <- grid$w
  avg_by_level <- function(Pmat, idx_list, w) {
    do.call(cbind, lapply(idx_list, function(idx) {
      num <- Pmat[, idx, drop = FALSE] %*% w[idx]
      den <- sum(w[idx])
      as.numeric(num / den)
    }))
  }
  pbar_draws <- avg_by_level(p_sim, idx_list, w_vec)
  colnames(pbar_draws) <- names(idx_list)
  
  # RR draws: level / reference
  rr_draws <- sapply(comp_levels, function(lvl) pbar_draws[, lvl] / pbar_draws[, exp_ref])
  
  ci_lo <- apply(rr_draws, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  ci_hi <- apply(rr_draws, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
  pval  <- apply(rr_draws, 2, function(x) { px <- mean(x <= 1, na.rm = TRUE); 2 * pmin(px, 1 - px) })
  
  out <- data.table(
    ExposureVar = exp_var,
    PainSite    = pain_site,
    Comparison  = paste0(comp_levels, " vs ", exp_ref),
    RR          = as.numeric(rr_point[comp_levels]),
    LowerCI     = as.numeric(ci_lo[comp_levels]),
    UpperCI     = as.numeric(ci_hi[comp_levels]),
    P_value     = as.numeric(pval[comp_levels]),
    TimeWindow  = as.character(target_tw),
    NSims       = NSIMS
  )
  return(out)
}

# --- Run per exposure; write one CSV per exposure ---
results_base <- "/Users/Patty/Desktop/LifespanPain/Results/Exposure_AgeAdjusted/"
dir.create(results_base, showWarnings = FALSE, recursive = TRUE)

for (exp_name in names(exposure_analyses)) {
  exp_cfg <- exposure_analyses[[exp_name]]
  
  res_list <- foreach(
    pain_site = pain_sites,
    .packages = c("lme4","dplyr","splines","data.table","tidyr","MASS")
  ) %dopar% {
    compute_rr_for_exposure(pain_site, exp_cfg, data)
  }
  
  res <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
  if (nrow(res)) {
    outfile <- file.path(results_base, paste0("AgeAdjusted_RR_", exp_cfg$file_suffix, ".csv"))
    fwrite(res, outfile)
    message("Saved: ", outfile)
  } else {
    message("No results for exposure: ", exp_name)
  }
}

stopCluster(cl)
rm(data); gc()
cat("Age-adjusted RR (BMI, Smoking, Income) completed.\n")
