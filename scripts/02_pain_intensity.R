# Load necessary libraries
library(lme4)
library(dplyr)
library(splines)
library(ggplot2)
library(merTools)
library(glmmTMB)
library(scales)
library(ggthemes)
library(gridExtra)
library(cowplot)
library(MASS)       # mvrnorm
library(data.table) # fwrite, fast ops

set.seed(123)

# Load your data
data <- read.csv("/Users/Patty/Desktop/R/LifespanPain/Data/PainIntensity_data.csv")

# Ensure 'Study' and 'Sex' are treated as factors
data$Study         <- factor(data$Study)
data$Sex           <- factor(data$Sex) # coded "1" = Male, "2" = Female in your file
data$Global_Region <- factor(data$Global_Region)

# Ensure intensity_question_window is a factor with desired levels
data$intensity_question_window <- as.character(data$intensity_question_window)
data$intensity_question_window[is.na(data$intensity_question_window)] <- "not_specified"
data$intensity_question_window <- factor(
  data$intensity_question_window,
  levels = c("past_week", "past_day", "past_month", "not_specified")
)

# Age grid and bins
AGE_GRID <- seq(5, 100, length.out = 100)
REP_AGES <- c(5, 18, 35, 50, 65, 80, 100)

SLOPE_WINDOWS <- list(
  `5-17`   = c(5, 17),
  `18-35`  = c(18, 35),
  `35-50`  = c(35, 50),
  `50-65`  = c(50, 65),
  `65-80`  = c(65, 80),
  `80-100` = c(80, 100)
)

BANDS <- list(
  `5-17`   = c(5, 17),
  `18-35`  = c(18, 35),
  `35-50`  = c(35, 50),
  `50-65`  = c(50, 65),
  `65-80`  = c(65, 80),
  `80-100` = c(80, 100)
)

AUC_RANGE <- c(5, 100)
N_SIM     <- 1000

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

trapz_between <- function(x, y, L, U) {
  xs <- c(L, x[x > L & x < U], U)
  ys <- approx(x, y, xout = xs, rule = 2)$y
  sum((head(ys, -1) + tail(ys, 1)) * diff(xs) / 2)
}

theme_Publication <- function(base_size = 6.25, base_family = "sans") {
  ggthemes::theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title      = element_blank(),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background  = element_rect(fill = "transparent", colour = NA),
      panel.border     = element_rect(colour = NA),
      axis.title       = element_blank(),
      axis.line        = element_line(color = "black", size = 0.25),
      axis.ticks       = element_line(color = "black", size = 0.25),
      axis.ticks.length = unit(1.75, "pt"),
      axis.text.x      = element_text(size = base_size, margin = margin(t = 1.5, unit = "pt")),
      axis.text.y      = element_text(size = base_size, margin = margin(r = 1.5, unit = "pt")),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "none",
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text       = element_blank()
    )
}

# ---------------------------------------------------------------------------
# DEPLOYMENT HELPER
# Mirrors Step 5e from the Global Trajectories script, adapted for glmmTMB.
#
# Arguments:
#   model         : fitted glmmTMB model
#   data_filtered : analytic sample used to fit the model
#   outcome_name  : string label used in file names (e.g. "AveragePainIntensity")
#   link          : "log" (Gamma) or "logit" (binomial)
#   fem_code      : level string for Female in data$Sex (e.g. "2")
#   mal_code      : level string for Male  in data$Sex (e.g. "1")
#   deploy_dir    : output directory for deployment package files
# ---------------------------------------------------------------------------
save_deployment_package <- function(model,
                                    data_filtered,
                                    outcome_name,
                                    link,          # "log" or "logit"
                                    fem_code,
                                    mal_code,
                                    deploy_dir) {
  
  dir.create(deploy_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---- 1. Tau (Standard and Robust) from Study random effect ---------------
  
  re_cond   <- ranef(model)$cond            # list keyed by grouping factor name
  vc_cond   <- VarCorr(model)$cond
  
  # Study random intercepts
  u_study   <- re_cond[["Study"]][, "(Intercept)"]
  sd_study  <- attr(vc_cond[["Study"]], "stddev")[[1]]
  
  # Global_Region random intercepts (if present)
  if ("Global_Region" %in% names(re_cond)) {
    u_region  <- re_cond[["Global_Region"]][, "(Intercept)"]
    sd_region <- attr(vc_cond[["Global_Region"]], "stddev")[[1]]
  } else {
    u_region  <- numeric(0)
    sd_region <- 0
  }
  
  # Standard tau (SD-based, pooled across grouping levels)
  tau_standard <- sqrt(sd_region^2 + sd_study^2)
  
  # Robust tau (MAD-based, resistant to outlier studies)
  mad_study  <- mad(u_study,  center = median(u_study),  constant = 1.4826)
  mad_region <- if (length(u_region) > 1)
    mad(u_region, center = median(u_region), constant = 1.4826)
  else 0
  tau_robust <- sqrt(mad_region^2 + mad_study^2)
  
  # ---- Residual CV (for individual-level percentile scoring) ----
  if (link == "log") {
    residual_cv <- sigma(model)
  } else {
    residual_cv <- NA_real_
  }
  
  # ---- 2. Sex weights for "Both" (observed ratio in analytic sample) -------
  
  sex_tab  <- prop.table(table(data_filtered$Sex))
  wt_fem   <- as.numeric(sex_tab[fem_code])
  wt_mal   <- as.numeric(sex_tab[mal_code])
  if (is.na(wt_fem)) wt_fem <- 0.5
  if (is.na(wt_mal)) wt_mal <- 0.5
  total    <- wt_fem + wt_mal
  wt_fem   <- wt_fem / total
  wt_mal   <- wt_mal / total
  
  # ---- 3. Heterogeneity diagnostic plot ------------------------------------
  
  png(filename = file.path(deploy_dir, paste0(outcome_name, "_Heterogeneity_Diagnostic.png")),
      width = 800, height = 600)
  
  par(mar = c(5, 4, 4, 2) + 0.1)
  hist(u_study, breaks = 30, probability = TRUE,
       col = "grey90", border = "grey60",
       main = paste0("Heterogeneity Distribution: ", outcome_name),
       xlab = paste0("Study Deviation (", link, " scale)"),
       sub  = "Bars = Study BLUPs | Blue = Robust (MAD) | Red = Standard (SD)")
  lines(density(u_study), lwd = 2, col = "black")
  abline(v = c(-tau_robust,   tau_robust),   col = "blue", lwd = 3, lty = 1)
  abline(v = c(-tau_standard, tau_standard), col = "red",  lwd = 2, lty = 2)
  legend("topright",
         legend = c(sprintf("Robust Tau (MAD): %.3f",   tau_robust),
                    sprintf("Standard Tau (SD): %.3f",  tau_standard)),
         col  = c("blue","red"), lty = c(1,2), lwd = c(3,2))
  dev.off()
  
  # ---- 4. Deployment params CSV -------------------------------------------
  
  param_dt <- data.table(
    outcome          = outcome_name,
    link             = link,
    tau_standard     = tau_standard,
    tau_robust       = tau_robust,
    residual_cv       = residual_cv,
    n_studies        = length(u_study),
    n_regions        = if (length(u_region) > 0) length(u_region) else NA_integer_,
    sex_weight_female = wt_fem,
    sex_weight_male   = wt_mal,
    note = paste0(
      "Reference table covers all observed intensity_question_window levels. ",
      "Sex = 'Both' is standardized to the observed sex ratio in the analytic sample. ",
      "Predictions are population-average (re.form = NA)."
    )
  )
  fwrite(param_dt,
         file = file.path(deploy_dir, paste0(outcome_name, "_Deployment_Params.csv")))
  
  # ---- 5. Universal reference table ----------------------------------------
  # Grid: every integer age 5-100 × all Sex levels × all window levels
  
  windows_in_model <- levels(model.frame(model)$intensity_question_window)
  sex_levels       <- levels(data_filtered$Sex)
  
  # helper: map codes -> labels for OUTPUT only
  sex_label <- function(x) {
    ifelse(x == mal_code, "Male",
           ifelse(x == fem_code, "Female",
                  ifelse(x == "Both", "Both", as.character(x))
           )
    )
  }
  
  ref_grid <- expand.grid(
    Age                      = 5:100,
    Sex                      = sex_levels,
    intensity_question_window = windows_in_model,
    stringsAsFactors = FALSE
  )
  
  # keep internal coding for prediction
  ref_grid$Sex <- factor(ref_grid$Sex, levels = levels(data_filtered$Sex))
  ref_grid$intensity_question_window <- factor(
    ref_grid$intensity_question_window,
    levels = levels(data_filtered$intensity_question_window)
  )
  
  ref_grid$Expected_Value <- predict(
    model,
    newdata  = ref_grid,
    type     = "response",
    re.form  = NA,
    allow.new.levels = TRUE
  )
  
  ref_dt <- as.data.table(ref_grid)
  
  # Add Sex = "Both" (sex-standardized over observed ratio)
  ref_both <- ref_dt[Sex %in% c(fem_code, mal_code)][
    ,
    .(Expected_Value = sum(Expected_Value * ifelse(Sex == fem_code, wt_fem, wt_mal))),
    by = .(Age, intensity_question_window)
  ]
  ref_both[, Sex := "Both"]
  
  ref_dt <- rbindlist(list(ref_dt, ref_both), use.names = TRUE, fill = TRUE)
  
  # OUTPUT: convert Sex codes to labels for the saved file
  ref_dt_out <- copy(ref_dt)
  ref_dt_out[, Sex := sex_label(as.character(Sex))]
  
  fwrite(
    ref_dt_out,
    file = file.path(deploy_dir, paste0(outcome_name, "_Deployment_RefTable.csv"))
  )
  
  message(sprintf(
    "[%s] Deployment package saved. Robust Tau = %.3f | Standard Tau = %.3f | Studies = %d",
    outcome_name, tau_robust, tau_standard, length(u_study)
  ))
  
  invisible(list(
    tau_standard  = tau_standard,
    tau_robust    = tau_robust,
    sex_wt_female = wt_fem,
    sex_wt_male   = wt_mal
  ))
}


# =========================================================================
# 1) Mean intensity (Gamma)
# =========================================================================
model_and_plot <- function(data,
                           save_path  = "~/Desktop/LifespanPain/Figures/Fig_2_Intensity_NumSites/",
                           csv_path   = "/Users/Patty/Desktop/LifespanPain/Results/PainIntensity_NumSites/",
                           deploy_dir = "/Users/Patty/Desktop/LifespanPain/Results/PainIntensity_NumSites/Deployment_Package/") {
  
  # Filter for sufferers only (1-10). 0s were set to NaN in Python.
  data_filtered <- data %>%
    filter(!is.na(AveragePainIntensity),
           AveragePainIntensity > 0,
           !is.na(Age), !is.na(Sex), !is.na(Study), !is.na(Global_Region))
  
  model_formula <- AveragePainIntensity ~ ns(Age, df = 3) * Sex +
    intensity_question_window + (1 | Study) + (1 | Global_Region)
  
  model <- glmmTMB(
    formula = model_formula,
    data    = data_filtered,
    family  = Gamma(link = "log")
  )
  
  sex_levels <- levels(data_filtered$Sex)
  predict_data <- expand.grid(
    Age = seq(5, 100, length.out = 100),
    Sex = sex_levels,
    intensity_question_window = factor("past_week",
                                       levels = levels(data_filtered$intensity_question_window))
  )
  
  pr <- predict(model, newdata = predict_data, type = "response",
                se.fit = TRUE, re.form = NA)
  predict_data$predicted_value <- pr$fit
  predict_data$lower_ci        <- pmax(0, pr$fit - 1.96 * pr$se.fit)
  predict_data$upper_ci        <-        pr$fit + 1.96 * pr$se.fit
  
  predict_data <- predict_data %>%
    group_by(Sex) %>%
    arrange(Age) %>%
    mutate(rate_of_change = c(NA, diff(predicted_value) / diff(Age))) %>%
    ungroup()
  
  peak_summary_for_plot <- predict_data %>%
    filter(!is.na(rate_of_change)) %>%
    group_by(Sex) %>%
    summarize(
      idx_max_roc          = which.max(rate_of_change),
      peak_age_display     = round(Age[idx_max_roc]),
      peak_value_at_max_roc = predicted_value[idx_max_roc],
      max_roc_rate         = rate_of_change[idx_max_roc],
      .groups = "drop"
    )
  
  if (!dir.exists(csv_path)) dir.create(csv_path, recursive = TRUE)
  predictions_df <- dplyr::select(
    predict_data, Age, Sex,
    Predicted_Value = predicted_value,
    Lower_CI = lower_ci, Upper_CI = upper_ci,
    Rate_of_Change = rate_of_change
  )
  write.csv(predictions_df,
            file = file.path(csv_path, "AveragePainIntensity_predictions_rate_of_change.csv"),
            row.names = FALSE)
  
  # ---- Manuscript stats via fixed-effects simulation ----
  fem_code <- if ("2" %in% sex_levels) "2" else "Female"
  mal_code <- if ("1" %in% sex_levels) "1" else "Male"
  
  grid <- expand.grid(
    Age = AGE_GRID, Sex = sex_levels,
    intensity_question_window = factor("past_week",
                                       levels = levels(data_filtered$intensity_question_window))
  )
  grid$Sex <- factor(grid$Sex, levels = levels(data_filtered$Sex))
  tm <- try(terms(model, component = "cond"), silent = TRUE)
  if (inherits(tm, "try-error")) tm <- terms(model)
  X <- model.matrix(stats::delete.response(tm), data = grid)
  
  beta_hat <- fixef(model)$cond
  V_beta   <- vcov(model)$cond
  draws    <- MASS::mvrnorm(n = N_SIM, mu = beta_hat, Sigma = V_beta)
  
  eta_mat <- X %*% t(draws)
  mu_mat  <- exp(eta_mat)   # Gamma(log)
  
  idx_f <- which(grid$Sex == fem_code)
  idx_m <- which(grid$Sex == mal_code)
  if (length(idx_f) == 0 || length(idx_m) == 0) stop("Sex levels not found as expected.")
  
  Af <- mu_mat[idx_f, , drop = FALSE]
  Bm <- mu_mat[idx_m, , drop = FALSE]
  
  RD_age_mat <- Af - Bm
  RR_age_mat <- Af / pmax(Bm, .Machine$double.eps)
  
  mu_pt  <- exp(X %*% beta_hat)
  A_pt   <- mu_pt[idx_f]; B_pt <- mu_pt[idx_m]
  RD_pt  <- A_pt - B_pt
  RR_pt  <- A_pt / pmax(B_pt, .Machine$double.eps)
  
  rd_overall_draw <- colMeans(Af) - colMeans(Bm)
  rr_overall_draw <- colMeans(Af) / colMeans(Bm)
  
  qfun <- function(x) quantile(x, c(.025, .975), na.rm = TRUE)
  rd_overall_ci <- qfun(rd_overall_draw)
  rr_overall_ci <- qfun(rr_overall_draw)
  
  band_rows <- list(); j <- 0L
  for (bn in names(BANDS)) {
    j <- j + 1L; br <- BANDS[[bn]]
    keep <- AGE_GRID >= br[1] & AGE_GRID <= br[2]
    rd_band_draw <- colMeans(Af[keep, , drop = FALSE]) - colMeans(Bm[keep, , drop = FALSE])
    rr_band_draw <- colMeans(Af[keep, , drop = FALSE]) / colMeans(Bm[keep, , drop = FALSE])
    band_rows[[j]] <- data.table(
      age_band    = bn,
      RD_points   = mean(RD_pt[keep]),
      RD_points_L = qfun(rd_band_draw)[1],
      RD_points_U = qfun(rd_band_draw)[2],
      RR   = mean(RR_pt[keep]),
      RR_L = qfun(rr_band_draw)[1],
      RR_U = qfun(rr_band_draw)[2]
    )
  }
  band_dt <- rbindlist(band_rows)
  diffs_dt <- rbindlist(list(
    data.table(age_band = "overall",
               RD_points   = mean(RD_pt), RD_points_L = rd_overall_ci[1],
               RD_points_U = rd_overall_ci[2],
               RR   = mean(RR_pt), RR_L = rr_overall_ci[1], RR_U = rr_overall_ci[2]),
    band_dt
  ))
  fwrite(diffs_dt, file.path(csv_path, "AveragePainIntensity_Differences.csv"))
  
  peak_rd_draw  <- apply(RD_age_mat, 2, function(v) v[which.max(abs(v))])
  peak_age_draw <- apply(RD_age_mat, 2, function(v) AGE_GRID[which.max(abs(v))])
  peak_dt <- data.table(
    Peak_RD_points   = RD_pt[which.max(abs(RD_pt))],
    Peak_RD_points_L = qfun(peak_rd_draw)[1],
    Peak_RD_points_U = qfun(peak_rd_draw)[2],
    Peak_Age         = AGE_GRID[which.max(abs(RD_pt))],
    Peak_Age_L       = qfun(peak_age_draw)[1],
    Peak_Age_U       = qfun(peak_age_draw)[2],
    Interaction_p    = NA_real_
  )
  RD_L <- apply(RD_age_mat, 1, function(x) qfun(x)[1])
  RD_U <- apply(RD_age_mat, 1, function(x) qfun(x)[2])
  RD_M <- rowMeans(RD_age_mat)
  excl0 <- (RD_L > 0 & RD_U > 0) | (RD_L < 0 & RD_U < 0)
  r <- rle(excl0)
  if (any(r$values)) {
    idx_run   <- which.max(ifelse(r$values, r$lengths, 0))
    end_pos   <- cumsum(r$lengths)[idx_run]
    start_pos <- end_pos - r$lengths[idx_run] + 1
    clear_start <- AGE_GRID[start_pos]; clear_end <- AGE_GRID[end_pos]
  } else { clear_start <- NA_real_; clear_end <- NA_real_ }
  peak_dt[, `:=`(ClearIntervalStart = clear_start, ClearIntervalEnd = clear_end)]
  fwrite(peak_dt,  file.path(csv_path, "AveragePainIntensity_Peaks_And_ClearInterval.csv"))
  fwrite(data.table(Age = AGE_GRID, RD = RD_M, RD_L = RD_L, RD_U = RD_U),
         file.path(csv_path, "AveragePainIntensity_RD_Bands.csv"))
  
  rep_rows <- lapply(REP_AGES, function(a) {
    Af_a <- apply(Af, 2, function(v) approx(AGE_GRID, v, xout = a, rule = 2)$y)
    Bm_a <- apply(Bm, 2, function(v) approx(AGE_GRID, v, xout = a, rule = 2)$y)
    RD_a <- Af_a - Bm_a
    RR_a <- Af_a / pmax(Bm_a, .Machine$double.eps)
    A_pt_a <- approx(AGE_GRID, A_pt, xout = a, rule = 2)$y
    B_pt_a <- approx(AGE_GRID, B_pt, xout = a, rule = 2)$y
    data.table(Age = a,
               RD_points   = A_pt_a - B_pt_a,
               RD_points_L = qfun(RD_a)[1], RD_points_U = qfun(RD_a)[2],
               RR   = A_pt_a / B_pt_a,
               RR_L = qfun(RR_a)[1], RR_U = qfun(RR_a)[2],
               GroupA = "Female", GroupB = "Male")
  })
  fwrite(rbindlist(rep_rows), file.path(csv_path, "AveragePainIntensity_RepAges.csv"))
  
  slope_rows <- list()
  for (nm in names(SLOPE_WINDOWS)) {
    win <- SLOPE_WINDOWS[[nm]]; L <- win[1]; H <- win[2]
    A_L <- approx(AGE_GRID, A_pt, xout = L, rule = 2)$y
    A_H <- approx(AGE_GRID, A_pt, xout = H, rule = 2)$y
    B_L <- approx(AGE_GRID, B_pt, xout = L, rule = 2)$y
    B_H <- approx(AGE_GRID, B_pt, xout = H, rule = 2)$y
    sA_pt <- (A_H - A_L)/(H-L)*10; sB_pt <- (B_H - B_L)/(H-L)*10; sD_pt <- sA_pt - sB_pt
    sA_bt <- sB_bt <- sD_bt <- numeric(N_SIM)
    for (i in seq_len(N_SIM)) {
      A_Li <- approx(AGE_GRID, Af[,i], xout = L, rule = 2)$y
      A_Hi <- approx(AGE_GRID, Af[,i], xout = H, rule = 2)$y
      B_Li <- approx(AGE_GRID, Bm[,i], xout = L, rule = 2)$y
      B_Hi <- approx(AGE_GRID, Bm[,i], xout = H, rule = 2)$y
      sA_bt[i] <- (A_Hi - A_Li)/(H-L)*10
      sB_bt[i] <- (B_Hi - B_Li)/(H-L)*10
      sD_bt[i] <- sA_bt[i] - sB_bt[i]
    }
    sA_ci <- qfun(sA_bt); sB_ci <- qfun(sB_bt); sD_ci <- qfun(sD_bt)
    slope_rows[[length(slope_rows)+1]] <- data.table(Window = nm, Age_L = L, Age_H = H,
                                                     Group = "Female",      Slope_points_per_decade = sA_pt, L = sA_ci[1], U = sA_ci[2])
    slope_rows[[length(slope_rows)+1]] <- data.table(Window = nm, Age_L = L, Age_H = H,
                                                     Group = "Male",        Slope_points_per_decade = sB_pt, L = sB_ci[1], U = sB_ci[2])
    slope_rows[[length(slope_rows)+1]] <- data.table(Window = nm, Age_L = L, Age_H = H,
                                                     Group = "Female-Male", Slope_points_per_decade = sD_pt, L = sD_ci[1], U = sD_ci[2])
  }
  fwrite(rbindlist(slope_rows), file.path(csv_path, "AveragePainIntensity_Slopes.csv"))
  
  L <- AUC_RANGE[1]; U <- AUC_RANGE[2]
  auc_pt <- trapz_between(AGE_GRID, RD_pt, L, U)
  auc_bt <- apply(RD_age_mat, 2, function(v) trapz_between(AGE_GRID, v, L, U))
  auc_ci <- qfun(auc_bt)
  fwrite(data.table(Age_Min = L, Age_Max = U,
                    AUC_RD_points_years   = auc_pt,
                    AUC_RD_points_years_L = auc_ci[1],
                    AUC_RD_points_years_U = auc_ci[2],
                    Mean_RD_points   = auc_pt/(U-L),
                    Mean_RD_points_L = auc_ci[1]/(U-L),
                    Mean_RD_points_U = auc_ci[2]/(U-L),
                    GroupA = "Female", GroupB = "Male"),
         file.path(csv_path, "AveragePainIntensity_IntegratedDifference.csv"))
  
  # Plots
  sex_colors <- c("1" = "#3AADBB", "2" = "#FB7D25")
  plot_intensity_trajectory <- ggplot(predict_data, aes(x = Age, y = predicted_value, color = Sex)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Sex), alpha = 0.1, linetype = 0) +
    geom_line(size = 1) +
    scale_x_continuous(breaks = seq(10, 100, by = 10), expand = c(0, 0)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(xlim = c(5, 100)) +
    scale_color_manual(values = sex_colors) +
    scale_fill_manual(values = sex_colors) +
    theme_Publication(base_size = 6.25) +
    geom_point(data = peak_summary_for_plot,
               aes(x = peak_age_display, y = peak_value_at_max_roc, color = Sex),
               size = 1.5, shape = 21, fill = "white", stroke = 1) +
    geom_text(data = peak_summary_for_plot,
              aes(x = peak_age_display, y = peak_value_at_max_roc,
                  label = peak_age_display, color = Sex),
              vjust = -1, size = 2.5)
  
  plot_intensity_roc <- ggplot(predict_data, aes(x = Age, y = rate_of_change, color = Sex)) +
    geom_segment(data = peak_summary_for_plot,
                 aes(x = peak_age_display, xend = peak_age_display,
                     y = 0, yend = max_roc_rate, color = Sex),
                 linetype = "11", size = 0.5) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray", size = 0.5) +
    geom_line(size = 0.5) +
    scale_x_continuous(breaks = seq(10, 100, by = 10), expand = c(0, 0)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
    coord_cartesian(xlim = c(5, 100)) +
    scale_color_manual(values = sex_colors) +
    theme_Publication(base_size = 6.25) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank(), axis.line.x  = element_blank())
  
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
  combined_plot <- plot_grid(plot_intensity_trajectory, plot_intensity_roc,
                             ncol = 1, align = "v", axis = "lr",
                             rel_heights = c(3, 1.15))
  ggsave(filename = file.path(save_path, "AveragePainIntensity_combined_plot.pdf"),
         plot = combined_plot, width = 2, height = 1.35, units = "in")
  
  # ---- DEPLOYMENT PACKAGE -------------------------------------------------
  save_deployment_package(
    model         = model,
    data_filtered = data_filtered,
    outcome_name  = "AveragePainIntensity",
    link          = "log",
    fem_code      = fem_code,
    mal_code      = mal_code,
    deploy_dir    = deploy_dir
  )
  
  cat("Saved combined plot for Average Pain Intensity and stats CSVs.\n")
  list(fit_model     = model,
       combined_plot = combined_plot,
       predictions_df = predictions_df)
}


# =========================================================================
# 2) Severe pain probability (>=7)
# =========================================================================
model_and_plot_severe <- function(data,
                                  save_path  = "~/Desktop/LifespanPain/Figures/Fig_2_Intensity_NumSites/",
                                  csv_path   = "/Users/Patty/Desktop/LifespanPain/Results/PainIntensity_NumSites/",
                                  deploy_dir = "/Users/Patty/Desktop/LifespanPain/Results/PainIntensity_NumSites/Deployment_Package/") {
  
  data_filtered <- data %>%
    filter(!is.na(SeverePain),
           !is.na(Age), !is.na(Sex), !is.na(Study), !is.na(Global_Region))
  
  model_formula <- SeverePain ~ ns(Age, df = 3) * Sex +
    intensity_question_window + (1 | Study) + (1 | Global_Region)
  
  model <- glmmTMB(
    formula = model_formula,
    data    = data_filtered,
    family  = binomial(link = "logit")
  )
  
  sex_levels <- levels(data_filtered$Sex)
  predict_data <- expand.grid(
    Age = seq(5, 100, length.out = 100),
    Sex = sex_levels,
    intensity_question_window = factor("past_week",
                                       levels = levels(data_filtered$intensity_question_window))
  )
  
  pr <- predict(model, newdata = predict_data, type = "link",
                se.fit = TRUE, re.form = NA)
  predict_data$prob      <- plogis(pr$fit)
  predict_data$lower_ci  <- plogis(pr$fit - 1.96 * pr$se.fit)
  predict_data$upper_ci  <- plogis(pr$fit + 1.96 * pr$se.fit)
  
  predict_data <- predict_data %>%
    group_by(Sex) %>%
    arrange(Age) %>%
    mutate(rate_of_change = c(NA_real_, diff(prob) / diff(Age))) %>%
    ungroup()
  
  peak_summary_for_plot <- predict_data %>%
    filter(!is.na(rate_of_change)) %>%
    group_by(Sex) %>%
    summarize(
      idx_max_roc           = which.max(rate_of_change),
      peak_age_display      = round(Age[idx_max_roc]),
      peak_prob_at_max_roc  = prob[idx_max_roc],
      max_roc_rate          = rate_of_change[idx_max_roc],
      .groups = "drop"
    )
  
  if (!dir.exists(csv_path)) dir.create(csv_path, recursive = TRUE)
  predictions_df <- dplyr::select(
    predict_data, Age, Sex,
    Predicted_Prob = prob,
    Lower_CI = lower_ci, Upper_CI = upper_ci,
    Rate_of_Change = rate_of_change
  )
  write.csv(predictions_df,
            file = file.path(csv_path, "SeverePain_gt7_predictions_rate_of_change.csv"),
            row.names = FALSE)
  
  # ---- Manuscript stats via fixed-effects simulation ----
  fem_code <- if ("2" %in% sex_levels) "2" else "Female"
  mal_code <- if ("1" %in% sex_levels) "1" else "Male"
  
  grid <- expand.grid(
    Age = AGE_GRID, Sex = sex_levels,
    intensity_question_window = factor("past_week",
                                       levels = levels(data_filtered$intensity_question_window))
  )
  grid$Sex <- factor(grid$Sex, levels = levels(data_filtered$Sex))
  tm <- try(terms(model, component = "cond"), silent = TRUE)
  if (inherits(tm, "try-error")) tm <- terms(model)
  X <- model.matrix(stats::delete.response(tm), data = grid)
  
  beta_hat <- fixef(model)$cond
  V_beta   <- vcov(model)$cond
  draws    <- MASS::mvrnorm(n = N_SIM, mu = beta_hat, Sigma = V_beta)
  
  eta_mat <- X %*% t(draws)
  mu_mat  <- plogis(eta_mat)   # Binomial(logit)
  
  idx_f <- which(grid$Sex == fem_code)
  idx_m <- which(grid$Sex == mal_code)
  Af <- mu_mat[idx_f, , drop = FALSE]
  Bm <- mu_mat[idx_m, , drop = FALSE]
  
  RD_age_mat <- Af - Bm
  RR_age_mat <- Af / pmax(Bm, .Machine$double.eps)
  
  mu_pt <- plogis(X %*% beta_hat)
  A_pt  <- mu_pt[idx_f]; B_pt <- mu_pt[idx_m]
  RD_pt <- A_pt - B_pt
  RR_pt <- A_pt / pmax(B_pt, .Machine$double.eps)
  
  rd_overall_draw <- colMeans(Af) - colMeans(Bm)
  rr_overall_draw <- colMeans(Af) / colMeans(Bm)
  qfun <- function(x) quantile(x, c(.025, .975), na.rm = TRUE)
  rd_overall_ci <- qfun(rd_overall_draw) * 100
  rr_overall_ci <- qfun(rr_overall_draw)
  
  band_rows <- list(); j <- 0L
  for (bn in names(BANDS)) {
    j <- j + 1L; br <- BANDS[[bn]]
    keep <- AGE_GRID >= br[1] & AGE_GRID <= br[2]
    rd_band_draw <- (colMeans(Af[keep, , drop = FALSE]) - colMeans(Bm[keep, , drop = FALSE])) * 100
    rr_band_draw <-  colMeans(Af[keep, , drop = FALSE]) / colMeans(Bm[keep, , drop = FALSE])
    band_rows[[j]] <- data.table(
      age_band = bn,
      RD_pp    = mean(RD_pt[keep]) * 100,
      RD_pp_L  = qfun(rd_band_draw)[1],
      RD_pp_U  = qfun(rd_band_draw)[2],
      RR   = mean(RR_pt[keep]),
      RR_L = qfun(rr_band_draw)[1],
      RR_U = qfun(rr_band_draw)[2]
    )
  }
  diffs_dt <- rbindlist(list(
    data.table(age_band = "overall",
               RD_pp   = mean(RD_pt) * 100,
               RD_pp_L = rd_overall_ci[1], RD_pp_U = rd_overall_ci[2],
               RR      = mean(RR_pt),
               RR_L    = rr_overall_ci[1], RR_U    = rr_overall_ci[2]),
    rbindlist(band_rows)
  ))
  fwrite(diffs_dt, file.path(csv_path, "SeverePain_gt7_Differences.csv"))
  
  peak_rd_draw  <- apply(RD_age_mat, 2, function(v) v[which.max(abs(v))])
  peak_age_draw <- apply(RD_age_mat, 2, function(v) AGE_GRID[which.max(abs(v))])
  RD_L <- apply(RD_age_mat, 1, function(x) qfun(x)[1])
  RD_U <- apply(RD_age_mat, 1, function(x) qfun(x)[2])
  RD_M <- rowMeans(RD_age_mat)
  excl0 <- (RD_L > 0 & RD_U > 0) | (RD_L < 0 & RD_U < 0)
  r <- rle(excl0)
  if (any(r$values)) {
    idx_run   <- which.max(ifelse(r$values, r$lengths, 0))
    end_pos   <- cumsum(r$lengths)[idx_run]
    start_pos <- end_pos - r$lengths[idx_run] + 1
    clear_start <- AGE_GRID[start_pos]; clear_end <- AGE_GRID[end_pos]
  } else { clear_start <- NA_real_; clear_end <- NA_real_ }
  
  model_add     <- update(model, . ~ . - ns(Age, df = 3):Sex)
  LRT           <- anova(model_add, model)
  p_interaction <- LRT$`Pr(>Chisq)`[2]
  
  peak_dt <- data.table(
    Peak_RD_pp   = RD_pt[which.max(abs(RD_pt))] * 100,
    Peak_RD_pp_L = qfun(peak_rd_draw * 100)[1],
    Peak_RD_pp_U = qfun(peak_rd_draw * 100)[2],
    Peak_Age     = AGE_GRID[which.max(abs(RD_pt))],
    Peak_Age_L   = qfun(peak_age_draw)[1],
    Peak_Age_U   = qfun(peak_age_draw)[2],
    ClearIntervalStart = clear_start,
    ClearIntervalEnd   = clear_end,
    Interaction_p      = p_interaction
  )
  fwrite(peak_dt, file.path(csv_path, "SeverePain_gt7_Peaks_And_Interaction.csv"))
  fwrite(data.table(Age = AGE_GRID, RD = RD_M*100, RD_L = RD_L*100, RD_U = RD_U*100),
         file.path(csv_path, "SeverePain_gt7_RD_Bands.csv"))
  
  rep_rows <- lapply(REP_AGES, function(a) {
    Af_a <- apply(Af, 2, function(v) approx(AGE_GRID, v, xout = a, rule = 2)$y)
    Bm_a <- apply(Bm, 2, function(v) approx(AGE_GRID, v, xout = a, rule = 2)$y)
    RD_a <- (Af_a - Bm_a) * 100
    RR_a <- Af_a / pmax(Bm_a, .Machine$double.eps)
    A_pt_a <- approx(AGE_GRID, A_pt, xout = a, rule = 2)$y
    B_pt_a <- approx(AGE_GRID, B_pt, xout = a, rule = 2)$y
    data.table(Age = a,
               RD_pp   = (A_pt_a - B_pt_a) * 100,
               RD_pp_L = qfun(RD_a)[1], RD_pp_U = qfun(RD_a)[2],
               RR   = A_pt_a / B_pt_a,
               RR_L = qfun(RR_a)[1], RR_U = qfun(RR_a)[2],
               GroupA = "Female", GroupB = "Male")
  })
  fwrite(rbindlist(rep_rows), file.path(csv_path, "SeverePain_gt7_RepAges.csv"))
  
  slope_rows <- list()
  for (nm in names(SLOPE_WINDOWS)) {
    win <- SLOPE_WINDOWS[[nm]]; L <- win[1]; H <- win[2]
    A_L <- approx(AGE_GRID, A_pt, xout = L, rule = 2)$y
    A_H <- approx(AGE_GRID, A_pt, xout = H, rule = 2)$y
    B_L <- approx(AGE_GRID, B_pt, xout = L, rule = 2)$y
    B_H <- approx(AGE_GRID, B_pt, xout = H, rule = 2)$y
    sA_pt <- (A_H - A_L)/(H-L)*10*100; sB_pt <- (B_H - B_L)/(H-L)*10*100; sD_pt <- sA_pt - sB_pt
    sA_bt <- sB_bt <- sD_bt <- numeric(N_SIM)
    for (i in seq_len(N_SIM)) {
      A_Li <- approx(AGE_GRID, Af[,i], xout = L, rule = 2)$y
      A_Hi <- approx(AGE_GRID, Af[,i], xout = H, rule = 2)$y
      B_Li <- approx(AGE_GRID, Bm[,i], xout = L, rule = 2)$y
      B_Hi <- approx(AGE_GRID, Bm[,i], xout = H, rule = 2)$y
      sA_bt[i] <- (A_Hi - A_Li)/(H-L)*10*100
      sB_bt[i] <- (B_Hi - B_Li)/(H-L)*10*100
      sD_bt[i] <- sA_bt[i] - sB_bt[i]
    }
    sA_ci <- qfun(sA_bt); sB_ci <- qfun(sB_bt); sD_ci <- qfun(sD_bt)
    slope_rows[[length(slope_rows)+1]] <- data.table(Window = nm, Age_L = L, Age_H = H,
                                                     Group = "Female",      Slope_pp_per_decade = sA_pt, L = sA_ci[1], U = sA_ci[2])
    slope_rows[[length(slope_rows)+1]] <- data.table(Window = nm, Age_L = L, Age_H = H,
                                                     Group = "Male",        Slope_pp_per_decade = sB_pt, L = sB_ci[1], U = sB_ci[2])
    slope_rows[[length(slope_rows)+1]] <- data.table(Window = nm, Age_L = L, Age_H = H,
                                                     Group = "Female-Male", Slope_pp_per_decade = sD_pt, L = sD_ci[1], U = sD_ci[2])
  }
  fwrite(rbindlist(slope_rows), file.path(csv_path, "SeverePain_gt7_Slopes.csv"))
  
  L <- AUC_RANGE[1]; U <- AUC_RANGE[2]
  auc_pt <- trapz_between(AGE_GRID, RD_pt * 100, L, U)
  auc_bt <- apply(RD_age_mat * 100, 2, function(v) trapz_between(AGE_GRID, v, L, U))
  auc_ci <- qfun(auc_bt)
  fwrite(data.table(Age_Min = L, Age_Max = U,
                    AUC_RD_pp_years   = auc_pt,
                    AUC_RD_pp_years_L = auc_ci[1],
                    AUC_RD_pp_years_U = auc_ci[2],
                    Mean_RD_pp   = auc_pt/(U-L),
                    Mean_RD_pp_L = auc_ci[1]/(U-L),
                    Mean_RD_pp_U = auc_ci[2]/(U-L),
                    GroupA = "Female", GroupB = "Male"),
         file.path(csv_path, "SeverePain_gt7_IntegratedDifference.csv"))
  
  # Plots
  sex_colors <- c("1" = "#3AADBB", "2" = "#FB7D25")
  plot_severe_traj <- ggplot(predict_data, aes(x = Age, y = prob, color = Sex)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Sex), alpha = 0.1, linetype = 0) +
    geom_line(size = 1) +
    scale_x_continuous(breaks = seq(10, 100, by = 10), expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    coord_cartesian(xlim = c(5, 100), ylim = c(0, NA)) +
    scale_color_manual(values = sex_colors) +
    scale_fill_manual(values = sex_colors) +
    theme_Publication(base_size = 6.25) +
    geom_point(data = peak_summary_for_plot,
               aes(x = peak_age_display, y = peak_prob_at_max_roc, color = Sex),
               size = 1.5, shape = 21, fill = "white", stroke = 1) +
    geom_text(data = peak_summary_for_plot,
              aes(x = peak_age_display, y = peak_prob_at_max_roc,
                  label = peak_age_display, color = Sex),
              vjust = -1, size = 2.5)
  
  plot_severe_roc <- ggplot(predict_data, aes(x = Age, y = rate_of_change, color = Sex)) +
    geom_segment(data = peak_summary_for_plot,
                 aes(x = peak_age_display, xend = peak_age_display,
                     y = 0, yend = max_roc_rate, color = Sex),
                 linetype = "11", size = 0.5) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray", size = 0.5) +
    geom_line(size = 0.5) +
    scale_x_continuous(breaks = seq(10, 100, by = 10), expand = c(0, 0)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
    coord_cartesian(xlim = c(5, 100)) +
    scale_color_manual(values = sex_colors) +
    theme_Publication(base_size = 6.25) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank(), axis.line.x  = element_blank())
  
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
  combined_plot <- plot_grid(plot_severe_traj, plot_severe_roc,
                             ncol = 1, align = "v", axis = "lr",
                             rel_heights = c(3, 1.15))
  ggsave(filename = file.path(save_path, "SeverePain_gt7_combined_plot.pdf"),
         plot = combined_plot, width = 2, height = 1.35, units = "in")
  
  # ---- DEPLOYMENT PACKAGE -------------------------------------------------
  save_deployment_package(
    model         = model,
    data_filtered = data_filtered,
    outcome_name  = "SeverePain_gt7",
    link          = "logit",
    fem_code      = fem_code,
    mal_code      = mal_code,
    deploy_dir    = deploy_dir
  )
  
  cat("Saved combined plot for Severe Pain and stats CSVs.\n")
  list(fit_model      = model,
       combined_plot  = combined_plot,
       predictions_df = predictions_df)
}


# =========================================================================
# Run both
# =========================================================================
results_intensity <- model_and_plot(data)
results_severe    <- model_and_plot_severe(data)