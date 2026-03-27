# ===================================================================
# Lifespan Pain Trajectories - IPD Model for GeneralizedPain_ACR2016_STRICT (single outcome)
# Updated: Adds case_definition covariate
#          Prediction standardized to case_definition = "general" (if present)
# (CR2 removed)
# ===================================================================

# --- Step 0: Libraries & Setup ---
library(lme4)
library(dplyr)
library(data.table)
library(splines)
library(ggplot2)
library(marginaleffects)

set.seed(123)
options(contrasts = c("contr.sum", "contr.poly"))

# --- Step 1: CONFIGURATION ---
outcome_name <- "GeneralizedPain_ACR2016_STRICT"   # expects 0/1/NA
data_path <- "/Users/Patty/Desktop/R/LifespanPain/Data/WidespreadPainData.csv"
results_path <- "/Users/Patty/Desktop/LifespanPain/Results/PainIntensity_NumSites/"

# case_definition settings (match your multi-site script)
CASE_LEVELS <- c("general", "chronic_explicit")
PRED_CASEDEF_DEFAULT <- "general"

REP_AGES <- c(5, 18, 35, 50, 65, 80, 100)
BANDS <- list(`5-17`=c(5,17), `18-35`=c(18,35), `35-50`=c(35,50),
              `50-65`=c(50,65), `65-80`=c(65,80), `80-100`=c(80,100))

dir.create(results_path, showWarnings = FALSE, recursive = TRUE)

# --- Step 2: Load & Prepare Data (blanket time_window_category) ---
req_cols <- c("Study","Sex","Age","Global_Region","time_window_category","case_definition", outcome_name)
stopifnot(file.exists(data_path))
data <- data.table::fread(data_path, select = req_cols)

data[, `:=`(
  Study = factor(Study),
  Sex = factor(as.integer(Sex), levels = c(1,2), labels = c("Male","Female")),
  Global_Region = factor(trimws(as.character(Global_Region))),
  time_window_category = factor(trimws(as.character(time_window_category))),
  case_definition = trimws(as.character(case_definition))
)]

# data <- data %>% sample_frac(.75)

# drop empty level if present
if (any(nchar(levels(data$time_window_category)) == 0)) {
  data$time_window_category <- droplevels(data$time_window_category)
}

# optional study exclusion (kept off by default)
# data <- data %>% filter(Study != "MCSS")

# Age knots (global)
BKN <- c(5, 100)
age_all <- data$Age[is.finite(data$Age)]
IK3 <- as.numeric(quantile(age_all, c(.33,.67), na.rm = TRUE))
IK3 <- pmin(pmax(IK3, BKN[1] + 0.01), BKN[2] - 0.01)
message(sprintf("Global knots: Boundary [%.2f, %.2f]; Internal df=3: %s",
                BKN[1], BKN[2], paste(round(IK3, 2), collapse = ", ")))

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
  message(sprintf("NOTE: '%s' not in levels(%s). Using '%s' for prediction.", desired, var, lv[1]))
  lv[1]
}

# --- Step 3: Prepare analysis data (GeneralizedPain_ACR2016_STRICT only) ---
stopifnot(outcome_name %in% names(data))

dt <- data[
  !is.na(get(outcome_name)) & !is.na(Age) & !is.na(Sex) &
    !is.na(Study) & !is.na(Global_Region) & !is.na(time_window_category) &
    !is.na(case_definition)
]

dt[[outcome_name]] <- as.numeric(as.character(dt[[outcome_name]]))

# Clean and restrict case_definition to allowed levels, mirror your other script
dt[, case_definition := trimws(as.character(case_definition))]
dt[case_definition == "", case_definition := NA_character_]
dt[!(case_definition %in% CASE_LEVELS), case_definition := NA_character_]
dt <- dt[!is.na(case_definition)]
dt[, case_definition := factor(case_definition)]
if ("general" %in% levels(dt$case_definition)) {
  dt[, case_definition := relevel(case_definition, ref = "general")]
}

dt <- dt[, `:=`(
  Sex = droplevels(Sex),
  Global_Region = droplevels(Global_Region),
  time_window_category = droplevels(time_window_category),
  case_definition = droplevels(case_definition),
  Study = droplevels(Study)
)]

dt[, Sex := relevel(factor(Sex), ref = "Female")]  # Female reference

# Relevel past_week if present
if ("past_week" %in% levels(dt$time_window_category)) {
  dt[, time_window_category := relevel(time_window_category, ref = "past_week")]
}

# --- Step 4: Fit IPD GLMM ---
form_full <- as.formula(paste0(
  outcome_name, " ~ ns(Age, df=3) * Sex + ",
  "time_window_category + case_definition + ",
  "(1 | Study) + (1 | Global_Region)"
))

fit <- glmer(
  formula = form_full, data = dt, family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))
)

# --- Step 5: Prediction grid (past_week, both sexes), standardized to case_definition ---
age_grid <- seq(BKN[1], BKN[2], length.out = 101)
gA <- "Female"; gB <- "Male"

CASE_PRED <- pick_level(fit, "case_definition", PRED_CASEDEF_DEFAULT)

pred_grid <- datagrid(
  model = fit,
  Age = age_grid,
  Sex = c(gA, gB),
  time_window_category = "past_week",
  case_definition = CASE_PRED
)

# --- File 1: Main predictions for plotting ---
preds <- avg_predictions(
  fit,
  newdata = pred_grid,
  by      = c("Age","Sex"),
  type    = "response",
  re.form = NA
)

preds_dt <- as.data.table(preds)[, .(
  Age, Sex,
  Predicted_Prob = estimate,
  SE = std.error,
  Lower_CI = conf.low,
  Upper_CI = conf.high
)]

preds_dt[, Rate_Of_Change_Annual := c(NA, diff(Predicted_Prob)) / c(NA, diff(Age)), by = Sex]
preds_dt[, Rate_Of_Change_Annual := nafill(Rate_Of_Change_Annual, type = "nocb"), by = Sex]
preds_dt[, Rate_Of_Change_Decade := Rate_Of_Change_Annual * 10]

fwrite(preds_dt, file = file.path(results_path, "GeneralizedPain_ACR2016_STRICT_Global_PastWeek_Results.csv"))

# --- File 2: RD bands (Male minus Female) over age ---
rd_curve <- comparisons(
  fit,
  newdata    = datagrid(
    model = fit,
    Age = age_grid,
    time_window_category = "past_week",
    case_definition = CASE_PRED
  ),
  variables  = "Sex",
  type       = "response",
  comparison = "difference",
  re.form    = NA
)

rd_bands_dt <- data.table(
  Age = rd_curve$Age,
  RD  = rd_curve$estimate,
  RD_L = rd_curve$conf.low,
  RD_U = rd_curve$conf.high
)

fwrite(rd_bands_dt, file = file.path(results_path, "GeneralizedPain_ACR2016_STRICT_Global_PastWeek_RD_Bands.csv"))

# --- File 3: Integrated Difference (AUC over 5 to 100) ---
auc_pt_prop <- trapz_between(rd_curve$Age, rd_curve$estimate, BKN[1], BKN[2])

auc_dt <- data.table(
  Age_Min = BKN[1], Age_Max = BKN[2],
  AUC_RD_pp_years = auc_pt_prop * 100,
  AUC_RD_pp_years_L = NA_real_, AUC_RD_pp_years_U = NA_real_,
  Mean_RD_pp = (auc_pt_prop * 100) / (BKN[2] - BKN[1]),
  Mean_RD_pp_L = NA_real_, Mean_RD_pp_U = NA_real_,
  GroupA = gA, GroupB = gB
)

fwrite(auc_dt, file = file.path(results_path, "GeneralizedPain_ACR2016_STRICT_Global_PastWeek_IntegratedDifference.csv"))

# --- File 4: Differences by age band (RD & RR) ---
diffs_list <- list()

rd_overall <- avg_comparisons(
  fit, variables = "Sex", type = "response",
  re.form = NA, comparison = "difference",
  newdata = datagrid(
    model = fit,
    Age = age_grid,
    time_window_category = "past_week",
    case_definition = CASE_PRED
  )
)

rr_overall <- avg_comparisons(
  fit, variables = "Sex", type = "response",
  re.form = NA, comparison = "ratio",
  newdata = datagrid(
    model = fit,
    Age = age_grid,
    time_window_category = "past_week",
    case_definition = CASE_PRED
  )
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
    fit, variables = "Sex", type = "response",
    re.form = NA, comparison = "difference",
    newdata = datagrid(
      model = fit,
      Age = br,
      time_window_category = "past_week",
      case_definition = CASE_PRED
    )
  )
  
  rr_band <- avg_comparisons(
    fit, variables = "Sex", type = "response",
    re.form = NA, comparison = "ratio",
    newdata = datagrid(
      model = fit,
      Age = br,
      time_window_category = "past_week",
      case_definition = CASE_PRED
    )
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
fwrite(diffs_dt, file = file.path(results_path, "GeneralizedPain_ACR2016_STRICT_Global_PastWeek_Differences.csv"))

# --- File 5: Rep ages (RD & RR at REP_AGES) ---
rd_ages <- comparisons(
  fit,
  newdata   = datagrid(
    model = fit,
    Age = REP_AGES,
    time_window_category = "past_week",
    case_definition = CASE_PRED
  ),
  variables = "Sex", type = "response",
  re.form   = NA, comparison = "difference"
)

rr_ages <- comparisons(
  fit,
  newdata   = datagrid(
    model = fit,
    Age = REP_AGES,
    time_window_category = "past_week",
    case_definition = CASE_PRED
  ),
  variables = "Sex", type = "response",
  re.form   = NA, comparison = "ratio"
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

fwrite(rep_dt, file = file.path(results_path, "GeneralizedPain_ACR2016_STRICT_Global_PastWeek_RepAges.csv"))

# --- File 6: Slopes by age band (pp/decade) ---
slope_rows <- list()

for (nm in names(BANDS)) {
  win <- BANDS[[nm]]
  
  slp <- avg_slopes(
    fit,
    newdata = datagrid(
      model = fit,
      Age = win,
      time_window_category = "past_week",
      case_definition = CASE_PRED
    ),
    by      = "Sex",
    type    = "response",
    re.form = NA,
    variable = "Age"
  )
  
  slp_diff <- avg_comparisons(
    fit,
    newdata = datagrid(
      model = fit,
      Age = win,
      time_window_category = "past_week",
      case_definition = CASE_PRED
    ),
    variables = list(Age = 1),
    by = "Sex",
    type = "response",
    re.form = NA
  )
  
  slope_rows[[length(slope_rows)+1]] <- data.table(
    Window=nm, Age_L=win[1], Age_H=win[2],
    Group=gA, Slope_pp_per_decade=slp[slp$Sex==gA,]$estimate*1000,
    L=slp[slp$Sex==gA,]$conf.low*1000, U=slp[slp$Sex==gA,]$conf.high*1000
  )
  
  slope_rows[[length(slope_rows)+1]] <- data.table(
    Window=nm, Age_L=win[1], Age_H=win[2],
    Group=gB, Slope_pp_per_decade=slp[slp$Sex==gB,]$estimate*1000,
    L=slp[slp$Sex==gB,]$conf.low*1000, U=slp[slp$Sex==gB,]$conf.high*1000
  )
  
  slope_rows[[length(slope_rows)+1]] <- data.table(
    Window=nm, Age_L=win[1], Age_H=win[2],
    Group=paste0(gA,"-",gB),
    Slope_pp_per_decade=slp_diff$estimate*1000,
    L=slp_diff$conf.low*1000, U=slp_diff$conf.high*1000
  )
}

slopes_dt <- rbindlist(slope_rows)
fwrite(slopes_dt, file = file.path(results_path, "GeneralizedPain_ACR2016_STRICT_Global_PastWeek_Slopes.csv"))

# --- File 7: Peak RD and Age x Sex interaction LRT ---
peak_info <- rd_bands_dt[which.max(abs(RD))]

excl0 <- (rd_bands_dt$RD_L > 0 & rd_bands_dt$RD_U > 0) | (rd_bands_dt$RD_L < 0 & rd_bands_dt$RD_U < 0)
r <- rle(excl0)

if (any(r$values)) {
  idx_run <- which.max(ifelse(r$values, r$lengths, 0))
  end_pos <- cumsum(r$lengths)[idx_run]
  start_pos <- end_pos - r$lengths[idx_run] + 1
  clear_start <- age_grid[start_pos]; clear_end <- age_grid[end_pos]
} else {
  clear_start <- NA_real_; clear_end <- NA_real_
}

# interaction test (remove Age x Sex)
form_add <- update(form_full, . ~ . - ns(Age, df=3):Sex)
fit_add <- update(fit, formula. = form_add)
p_interaction <- as.numeric(anova(fit_add, fit, test = "LRT")$`Pr(>Chisq)`[2])

peaks_dt <- data.table(
  Peak_RD_pp = peak_info$RD * 100,
  Peak_RD_pp_L = peak_info$RD_L * 100,
  Peak_RD_pp_U = peak_info$RD_U * 100,
  Peak_Age = peak_info$Age,
  Peak_Age_L = NA_real_, Peak_Age_U = NA_real_,
  ClearIntervalStart = clear_start, ClearIntervalEnd = clear_end,
  Interaction_p = p_interaction, Boots_Successful = NA_integer_
)

fwrite(peaks_dt, file = file.path(results_path, "_Global_PastWeek_Peaks_And_Interaction.csv"))

message(sprintf(
  "SUCCESS: GeneralizedPain_ACR2016_STRICT IPD analysis complete. Standardized predictions use case_definition = '%s'. Files saved to PainIntensity_NumSites.",
  CASE_PRED
))

# ===================================================================
# PLOTTING SCRIPT (unchanged aesthetics, now reading the updated results file)
# ===================================================================

# --- Libraries ---
library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(ggthemes)
library(cowplot)

# --- Theme (transparent) ---
theme_Publication <- function(base_size = 6.25, base_family = "sans") {
  theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_blank(),
      text = element_blank(),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background  = element_rect(fill = "transparent", colour = NA),
      panel.border     = element_rect(colour = NA),
      axis.title       = element_blank(),
      axis.line        = element_line(color = "black", size = 0.25),
      axis.ticks       = element_line(color = "black", size = 0.25),
      axis.ticks.length= unit(1.75, "pt"),
      axis.text.x      = element_text(size = 6.25, margin = margin(t = 1.5, unit = "pt")),
      axis.text.y      = element_text(size = 6.25, margin = margin(r = 1.5, unit = "pt")),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "none",
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text       = element_blank()
    )
}

# --- Paths ---
infile  <- "/Users/Patty/Desktop/LifespanPain/Results/PainIntensity_NumSites/GeneralizedPain_ACR2016_STRICT_Global_PastWeek_Results.csv"
outfile <- "/Users/Patty/Desktop/LifespanPain/Figures/Fig_2_Intensity_NumSites/GeneralizedPain_ACR2016_STRICT_PastWeek_combined_plot.pdf"
dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)

# --- Load ---
dat <- fread(infile)
dat$Sex <- factor(dat$Sex)

# --- Peak rate-of-change per sex (for markers) ---
peak_roc <- dat %>%
  group_by(Sex) %>%
  arrange(Age) %>%
  summarize(
    peak_age  = round(Age[which.max(Rate_Of_Change_Annual)]),
    peak_prob = Predicted_Prob[which.max(Rate_Of_Change_Annual)],
    peak_rate = max(Rate_Of_Change_Annual, na.rm = TRUE),
    .groups = "drop"
  )

# --- Probability trajectory ---
plot_prob <- ggplot(dat, aes(x = Age, y = Predicted_Prob, color = Sex)) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI, fill = Sex), alpha = 0.1, linetype = 0) +
  geom_line(size = 1) +
  labs(title = "Generalized pain", x = "Age (Years)", y = "Probability of Generalized Pain") +
  scale_x_continuous(breaks = seq(10, 100, by = 10), expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  coord_cartesian(xlim = c(5, 100)) +
  scale_color_manual(values = c("#FB7D25", "#3AADBB")) +
  scale_fill_manual(values = c("#FB7D25", "#3AADBB")) +
  theme_Publication(base_size = 6.25) +
  geom_point(
    data = peak_roc,
    aes(x = peak_age, y = peak_prob, color = Sex),
    size = 1.5, shape = 21, fill = "white", stroke = 1
  ) +
  geom_text(
    data = peak_roc,
    aes(x = peak_age, y = peak_prob, label = round(peak_age, 1), color = Sex),
    vjust = -1, size = 2.5
  )

# --- Rate of change ---
plot_rate <- ggplot(dat, aes(x = Age, y = Rate_Of_Change_Annual, color = Sex)) +
  geom_segment(
    data = subset(peak_roc, Sex == levels(dat$Sex)[1]),
    aes(x = peak_age, xend = peak_age, y = 0, yend = peak_rate),
    linetype = "21", color = "#FB7D25", size = 0.5
  ) +
  geom_segment(
    data = subset(peak_roc, Sex == levels(dat$Sex)[2]),
    aes(x = peak_age, xend = peak_age, y = 0, yend = peak_rate),
    linetype = "21", color = "#3AADBB", size = 0.5
  ) +
  geom_hline(yintercept = 0, color = "gray", size = 0.5) +
  geom_line(size = 0.5) +
  labs(title = "Rate of Change", x = "", y = "Rate of Change") +
  scale_x_continuous(breaks = seq(10, 100, by = 10), expand = c(0, 0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3), labels = percent_format(accuracy = 0.1)) +
  coord_cartesian(xlim = c(5, 100)) +
  scale_color_manual(values = c("#FB7D25", "#3AADBB")) +
  theme_Publication(base_size = 6.25) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x  = element_blank()
  )

# --- Combine & save ---
combined_plot <- plot_grid(
  plot_prob, plot_rate,
  ncol = 1, align = "v", axis = "lr",
  rel_heights = c(3, 1.15)
)

ggsave(
  filename = outfile,
  plot = combined_plot,
  width = 2, height = 1.35, units = "in",
  bg = "transparent"
)

cat("Saved GeneralizedPain_ACR2016_STRICT combined plot to:", outfile, "\n")
