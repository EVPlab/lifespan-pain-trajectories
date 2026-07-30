# =====================================================================
# run_demo.R
# Minimal, self-contained demonstration of the core age-trajectory model
# used in "Pain across the lifespan". Runs on the synthetic dataset so no
# real (restricted) data are needed.
#
# Steps:
#   1. Load ./data/PainData_CS_synthetic.csv (run generate_synthetic_data.R first)
#   2. Fit the primary GLMM for one pain site (BackPain):
#        outcome ~ ns(Age, df = 3) * Sex + recall + case_def + (1|Study) + (1|Region)
#   3. Predict age-specific prevalence by sex and save a table + a figure
#
# Usage:  Rscript demo/run_demo.R
# Expected run time: ~1-3 minutes on a normal desktop.
# =====================================================================

suppressPackageStartupMessages({
  library(data.table); library(lme4); library(splines); library(ggplot2)
})

site   <- "BackPain"
infile <- "data/PainData_CS_synthetic.csv"
outdir <- "demo/output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(infile))
  stop("Synthetic data not found. Run: Rscript demo/generate_synthetic_data.R")

d <- fread(infile)
d[, Sex := factor(Sex, levels = c(1, 2), labels = c("Male", "Female"))]
d[, Sex := relevel(Sex, ref = "Female")]
for (v in c("Study","Global_Region","time_window_category","case_definition"))
  d[[v]] <- as.factor(d[[v]])
d[, time_window_category := relevel(time_window_category, ref = "past_month")]
d[, case_definition := relevel(case_definition, ref = "general")]

BKN <- c(5, 100)
cat(sprintf("Fitting GLMM for %s on %d synthetic rows ...\n", site, nrow(d)))
form <- as.formula(paste0(
  site, " ~ ns(Age, df = 3, Boundary.knots = BKN) * Sex + ",
  "time_window_category + case_definition + (1 | Study) + (1 | Global_Region)"))
fit <- glmer(form, data = d, family = binomial("logit"),
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

# Predict marginal prevalence over age, standardized to past_month / general
grid <- CJ(Age = seq(5, 100, by = 1), Sex = c("Female", "Male"))
grid[, `:=`(time_window_category = "past_month", case_definition = "general")]
grid[, Sex := relevel(factor(Sex, levels = levels(d$Sex)), ref = "Female")]
grid[, Predicted_Prevalence := predict(fit, newdata = grid, re.form = NA, type = "response")]

fwrite(grid[, .(Age, Sex, Predicted_Prevalence)],
       file.path(outdir, paste0(site, "_demo_trajectory.csv")))

p <- ggplot(grid, aes(Age, Predicted_Prevalence, colour = Sex)) +
  geom_line(linewidth = 1.1) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = paste0("DEMO (synthetic data): ", site, " prevalence across the lifespan"),
       y = "Predicted prevalence", x = "Age (years)") +
  theme_classic(base_size = 12)
ggsave(file.path(outdir, paste0(site, "_demo_trajectory.png")), p, width = 7, height = 4.5, dpi = 150)

cat("Done.\nExpected output:\n",
    "  - demo/output/", site, "_demo_trajectory.csv  (192 rows: Age x Sex predicted prevalence)\n",
    "  - demo/output/", site, "_demo_trajectory.png  (an inverted-U-like age curve by sex)\n", sep = "")
