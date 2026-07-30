# =====================================================================
# generate_synthetic_data.R
# Creates a SMALL, FULLY SYNTHETIC dataset that mirrors the column
# structure expected by the analysis scripts in this repository.
#
# IMPORTANT: This file contains NO real participant data. Every value is
# randomly generated. It exists only so that reviewers and readers can run
# the analysis pipeline end-to-end. It is NOT derived from, or a sample of,
# any contributing study, and therefore does not fall under any data-use
# agreement.
#
# Output: ./data/PainData_CS_synthetic.csv
# Usage:  Rscript demo/generate_synthetic_data.R
# =====================================================================

set.seed(2026)
suppressPackageStartupMessages(library(data.table))

out_dir <- "data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Design of the synthetic "studies" ------------------------------
regions <- c("Central/Southern Asia","Eastern/South-Eastern Asia","Europe/Northern America",
             "Latin America/Caribbean","Northern Africa/Western Asia","Oceania",
             "Sub-Saharan Africa","Australia/New Zealand","Eastern Asia")
windows <- c("past_week","past_month","past_3_12months")
casedef <- c("general","chronic_explicit")

n_studies <- 12
studies <- data.table(
  Study         = sprintf("SYN%02d", seq_len(n_studies)),
  Global_Region = sample(regions, n_studies, replace = TRUE),
  time_window_category = sample(windows, n_studies, replace = TRUE, prob = c(.2,.6,.2)),
  case_definition = sample(casedef, n_studies, replace = TRUE, prob = c(.85,.15)),
  HDI  = round(runif(n_studies, 0.45, 0.95), 3),
  Year = sample(1995:2024, n_studies, replace = TRUE),
  study_intercept = rnorm(n_studies, 0, 0.35),   # study-level random effect (logit)
  n   = sample(800:2500, n_studies, replace = TRUE)
)

# ---- Expand to participants -----------------------------------------
mk <- function(s) {
  n <- s$n
  age <- pmin(pmax(round(rgamma(n, shape = 6, scale = 8) + 5), 5), 95)
  sex <- sample(c(1L, 2L), n, replace = TRUE)              # 1 = Male, 2 = Female
  data.table(
    Study = s$Study, Global_Region = s$Global_Region,
    time_window_category = s$time_window_category,
    case_definition = s$case_definition,
    HDI = s$HDI, Year = s$Year,
    Age = age, Sex = sex,
    BMI_Category = sample(c("Normal","Overweight","Obese"), n, replace = TRUE, prob = c(.45,.35,.20)),
    SmokingStatus = sample(c("Non-Smoker","Ex-Smoker","Smoker"), n, replace = TRUE, prob = c(.55,.20,.25)),
    Income_class_quintile = sample(paste0("Q", 1:5), n, replace = TRUE),
    .study_int = s$study_intercept
  )
}
dt <- rbindlist(lapply(seq_len(n_studies), function(i) mk(studies[i])))
dt[, Smoking_category := fifelse(SmokingStatus == "Smoker", "Smoker", "Non-Smoker")]

# ---- Simulate age-dependent pain outcomes (illustrative shapes) -----
z  <- function(x) (x - 50) / 25
fem <- as.integer(dt$Sex == 2L)
a   <- z(dt$Age); ri <- dt$.study_int
inv <- plogis
# site-specific logit models: intercept + linear/quadratic age + sex + study RE
sites <- list(
  Pain                 = -0.2 + 0.35*a + 0.10*fem + ri,
  JointPain            = -0.8 + 0.85*a + 0.10*fem + ri,
  Headache             = -0.9 - 0.35*a - 0.30*a^2 + 0.45*fem + ri,
  FacialPain           = -3.1 - 0.10*a + 0.20*fem + ri,
  NeckShoulderPain     = -1.2 + 0.45*a + 0.30*fem + ri,
  ChestPain            = -2.6 + 0.30*a + 0.05*fem + ri,
  StomachAbdominalPain = -1.9 - 0.15*a + 0.25*fem + ri,
  BackPain             = -0.7 + 0.55*a + 0.15*fem + ri,
  HipPain              = -2.2 + 0.80*a + 0.20*fem + ri,
  KneePain             = -1.6 + 0.95*a + 0.15*fem + ri,
  HandPain             = -2.0 + 0.70*a + 0.25*fem + ri,
  ElbowPain            = -2.7 + 0.35*a + 0.05*fem + ri,
  FootPain             = -2.1 + 0.45*a + 0.15*fem + ri
)
for (nm in names(sites)) dt[[nm]] <- rbinom(nrow(dt), 1L, inv(sites[[nm]]))

dt[, .study_int := NULL]

# Column order roughly matching the real analytic file
setcolorder(dt, c("Study","Sex","Age","Global_Region","time_window_category","case_definition",
                  "HDI","Year","BMI_Category","SmokingStatus","Smoking_category","Income_class_quintile",
                  names(sites)))

outfile <- file.path(out_dir, "PainData_CS_synthetic.csv")
fwrite(dt, outfile)
cat(sprintf("Wrote %s  (%d rows, %d synthetic studies)\n", outfile, nrow(dt), n_studies))
