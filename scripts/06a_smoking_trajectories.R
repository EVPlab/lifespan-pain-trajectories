# -------------------------------
# Smoking (continuous age): 3 levels (Non-Smoker, Ex-Smoker, Smoker)
# Per-site TW (fixed to past_month), Study+Region RE, equal-sex averaging,
# 5–100 domain, marginaleffects
# -------------------------------

# Libraries
library(lme4)
library(dplyr)
library(data.table)
library(splines)
library(ggplot2)
library(scales)
library(foreach)
library(doParallel)
library(ggthemes)
library(marginaleffects)
library(tidyr)
library(RColorBrewer)

set.seed(123)
options(contrasts = c("contr.sum", "contr.poly"))

# Config
tw_levels <- c("past_week","past_month","past_3_plus_months","chronic_explicit")
AGE_MIN <- 13
AGE_MAX <- 100

pain_sites <- c('Pain','Headache','FacialPain','NeckShoulderPain','ChestPain',
                'StomachAbdominalPain','BackPain','HipPain','KneePain',
                'HandPain', 'FootPain', 'ElbowPain')

# pain_sites <- c('FacialPain')  # for testing

# Columns to load
time_window_cols <- paste0("time_window_category_", pain_sites)
required_columns <- c("Study","Global_Region","Age","Sex",
                      "SmokingStatus","Smoking_category",   # <- added SmokingStatus (+ fallback)
                      "time_window_category", time_window_cols, pain_sites)

data <- fread("/Users/Patty/Desktop/R/LifespanPain/Data/PainData_CS.csv",
              select = required_columns)

# ---------------- Harmonize Smoking to 3 levels ----------------
# Prefer 3-level SmokingStatus when present; otherwise map 2-level Smoking_category
data <- data %>%
  mutate(
    Smoking_category = trimws(as.character(Smoking_category)),
    SmokingStatus    = trimws(as.character(SmokingStatus)),
    SmokingStatus = case_when(
      !is.na(SmokingStatus) ~ SmokingStatus,                # keep 3-level when present
      Smoking_category %in% c("Smoker") ~ "Smoker", # common encodings (adjust if needed)
      Smoking_category %in% c("Non-Smoker") ~ "Non-Smoker",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    Study         = factor(Study),
    Global_Region = factor(trimws(as.character(Global_Region))),
    Sex           = factor(Sex),
    SmokingStatus = factor(SmokingStatus,
                           levels = c("Non-Smoker","Ex-Smoker","Smoker"))
  )

# Spline domain (fixed)
age_all <- data$Age[is.finite(data$Age)]
BKN <- c(AGE_MIN, AGE_MAX)

gc()

# Parallel
total_cores <- detectCores()
num_cores   <- max(1, min(total_cores - 1, 4))
cl          <- makeCluster(num_cores)
registerDoParallel(cl)
clusterExport(cl, varlist = c("BKN","tw_levels","AGE_MIN","AGE_MAX"), envir = environment())

# Theme
theme_Publication <- function(base_size = 7, base_family = "sans") {
  theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_blank(), text = element_blank(),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background  = element_rect(fill = "transparent", colour = NA),
      panel.border     = element_rect(colour = NA),
      axis.title = element_blank(),
      axis.line  = element_line(color = "black", size = 0.25),
      axis.ticks = element_line(color = "black", size = 0.25),
      axis.ticks.length = unit(1.75, "pt"),
      axis.text.x = element_text(size = 7, margin = margin(t = 1.5, unit = "pt")),
      axis.text.y = element_text(size = 7, margin = margin(r = 1.5, unit = "pt")),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      legend.position = "none",
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_blank()
    )
}

# Model + plot
model_and_plot <- function(
    pain_site,
    data,
    save_path = "~/Desktop/LifespanPain/Figures/Fig_3_Health/Exposure_trajectories/Smoking/",
    csv_path  = "/Users/Patty/Desktop/LifespanPain/Results/Exposure_trajectories/"
) {
  library(lme4); library(data.table); library(dplyr); library(splines)
  library(ggplot2); library(ggthemes); library(marginaleffects); library(tidyr); library(RColorBrewer)
  
  # Per-site time window
  tw_col <- paste0("time_window_category_", pain_site)
  if (!tw_col %in% names(data)) tw_col <- "time_window_category"
  
  # Filter rows for this site
  df <- as.data.table(data)[
    !is.na(get(pain_site)) & is.finite(Age) &
      !is.na(Study) & !is.na(SmokingStatus) &
      !is.na(Sex) & !is.na(Global_Region) & !is.na(get(tw_col))
  ]
  
  # Site-specific TW levels actually present (avoid empty levels)
  tw_fit <- tw_levels[tw_levels %in% unique(as.character(df[[tw_col]]))]
  df[, tw := factor(get(tw_col), levels = tw_fit)]
  
  # Drop unused levels on key factors
  df[, `:=`(
    SmokingStatus = droplevels(SmokingStatus),
    Sex           = droplevels(Sex)
  )]
  
  # Ensure binary numeric outcome (0/1)
  if (!is.numeric(df[[pain_site]])) {
    df[[pain_site]] <- as.numeric(as.character(df[[pain_site]]))
  }
  
  # GLMM: ns(Age)*SmokingStatus + Sex + tw + RE(Study, Region)
  # form <- as.formula(paste0(
  #   "`", pain_site, "` ~ ns(Age, df = 2) * SmokingStatus + ",
  #   "Sex + tw + (1 | Study) + (1 | Global_Region) + (1 | Study:Global_Region)"
  # ))
  
  form <- as.formula(paste0(
    "`", pain_site, "` ~ ns(Age, df = 2) * SmokingStatus + ",
    "Sex + tw + (1 | Global_Region) + (1 | Study)"
  ))
  
  fit <- tryCatch(
    glmer(form, data = df, family = binomial("logit"),
          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))),
    error = function(e) { message("glmer error: ", e$message); NULL }
  )
  if (is.null(fit)) return(NULL)
  
  # Prediction grid: tw = past_month (if present); equal-sex averaging
  age_seq     <- seq(AGE_MIN, AGE_MAX, length.out = 150)
  sex_levels  <- levels(df$Sex)
  smoke_lvls  <- levels(df$SmokingStatus)  # c("Non-Smoker","Ex-Smoker","Smoker") subset actually present
  target_tw   <- if ("past_month" %in% levels(df$tw)) "past_month" else levels(df$tw)[1]
  
  grid <- tidyr::crossing(
    Age = age_seq,
    SmokingStatus = factor(smoke_lvls, levels = smoke_lvls),
    Sex = factor(sex_levels, levels = sex_levels),
    tw  = factor(target_tw, levels = levels(df$tw))
  )
  
  # Equal-sex weights within Age × SmokingStatus
  sex_w <- setNames(rep(1 / length(sex_levels), length(sex_levels)), sex_levels)
  grid <- grid |>
    mutate(w = sex_w[as.character(Sex)]) |>
    group_by(Age, SmokingStatus) |>
    mutate(w = if (sum(w) > 0) w / sum(w) else 1 / n()) |>
    ungroup()
  
  # Population-average predictions
  ap <- marginaleffects::avg_predictions(
    fit,
    newdata = grid,
    by = c("Age","SmokingStatus"),
    type = "response",
    re.form = NA,
    wts = "w"
  )
  res <- as.data.table(ap)
  setnames(res, c("estimate","conf.low","conf.high","std.error"),
           c("Predicted_Prob","Lower_CI","Upper_CI","SE"))
  setorder(res, SmokingStatus, Age)
  
  # Rate of change (pp/year)
  res[, Rate_Of_Change := c(NA_real_, diff(Predicted_Prob) / diff(Age)), by = SmokingStatus]
  
  # Save CSV
  dir.create(csv_path, showWarnings = FALSE, recursive = TRUE)
  out_csv <- file.path(csv_path, paste0(pain_site, "_Smoking_continuous_predictions.csv"))
  fwrite(res[, .(Age, SmokingStatus, Predicted_Prob, Lower_CI, Upper_CI, Rate_Of_Change)], out_csv)
  message("Saved: ", out_csv)
  
  # --- Age Band RR vs. Non-Smoker ---
  message("Calculating: Smoking differences by age band for ", pain_site)
  
  age_bands <- list(
    "overall" = c(AGE_MIN, AGE_MAX),
    "13-17"   = c(13, 17),
    "18-35"   = c(18, 35),
    "35-50"   = c(35, 50),
    "50-65"   = c(50, 65),
    "65-80"   = c(65, 80),
    "80-100"  = c(80, 100)
  )
  
  all_diffs <- list()
  
  for (band_name in names(age_bands)) {
    df_band <- df[Age >= age_bands[[band_name]][1] & Age <= age_bands[[band_name]][2]]
    if (nrow(df_band) == 0) next
    
    df_band$SmokingStatus <- factor(df_band$SmokingStatus,
                                    levels = c("Non-Smoker", "Ex-Smoker", "Smoker"))
    
    comp <- avg_comparisons(
      fit,
      newdata    = df_band,
      variables  = "SmokingStatus",
      comparison = "ratio",
      type       = "response",
      re.form    = NA
    )
    comp_dt <- as.data.table(comp)
    comp_dt[, age_band := band_name]
    all_diffs[[band_name]] <- comp_dt[, .(age_band, contrast, RR = estimate, RR_L = conf.low, RR_U = conf.high)]
  }
  
  rr_final <- rbindlist(all_diffs)
  diff_csv  <- file.path(csv_path, paste0(pain_site, "_Smoking_PastMonth_Differences.csv"))
  fwrite(rr_final, file = diff_csv)
  message("Saved: ", diff_csv)
  
  # after building `res` from avg_predictions(...), keep only two lines to plot
  res_plot <- res %>% dplyr::filter(SmokingStatus %in% c("Non-Smoker","Smoker"))
  
  colors <- RColorBrewer::brewer.pal(5, "Oranges")
  color_map <- c("Non-Smoker" = colors[2], "Smoker" = colors[5])
  
  plot_prob <- ggplot(res_plot, aes(Age, Predicted_Prob, color = SmokingStatus, fill = SmokingStatus)) +
    geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.10, linetype = 0) +
    geom_line(linewidth = 0.8) +
    scale_x_continuous(breaks = seq(5, 100, by = 10), expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    coord_cartesian(xlim = c(5, 100)) +
    scale_color_manual(values = color_map, drop = TRUE) +
    scale_fill_manual(values = color_map, drop = TRUE) +
    theme_Publication(base_size = 7) +
    theme(aspect.ratio = 2/5)
  
  # rate-of-change
  res_plot_roc <- res_plot %>%
    group_by(SmokingStatus) %>%
    mutate(Rate_Of_Change = c(NA_real_, diff(Predicted_Prob) / diff(Age)))
  
  plot_roc <- ggplot(res_plot_roc, aes(Age, Rate_Of_Change, color = SmokingStatus)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray", linewidth = 0.7) +
    geom_line(linewidth = 0.8) +
    scale_x_continuous(breaks = seq(5, 100, by = 10), expand = c(0, 0)) +
    scale_y_continuous(labels = function(y) paste0(round(100*y, 1), " pp/yr")) +
    scale_color_manual(values = color_map, drop = TRUE) +
    theme_Publication(base_size = 7) +
    theme(aspect.ratio = 2/5)
  
  
  dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(save_path, paste0(pain_site, "_probability_plot_Smoking.pdf")),
         plot_prob, width = 2, height = 1, units = "in", dpi = 300)
  ggsave(file.path(save_path, paste0(pain_site, "_rate_of_change_plot_Smoking.pdf")),
         plot_roc, width = 2, height = 1, units = "in", dpi = 300)
  
  invisible(NULL)
}

# Parallel run
foreach(
  pain_site = pain_sites,
  .packages = c("lme4","dplyr","splines","ggplot2","data.table",
                "ggthemes","marginaleffects","tidyr","scales","RColorBrewer")
) %dopar% {
  model_and_plot(pain_site, data)
}

stopCluster(cl)
rm(data); gc()
cat("Smoking (3-level model) modeling and plotting completed.\n")
