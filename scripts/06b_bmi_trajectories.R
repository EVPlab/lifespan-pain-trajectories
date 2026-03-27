# Load necessary libraries
library(lme4)
library(dplyr)
library(data.table)  # More memory-efficient
library(splines)
library(ggplot2)
library(scales)      # For percentage formatting
library(foreach)
library(readr)
library(parallel)
library(ggthemes)
library(doParallel)
library(marginaleffects)

options(contrasts = c("contr.sum", "contr.poly"))

tw_levels <- c("past_week","past_month","past_3_plus_months","chronic_explicit")

pain_sites <- c('Headache','FacialPain','NeckShoulderPain','ChestPain',
                'StomachAbdominalPain','BackPain','HipPain','KneePain',
                'HandPain', 'FootPain', 'ElbowPain')

# Columns
time_window_cols <- paste0("time_window_category_", pain_sites)
required_columns <- c("Study","Global_Region","Age","BMI_Category","Sex",
                      "time_window_category", time_window_cols, pain_sites)

data <- data.table::fread("/Users/Patty/Desktop/R/LifespanPain/Data/PainData_CS.csv",
                          select = required_columns)
# data <- data %>% sample_frac(.5)

# Global spline knots (fixed across all fits)
age_all <- data$Age[is.finite(data$Age)]
BKN <- c(5, 100)
IK2 <- as.numeric(quantile(age_all, 0.50, na.rm = TRUE))          # 1 internal knot (median)
message(sprintf("Knots: Boundary [%.2f, %.2f]; Internal: %.2f", BKN[1], BKN[2], IK2))


# Convert columns to factors (assuming your BMI_Category are already labeled)
# Example levels: Normal, Overweight, Obese
data[, Study := factor(Study)]
data[, BMI_Category := factor(
  BMI_Category, 
  levels = c("Normal", "Overweight", "Obese")
)]
data[, Sex := factor(Sex)]

gc()  # Garbage collection

# -------------------------------
# 2. Parallel Setup
# -------------------------------
total_cores <- detectCores()
num_cores   <- min(total_cores - 1, 4)
cl          <- makeCluster(num_cores)
registerDoParallel(cl)

# -------------------------------
# 3. Model and Plot Function (Controlling for Sex)
# -------------------------------
model_and_plot <- function(
    pain_site,
    data,
    save_path = "~/Desktop/LifespanPain/Figures/Fig_3_Health/Exposure_trajectories/BMI/",
    csv_path  = "/Users/Patty/Desktop/LifespanPain/Results/Exposure_trajectories/"
) {
  library(lme4); library(data.table); library(ggplot2); library(dplyr); library(splines); library(ggthemes)
  
  # Theme (unchanged)
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
        legend.position = "none", strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
        strip.text = element_blank()
      )
  }
  
  # --- Filter & construct tw ---
  tw_col <- paste0("time_window_category_", pain_site)
  if (!tw_col %in% names(data)) tw_col <- "time_window_category"
  
  df <- as.data.table(data)[
    !is.na(get(pain_site)) & !is.na(Age) & !is.na(Study) &
      !is.na(BMI_Category) & !is.na(Sex) & !is.na(Global_Region) & !is.na(get(tw_col))
  ][
    , `:=`(
      Study = factor(Study),
      Global_Region = factor(trimws(as.character(Global_Region))),
      BMI_Category = factor(BMI_Category, levels = c("Normal","Overweight","Obese")),
      Sex = factor(Sex)
    )
  ]
  
  # Site-specific TW levels actually present (keeps original order)
  tw_fit <- tw_levels[tw_levels %in% unique(as.character(df[[tw_col]]))]
  # Make tw using only available levels
  df[, tw := factor(get(tw_col), levels = tw_fit)]
  
  # Ensure binary numeric outcome
  if (!is.numeric(df[[pain_site]])) {
    df[[pain_site]] <- as.numeric(as.character(df[[pain_site]]))
  }
  df <- df[is.finite(Age)]
  
  # --- Fit GLMM: ns(Age) * BMI + Sex + tw + (1|Study) + (1|Region) ---

  # form <- as.formula(paste0(
  #   "`", pain_site, "` ~ ns(Age, df = 2) * BMI_Category + ",
  #   "Sex + tw + (1 | Study) + (1 | Global_Region) + (1 | Study:Global_Region)"
  # ))
  
  form <- as.formula(paste0(
    "`", pain_site, "` ~ ns(Age, df = 2) * BMI_Category + ",
    "Sex + tw + (1 | Global_Region) + (1 | Study)"
  ))
  
  fit <- tryCatch(
    glmer(form, data = df, family = binomial(link="logit"),
          control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1000))),
    error = function(e) { message("glmer error: ", e$message); NULL }
  )
  if (is.null(fit)) return(NULL)
  
  # --- Prediction grid: tw fixed to past_month; NO region in grid; equal-sex averaging ---
  age_seq <- seq(5, 100, length.out = 150)
  sex_levels <- levels(df$Sex)
  bmi_levels <- levels(df$BMI_Category)
  target_tw <- if ("past_month" %in% levels(df$tw)) "past_month" else {
    df[, .N, by = tw][order(-N)][["tw"]][1]
  }
  
  grid <- tidyr::crossing(
    Age = age_seq,
    BMI_Category = factor(bmi_levels, levels = bmi_levels),
    Sex = factor(sex_levels, levels = sex_levels),
    tw  = factor(target_tw, levels = tw_fit)
  )
  
  # Equal-sex weights within Age × BMI; normalize locally
  sex_w <- setNames(rep(1/length(sex_levels), length(sex_levels)), sex_levels)
  grid <- grid |>
    mutate(w = sex_w[as.character(Sex)]) |>
    group_by(Age, BMI_Category) |>
    mutate(w = if (sum(w) > 0) w/sum(w) else 1/n()) |>
    ungroup()
  
  # Population-average predictions (exclude RE)
  ap <- marginaleffects::avg_predictions(
    fit,
    newdata = grid,
    by = c("Age","BMI_Category"),
    type = "response",
    re.form = NA,
    wts = "w"
  )
  res <- as.data.table(ap)
  setnames(res, c("estimate","conf.low","conf.high","std.error"),
           c("Predicted_Prob","Lower_CI","Upper_CI","SE"))
  setorder(res, BMI_Category, Age)
  
  # Rafgte of change (pp per year)
  res[, Rate_Of_Change := c(NA_real_, diff(Predicted_Prob) / diff(Age)), by = BMI_Category]
  
  # Save CSV
  out_csv <- file.path(csv_path, paste0(pain_site, "_BMI_continuous_predictions.csv"))
  fwrite(res[, .(Age, BMI_Category,
                 Predicted_Prob, Lower_CI, Upper_CI, Rate_Of_Change)], out_csv)
  message("Saved: ", out_csv)
  
  # --- Age Band RR vs. Normal BMI ---
  message("Calculating: BMI differences by age band for ", pain_site)
  
  age_bands <- list(
    "overall" = c(5, 100),
    "5-17"    = c(5, 17),
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
    
    df_band$BMI_Category <- factor(df_band$BMI_Category,
                                   levels = c("Normal", "Overweight", "Obese"))
    
    comp <- avg_comparisons(
      fit,
      newdata    = df_band,
      variables  = "BMI_Category",
      comparison = "ratio",
      type       = "response",
      re.form    = NA
    )
    comp_dt <- as.data.table(comp)
    comp_dt[, age_band := band_name]
    all_diffs[[band_name]] <- comp_dt[, .(age_band, contrast, RR = estimate, RR_L = conf.low, RR_U = conf.high)]
  }
  
  rr_final <- rbindlist(all_diffs)
  diff_csv  <- file.path(csv_path, paste0(pain_site, "_BMI_PastMonth_Differences.csv"))
  fwrite(rr_final, file = diff_csv)
  message("Saved: ", diff_csv)
  
  # --- Plotting (same aesthetics) ---
  # Color map
  suppressPackageStartupMessages(library(RColorBrewer))
  orange_palette <- RColorBrewer::brewer.pal(5, "Oranges")
  color_map_plot <- c("Obese" = orange_palette[5],
                      "Overweight" = orange_palette[3],
                      "Normal" = orange_palette[2])
  
  plot_prob <- ggplot(res, aes(Age, Predicted_Prob, color = BMI_Category)) +
    geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI, fill = BMI_Category),
                alpha = 0.10, linetype = 0) +
    geom_line(linewidth = 0.8) +
    scale_x_continuous(breaks = seq(10, 100, by = 10), expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    coord_cartesian(xlim = c(5, 92)) +
    scale_color_manual(values = color_map_plot) +
    scale_fill_manual(values = color_map_plot) +
    theme_Publication(base_size = 7) +
    theme(aspect.ratio = 2/5)
  
  # Peak ROC markers (optional; uses finite differences)
  peak_roc <- res[ , .SD[which.max(Rate_Of_Change)], by = BMI_Category]
  plot_prob <- plot_prob +
    geom_point(data = peak_roc,
               aes(x = Age, y = Predicted_Prob, color = BMI_Category),
               size = 1.5, shape = 21, fill = "white", stroke = 1)
  
  plot_roc <- ggplot(res, aes(Age, Rate_Of_Change, color = BMI_Category)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray", linewidth = 0.7) +
    geom_line(linewidth = 0.8) +
    scale_x_continuous(breaks = seq(10, 100, by = 10), expand = c(0, 0)) +
    scale_y_continuous(labels = function(y) paste0(round(100*y, 1), " pp/yr")) +
    scale_color_manual(values = color_map_plot) +
    theme_Publication(base_size = 7) +
    theme(aspect.ratio = 2/5)
  
  dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(save_path, paste0(pain_site, "_probability_plot_BMI.pdf")),
         plot_prob, width = 2, height = 1, units = "in", dpi = 300)
  ggsave(file.path(save_path, paste0(pain_site, "_rate_of_change_plot_BMI.pdf")),
         plot_roc, width = 2, height = 1, units = "in", dpi = 300)
  
  invisible(NULL)
}

# -------------------------------
# 4. Parallel Processing
# -------------------------------
total_cores <- parallel::detectCores()
num_cores   <- min(total_cores - 1, 4)
cl          <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

parallel::clusterExport(cl, varlist = c("BKN","IK2","tw_levels"),
                        envir = environment())

foreach::foreach(
  pain_site = pain_sites,
  .packages = c("lme4","dplyr","splines","ggplot2","data.table","ggthemes","marginaleffects","scales")
) %dopar% {
  model_and_plot(pain_site, data)
}

parallel::stopCluster(cl)
rm(data); gc()
cat("Modeling and plotting completed successfully.\n")

