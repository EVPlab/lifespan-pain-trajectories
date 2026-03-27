# Regional Script with 1 month reference + Annual Percent Change (APC)
# Updated: Adds case_definition covariate (with per-site override)
#          Prediction standardized to case_definition = "general" when available

# --- Step 0: Load Libraries and Configure Environment ---
library(lme4); library(dplyr); library(data.table); library(splines)
library(ggplot2); library(ggthemes); library(scales); library(forcats)
library(emmeans)
library(marginaleffects)   # for standardized averaging and derivatives

set.seed(123)
options(contrasts = c("contr.sum", "contr.poly"))

# --- CONFIGURATION ---
pain_sites_to_process <- c(
  "Pain", "JointPain", "BackPain", "Headache"
)
# case_definition settings (mirrors IPD script)
CASE_LEVELS <- c("general", "chronic_explicit")
PRED_CASEDEF_DEFAULT <- "general"

# PLEASE UPDATE these paths
data_path <- "/Users/Patty/Desktop/R/LifespanPain/Data/PainData_CS.csv"
results_output_path <- "/Users/Patty/Desktop/LifespanPain/Results/Regional_trajectories/"

main_regions_to_plot <- c(
  "East and SE Asia", 
  "Northern America", 
  "Western Europe", 
  "Eastern Europe",
  "Latin America", 
  "North Africa and West Asia", 
  "Sub-saharan Africa",
  "Central Asia"
)

la_missing_sites <- c("FacialPain")

region_colors <- c(
  "East and SE Asia"           = "#4B95C8",
  "Western Europe"             = "#49BC41",
  "Eastern Europe"             = "#BC417E",
  "Northern America"           = "#E03A3C",
  "Latin America"              = "#9241BC",
  "North Africa and West Asia" = "#BC9A41",
  "Sub-saharan Africa"         = "#1B9E77",
  "Central Asia"               = "#D95F02",
  "Global Average"             = "black"
)

# --- Step 1: Data Loading and Site-Agnostic Preparation ---
message("--- Loading and preparing data... ---")

time_window_cols <- paste0("time_window_category_", pain_sites_to_process)
case_def_cols    <- paste0("case_definition_", pain_sites_to_process)

required_columns <- c(
  "Study", "Sex", "Age", "Global_Region",
  "time_window_category",
  "case_definition",
  time_window_cols,
  case_def_cols,
  pain_sites_to_process
)

# Robust selection (skips columns that do not exist)
header_names <- names(fread(data_path, nrows = 0))
missing_cols <- setdiff(required_columns, header_names)
if (length(missing_cols) > 0) {
  message("NOTE: These requested columns were not found in the CSV and will be skipped:\n  ",
          paste(missing_cols, collapse = ", "))
}
select_cols <- intersect(required_columns, header_names)

data <- fread(data_path, select = select_cols)

# Ensure global case_definition exists even if not in file
if (!"case_definition" %in% names(data)) {
  data[, case_definition := NA_character_]
}

# Group regions
data[, Global_Region_Grouped := forcats::fct_other(
  trimws(as.character(Global_Region)),
  keep = main_regions_to_plot,
  other_level = "Other Regions"
)]
new_region_levels <- c(main_regions_to_plot, "Other Regions")
data[, Global_Region_Grouped := factor(Global_Region_Grouped, levels = new_region_levels)]

# Prepare base types
data[, `:=`(
  Study = as.factor(Study),
  Sex = as.factor(Sex),
  time_window_category = as.factor(trimws(as.character(time_window_category))),
  case_definition = trimws(as.character(case_definition))
)]

# Drop truly empty time_window_category levels if present
if ("time_window_category" %in% names(data)) {
  if (any(nchar(levels(data$time_window_category)) == 0)) {
    data$time_window_category <- forcats::fct_drop(data$time_window_category)
  }
}

# Create the results directory if it doesn't exist
dir.create(results_output_path, showWarnings = FALSE, recursive = TRUE)

# --- GLOBAL KNOTS (computed once, reused everywhere) ---
age_all <- data$Age[is.finite(data$Age)]
BKN <- range(age_all, na.rm = TRUE)
IK2 <- as.numeric(quantile(age_all, 0.50, na.rm = TRUE))
message(sprintf("Global spline knots set. Boundary: [%.2f, %.2f]; Internal(df=2): %.2f",
                BKN[1], BKN[2], IK2))

# --- Step 2: Define the Analysis, Plotting, and Saving Function ---
generate_prevalence_results <- function(pain_site, full_data, results_path) {
  message(paste("\n--- Processing Pain Site:", pain_site, "---"))
  
  tryCatch({
    
    # Basic filtering
    analysis_data <- full_data[
      !is.na(get(pain_site)) & !is.na(Age) &
        !is.na(Sex) & !is.na(Global_Region_Grouped) & !is.na(Study)
    ]
    
    # Region exclusions for specific sites
    if (pain_site %in% la_missing_sites) {
      analysis_data <- analysis_data[
        !Global_Region_Grouped %in% c("Latin America", "North Africa and West Asia")
      ]
    }
    
    # Drop unused factor levels
    analysis_data[, `:=`(
      Global_Region_Grouped = droplevels(Global_Region_Grouped),
      Study = droplevels(Study),
      Sex = droplevels(Sex),
      time_window_category = droplevels(time_window_category)
    )]
    
    # Per-site override for time_window_category
    win_col <- paste0("time_window_category_", pain_site)
    if (win_col %in% names(analysis_data)) {
      analysis_data[, time_window_category := ifelse(
        !is.na(get(win_col)),
        as.character(get(win_col)),
        as.character(time_window_category)
      )]
    } else {
      warning("Missing per-site window column: ", win_col, " (using study default).")
    }
    
    analysis_data[, time_window_category := factor(
      trimws(as.character(time_window_category)),
      levels = c("past_week", "past_month", "past_3_plus_months", "chronic_explicit")
    )]
    
    # ---------------- case_definition (mirrors IPD logic) ----------------
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
    } else {
      warning("Missing per-site case_definition column: ", case_col, " (using global default).")
    }
    
    analysis_data[case_definition == "", case_definition := NA_character_]
    analysis_data[!(case_definition %in% CASE_LEVELS), case_definition := NA_character_]
    
    # Drop missing covariates required for adjustment/standardization
    analysis_data <- analysis_data[!is.na(time_window_category)]
    analysis_data <- analysis_data[!is.na(case_definition)]
    analysis_data[, case_definition := factor(case_definition)]
    if ("general" %in% levels(analysis_data$case_definition)) {
      analysis_data[, case_definition := relevel(case_definition, ref = "general")]
    }
    # --------------------------------------------------------------------
    
    if (nrow(analysis_data) == 0) {
      warning(paste("Skipping", pain_site, "- No data available after filtering."))
      return(NULL)
    }
    
    # Ensure response is numeric 0/1-ish
    if (!is.numeric(analysis_data[[pain_site]])) {
      analysis_data[[pain_site]] <- as.numeric(as.character(analysis_data[[pain_site]]))
    }
    
    # Model with case_definition adjustment
    model_formula <- as.formula(paste0(
      pain_site, " ~ ",
      "ns(Age, df = 2) * Global_Region_Grouped + ",
      "ns(Age, df = 2) *  Sex +",  
      "time_window_category + case_definition + ",
      "(1 | Study) + (1 | Global_Region_Grouped:Study)"
    ))
    
    model <- glmer(
      model_formula,
      data = analysis_data,
      family = binomial(link = "logit"),
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
    )
    
    if (isSingular(model, tol = 1e-4)) {
      message(">> Warning: model is near-singular. Consider simplifying interactions or df if needed.")
    }
    
    # Age grid within boundaries
    min_age_analysis <- max(min(analysis_data$Age, na.rm = TRUE), BKN[1])
    max_age_analysis <- min(max(analysis_data$Age, na.rm = TRUE), BKN[2])
    if (!is.finite(min_age_analysis) || !is.finite(max_age_analysis) || min_age_analysis >= max_age_analysis) {
      warning("Age range invalid after applying boundary knots; skipping.")
      return(NULL)
    }
    age_grid_orig <- seq(min_age_analysis, max_age_analysis, length.out = 100)
    
    # Standardization
    sex_levels  <- levels(analysis_data$Sex)
    reg_levels  <- levels(analysis_data$Global_Region_Grouped)
    win_levels  <- levels(analysis_data$time_window_category)
    case_levels <- levels(analysis_data$case_definition)
    
    # Window standardization to past_month when available
    target_window <- "past_month"
    if (!(target_window %in% win_levels)) {
      tmpw <- analysis_data[, .N, by = time_window_category][order(-N)]
      target_window <- as.character(tmpw$time_window_category[1])
      warning(sprintf("'%s' not present; using '%s' for %s.",
                      "past_month", target_window, pain_site))
    }
    
    # Case-definition standardization to general when available
    target_case <- PRED_CASEDEF_DEFAULT
    if (!(target_case %in% case_levels)) {
      tmpc <- analysis_data[, .N, by = case_definition][order(-N)]
      target_case <- as.character(tmpc$case_definition[1])
      warning(sprintf("'%s' not present; using '%s' for %s.",
                      PRED_CASEDEF_DEFAULT, target_case, pain_site))
    }
    
    # Equal sex weights
    sex_w <- setNames(rep(1 / length(sex_levels), length(sex_levels)), sex_levels)
    
    # Region-specific grid (standardized to one window and one case_definition)
    region_grid <- CJ(
      Age = age_grid_orig,
      Global_Region_Grouped = reg_levels,
      Sex = sex_levels,
      time_window_category = factor(target_window, levels = win_levels),
      case_definition = factor(target_case, levels = case_levels),
      sorted = FALSE
    )
    region_grid[, w := sex_w[Sex]]
    
    regional_std <- avg_predictions(
      model,
      newdata = region_grid,
      by = c("Age", "Global_Region_Grouped"),
      type = "response",
      re.form = NA,
      wts = region_grid$w
    )
    regional_results_df <- as.data.frame(regional_std) |>
      dplyr::rename(prob = estimate, asymp.LCL = conf.low, asymp.UCL = conf.high, SE = std.error)
    regional_results_df$df <- NA_real_
    
    # Global weights (region mix from analysis_data)
    tmp_reg <- analysis_data[, .N, by = Global_Region_Grouped]
    tmp_reg[, w := N / sum(N)]
    REG_W <- setNames(tmp_reg$w, as.character(tmp_reg$Global_Region_Grouped))
    REG_W <- REG_W[reg_levels]
    REG_W <- REG_W / sum(REG_W)
    
    global_grid <- CJ(
      Age = age_grid_orig,
      Global_Region_Grouped = reg_levels,
      Sex = sex_levels,
      time_window_category = factor(target_window, levels = win_levels),
      case_definition = factor(target_case, levels = case_levels),
      sorted = FALSE
    )
    global_grid[, w := REG_W[as.character(Global_Region_Grouped)] * sex_w[Sex]]
    
    global_std <- avg_predictions(
      model,
      newdata = global_grid,
      by = "Age",
      type = "response",
      re.form = NA,
      wts = global_grid$w
    )
    global_results_df <- as.data.frame(global_std) |>
      dplyr::rename(prob = estimate, asymp.LCL = conf.low, asymp.UCL = conf.high, SE = std.error)
    global_results_df$df <- NA_real_
    global_results_df$Global_Region_Grouped <- "Global Average"
    
    # Combine standardized regional + global
    combined_results <- rbind(regional_results_df, global_results_df)
    setDT(combined_results)
    
    results_to_save <- combined_results[, .(
      Pain_Site = pain_site,
      Region = Global_Region_Grouped,
      Age = round(Age, 2),
      Predicted_Probability = prob,
      Lower_CI = asymp.LCL,
      Upper_CI = asymp.UCL,
      SE = SE,
      df = df
    )]
    
    # ===== Annual Percent Change (APC) =====
    # Prevalence for joins
    p_region <- as.data.frame(regional_std) |>
      dplyr::rename(prob = estimate) |>
      dplyr::select(Age, Global_Region_Grouped, prob)
    
    p_global <- as.data.frame(global_std) |>
      dplyr::rename(prob = estimate) |>
      dplyr::mutate(Global_Region_Grouped = "Global Average") |>
      dplyr::select(Age, Global_Region_Grouped, prob)
    
    # Derivative on probability scale
    d_region <- slopes(
      model, variables = "Age",
      newdata = region_grid,
      by = c("Age", "Global_Region_Grouped"),
      re.form = NA, wts = region_grid$w,
      type = "response"
    ) |>
      as.data.frame() |>
      dplyr::rename(dp = estimate, dp_L = conf.low, dp_U = conf.high, dp_SE = std.error)
    
    d_global <- slopes(
      model, variables = "Age",
      newdata = global_grid,
      by = "Age",
      re.form = NA, wts = global_grid$w,
      type = "response"
    ) |>
      as.data.frame() |>
      dplyr::rename(dp = estimate, dp_L = conf.low, dp_U = conf.high, dp_SE = std.error) |>
      dplyr::mutate(Global_Region_Grouped = "Global Average")
    
    # APC = 100 * dp/p; set to NA when p <= 0
    apc_region <- d_region |>
      dplyr::left_join(p_region, by = c("Age", "Global_Region_Grouped")) |>
      dplyr::mutate(
        den = ifelse(prob > 0, prob, NA_real_),
        APC = 100 * dp / den,
        APC_LCL = 100 * dp_L / den,
        APC_UCL = 100 * dp_U / den
      ) |>
      dplyr::transmute(
        Pain_Site = pain_site,
        Region = Global_Region_Grouped,
        Age,
        PercentChangePerYear = APC,
        Lower_CI = APC_LCL,
        Upper_CI = APC_UCL,
        SE = dp_SE,
        df = NA_real_
      )
    
    apc_global <- d_global |>
      dplyr::left_join(p_global, by = c("Age", "Global_Region_Grouped")) |>
      dplyr::mutate(
        den = ifelse(prob > 0, prob, NA_real_),
        APC = 100 * dp / den,
        APC_LCL = 100 * dp_L / den,
        APC_UCL = 100 * dp_U / den
      ) |>
      dplyr::transmute(
        Pain_Site = pain_site,
        Region = Global_Region_Grouped,
        Age,
        PercentChangePerYear = APC,
        Lower_CI = APC_LCL,
        Upper_CI = APC_UCL,
        SE = dp_SE,
        df = NA_real_
      )
    
    rates_to_save <- dplyr::bind_rows(apc_region, apc_global)
    
    # Save CSV per site (prevalence)
    output_filename <- file.path(results_path, paste0("Predicted_Prevalence_", pain_site, ".csv"))
    fwrite(results_to_save, file = output_filename)
    message(paste("--> Results for", pain_site, "saved to", output_filename))
    
    # Plot
    regional_plot_data <- results_to_save[Region %in% main_regions_to_plot]
    global_plot_data <- results_to_save[Region == "Global Average"]
    
    current_regions_in_plot <- unique(regional_plot_data$Region)
    plot_colors <- region_colors[c(as.character(current_regions_in_plot), "Global Average")]
    
    linetype_values <- c(
      setNames(rep("solid", length(current_regions_in_plot)), current_regions_in_plot),
      "Global Average" = "dashed"
    )
    
    prevalence_plot <- ggplot() +
      geom_ribbon(
        data = regional_plot_data,
        aes(x = Age, ymin = Lower_CI, ymax = Upper_CI, fill = Region),
        alpha = 0.05, linetype = 0
      ) +
      geom_line(
        data = regional_plot_data,
        aes(x = Age, y = Predicted_Probability, color = Region, linetype = Region),
        linewidth = 1.2
      ) +
      geom_ribbon(
        data = global_plot_data,
        aes(x = Age, ymin = Lower_CI, ymax = Upper_CI),
        fill = "black", alpha = 0.05, inherit.aes = FALSE
      ) +
      geom_line(
        data = global_plot_data,
        aes(x = Age, y = Predicted_Probability, linetype = "Global Average"),
        color = "black", linewidth = 1.0, inherit.aes = FALSE
      ) +
      scale_y_continuous(
        name = "Predicted Prevalence of Pain",
        labels = scales::percent_format(accuracy = 1),
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.05))
      ) +
      scale_x_continuous(
        name = "Age (Years)",
        breaks = seq(round(min_age_analysis / 10) * 10, round(max_age_analysis / 10) * 10, by = 10)
      ) +
      scale_color_manual(values = plot_colors, name = "Region") +
      scale_fill_manual(values = plot_colors, name = "Region", guide = "none") +
      scale_linetype_manual(values = linetype_values, name = "Region") +
      labs(
        title = paste("Predicted Lifespan Prevalence of", pain_site),
        subtitle = paste0(
          "Past-month reference; Sex averaged equally; Global curve uses region data-mix weights; ",
          "Case-definition standardized to: ", target_case
        ),
        color = "Region", linetype = "Region"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        legend.box = "horizontal"
      )
    
    return(list(plot = prevalence_plot, data = results_to_save, rate_data = rates_to_save))
    
  }, error = function(e) {
    message(paste("<<<< ERROR processing", pain_site, ">>>>"))
    message(e$message)
    return(NULL)
  })
}

# --- Step 3: Execute Analysis Serially using a for-loop ---
message("\n--- Starting serial analysis for all specified pain sites ---")
all_results <- vector("list", length(pain_sites_to_process))
names(all_results) <- pain_sites_to_process

for (site in pain_sites_to_process) {
  all_results[[site]] <- generate_prevalence_results(
    pain_site = site,
    full_data = data,
    results_path = results_output_path
  )
}

# --- Step 4: Process, Display, and Save All Results ---
message("\n--- All analyses complete. Processing and displaying results... ---")
successful <- vapply(all_results, function(x) !is.null(x), logical(1))

if (any(successful)) {
  for (site_name in names(all_results)[successful]) {
    message(paste0("\nDisplaying plot for: ", site_name))
    print(all_results[[site_name]]$plot)
  }
  
  # Master prevalence
  all_data_to_save <- rbindlist(lapply(all_results[successful], `[[`, "data"), use.names = TRUE, fill = TRUE)
  master_csv_path <- file.path(results_output_path, "Predicted_Prevalence_ALL_SITES.csv")
  fwrite(all_data_to_save, file = master_csv_path)
  message(paste0("\n--- Master CSV with all results saved to: ", master_csv_path, " ---"))
  
  # Master APC (Rate of Change)
  all_rates_to_save <- rbindlist(lapply(all_results[successful], `[[`, "rate_data"),
                                 use.names = TRUE, fill = TRUE)
  rate_csv_path <- file.path(results_output_path, "RateOfChange_ALL_SITES.csv")
  fwrite(all_rates_to_save, file = rate_csv_path)
  message(paste0("\n--- Master APC CSV saved to: ", rate_csv_path, " ---"))
  
  status_dt <- data.table(
    Pain_Site = names(all_results),
    Fit_Success = successful
  )
  fwrite(status_dt, file = file.path(results_output_path, "Model_Fit_Status.csv"))
} else {
  message("No results were generated successfully.")
}

message("\n--- Script finished. ---")
