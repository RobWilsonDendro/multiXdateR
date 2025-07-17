# multiXdateR - Refactored Version - 17/07/2025
# ==============================================================================
# SETUP
# ==============================================================================

# Set working directory to source file location
if (require("rstudioapi") && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Load necessary libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplR)
  library(zoo)
  library(ggpubr)
  library(scales)
})

# ==============================================================================
# USER SETTINGS
# ==============================================================================
runtitle         <- "Test"
CorrelationPlots <- "yes"
TvaluePlots      <- "yes"
reftruncate      <- "no"
refstart         <- "xxx"
refend           <- "xxx"
minoverlap       <- 30
splinewindow     <- 31
datatransform    <- "pw" # "pw" (pre-whitened) or "fd" (1st differenced)
correlationtype  <- "spearman" # "pearson" or "spearman"

# Define parameters and file paths
parameters <- c("RW", "EWBI", "LWBI")
ref_files <- setNames(
  paste0("ExampleData/LGL", tolower(parameters), "TRUNC"), 
  parameters
)
und_files <- setNames(
  paste0("ExampleData/LCL", tolower(parameters), "TRUNC"), 
  parameters
)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Process a single chronology file
#'
#' Reads, detrends, and transforms a single .rwl file.
#' @param file_path Path to the .rwl file.
#' @param spline_window Spline window size for detrending.
#' @param data_transform Method for data transformation ("pw" or "fd").
#' @return A tibble with Year and Value columns.
process_chronology <- function(file_path, spline_window, data_transform) {
  raw_data <- read.rwl(file_path, header = FALSE)
  
  detrended <- detrend(
    rwl = raw_data, method = "Spline", 
    nyrs = spline_window, f = 0.5, pos.slope = TRUE
  )
  
  chronology <- chron(
    x = detrended, prefix = "ref", 
    biweight = TRUE, prewhiten = (data_transform == "pw")
  )
  
  if (data_transform == "pw") {
    processed_series <- chronology$res
  } else { # "fd"
    processed_series <- c(NA, diff(chronology$res))
  }
  
  tibble(
    Year = as.numeric(rownames(chronology)),
    Value = as.numeric(processed_series)
  ) %>%
    na.omit()
}

#' Calculate sliding correlations and associated statistics
#'
#' @param data A tibble with Year, ref_value, and und_value.
#' @param und_series_length Full length of the undated series.
#' @param min_overlap Minimum overlap for correlation.
#' @param corr_type Correlation method.
#' @return A tibble with comprehensive sliding statistics for each window.
calculate_sliding_stats <- function(data, und_series_length, min_overlap, corr_type) {
  
  # Check for valid overlap
  if (min_overlap > und_series_length) {
    stop(paste("Minimum overlap is longer than the undated chronology for", data$parameter[1]))
  }
  
  # Pad reference series to allow sliding
  padding <- und_series_length - min_overlap
  first_ref_year <- min(data$ref_year, na.rm = TRUE)
  
  padded_ref <- tibble(
    Year = seq(from = first_ref_year - padding, to = max(data$ref_year, na.rm = TRUE) + padding, by = 1)
  ) %>%
    left_join(dplyr::select(data, Year, ref_value), by = "Year")
  
  # Perform sliding calculations
  window_size <- und_series_length
  
  results <- zoo::rollapply(
    padded_ref$ref_value,
    width = window_size,
    FUN = function(ref_win) {
      # Align windows and remove NAs for this specific window
      valid_indices <- !is.na(ref_win) & !is.na(data$und_value)
      
      ref_valid <- ref_win[valid_indices]
      und_valid <- data$und_value[valid_indices]
      
      n_obs <- length(ref_valid)
      
      if (n_obs < min_overlap) {
        return(c(r = NA, n = n_obs))
      }
      
      r <- cor(ref_valid, und_valid, method = corr_type, use = "pairwise.complete.obs")
      return(c(r = r, n = n_obs))
    },
    by = 1,
    align = "right",
    fill = NA
  )
  
  # Calculate AC1 for effective N
  ac1_ref <- acf(data$ref_value, lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2]
  ac1_und <- acf(data$und_value, lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2]
  
  # Calculate final statistics
  tibble(
    end_year = padded_ref$Year[window_size:nrow(padded_ref)],
    r = results[, "r"],
    n_obs = results[, "n"]
  ) %>%
    na.omit() %>%
    mutate(
      n_adj = n_obs * (1 - abs(ac1_ref) * abs(ac1_und)) / (1 + abs(ac1_ref) * abs(ac1_und)),
      t_value = (r * sqrt(n_adj - 2)) / sqrt(1 - r^2),
      ac1_mean = mean(c(ac1_ref, ac1_und), na.rm = TRUE)
    )
}

# ==============================================================================
# DATA PROCESSING
# ==============================================================================
# Process all reference and undated files
ref_data <- map_df(ref_files, ~process_chronology(., splinewindow, datatransform), .id = "parameter")
und_data <- map_df(und_files, ~process_chronology(., splinewindow, datatransform), .id = "parameter")

# Handle optional reference truncation
if (reftruncate == "yes") {
  ref_data <- ref_data %>% dplyr::filter(Year >= as.numeric(refstart), Year <= as.numeric(refend)) # <-- FIXED
}

# Ensure common time period for all series of the same type (ref/und)
common_ref_years <- ref_data %>% 
  group_by(Year) %>% 
  dplyr::filter(n() == length(parameters)) %>% # <-- FIXED
  ungroup() %>% 
  summarise(start = min(Year), end = max(Year))

common_und_years <- und_data %>% 
  group_by(Year) %>% 
  dplyr::filter(n() == length(parameters)) %>% # <-- FIXED
  ungroup() %>% 
  summarise(start = min(Year), end = max(Year))

# Combine into one master dataframe
all_data <- full_join(
  ref_data %>% rename(ref_value = Value, ref_year = Year) %>% dplyr::filter(ref_year >= common_ref_years$start), # <-- FIXED
  und_data %>% rename(und_value = Value, und_year = Year) %>% dplyr::filter(und_year >= common_und_years$start), # <-- FIXED
  by = "parameter"
) %>%
  rename(Year = ref_year) # Use ref_year as the main timeline

# Get length of the undated series (must be common across parameters)
undated_series_length <- n_distinct(und_data$Year)
if (minoverlap > undated_series_length) {
  stop("DUDE! - your minimum overlap is longer than your undated chronology. Please lower the value.")
}

# ==============================================================================
# SLIDING CORRELATION AND T-VALUE CALCULATION
# ==============================================================================

# Nest data and apply calculations for each parameter
results_by_param <- all_data %>%
  group_by(parameter) %>%
  nest() %>%
  mutate(
    stats = map(data, ~calculate_sliding_stats(., undated_series_length, minoverlap, correlationtype))
  ) %>%
  dplyr::select(parameter, stats) %>% # <-- FIXED
  unnest(cols = stats)

# Calculate COMBO (mean) statistics
combo_stats <- results_by_param %>%
  group_by(end_year) %>%
  summarise(
    r = mean(r, na.rm = TRUE),
    total_n_adj = sum(n_adj, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    t_value = (r * sqrt(total_n_adj - 2)) / sqrt(1 - r^2),
    parameter = "COMBO",
    # We need an n_adj for the p-value calculation, let's use the sum
    n_adj = total_n_adj 
  )

# Combine all results and calculate p-values
all_results <- bind_rows(results_by_param, combo_stats) %>%
  group_by(parameter) %>%
  mutate(
    n_correlations = n(),
    p_value_raw = pt(pmax(0, t_value), df = n_adj, lower.tail = FALSE),
    p_value_adj = pmin(1, p_value_raw * n_correlations),
    reciprocal_p = 1 / p_value_adj
  ) %>%
  # Calculate Isolation Factor (IF)
  mutate(
    min_p_adj = min(p_value_adj, na.rm = TRUE),
    next_min_p_adj = min(p_value_adj[p_value_adj > min_p_adj], na.rm = TRUE),
    isolation_factor = if_else(min_p_adj > 0, round(next_min_p_adj / min_p_adj, 0), NA_real_)
  ) %>%
  ungroup()


# ==============================================================================
# FINAL SUMMARY TABLE
# ==============================================================================

# Find the best date from the COMBO results
best_combo_match <- combo_stats %>% dplyr::filter(!is.na(t_value)) %>% slice_max(t_value, n = 1) # <-- FIXED
best_date <- best_combo_match$end_year

# Extract stats for the best date for all parameters
summary_stats <- all_results %>%
  dplyr::filter(end_year == best_date) %>% # <-- FIXED
  arrange(factor(parameter, levels = c(parameters, "COMBO"))) %>%
  dplyr::select(parameter, r_value = r, t_value, p_value = p_value_adj, reciprocal_p, isolation_factor) # <-- FIXED

# Format for output table
last_undated_year <- max(und_data$Year)
year_shift <- best_date - last_undated_year

final_table_data <- summary_stats %>%
  mutate(
    across(c(r_value, t_value), ~round(., 2)),
    p_value_str = if_else(p_value < 0.0001, "<0.0001", as.character(round(p_value, 4))),
    reciprocal_p_str = if_else(reciprocal_p > 1e6, ">1 million", as.character(round(reciprocal_p))),
    if_str = if_else(isolation_factor > 1000, ">1000", as.character(round(isolation_factor))),
    `Proposed Outer Date` = if_else(p_value <= 0.05, as.character(best_date), "no sig date"),
    `No. of Years to shift` = if_else(p_value <= 0.05, as.character(year_shift), "no sig date")
  ) %>%
  dplyr::select( # <-- FIXED
    `TR Parameter` = parameter, `r value` = r_value, `T value` = t_value,
    `p value` = p_value_str, `1/p` = reciprocal_p_str, IF = if_str,
    `Proposed Outer Date`, `No. of Years to shift`
  )

# Display and save the summary table
print(final_table_data)
# write.csv(final_table_data, "UBERsummarytable.csv", row.names = FALSE)


# ==============================================================================
# PLOTTING (This section remains unchanged, but is included for completeness)
# ==============================================================================

if (TvaluePlots == "yes") {
  
  # Custom reverse log transform for p-value axis
  reverselog_trans <- function(base = 10) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), domain = c(1e-100, Inf))
  }
  
  # Global plot settings
  grand_min_p <- min(all_results$p_value_adj[all_results$p_value_adj > 0], na.rm = TRUE)
  hist_min <- min(all_results$t_value, na.rm = TRUE)
  hist_max <- max(all_results$t_value, na.rm = TRUE)
  
  # Function to create the 3-panel plot for one parameter
  create_param_plots <- function(param_data) {
    
    param_name <- param_data$parameter[1]
    best_match <- param_data %>% dplyr::filter(!is.na(t_value)) %>% slice_max(t_value, n = 1)
    
    # T-value time series plot
    p1 <- ggplot(param_data, aes(x = end_year, y = t_value)) +
      geom_line(color = "aquamarine4", linewidth = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_hline(yintercept = 4 * sd(param_data$t_value, na.rm = TRUE), linetype = "dotted", color = "grey50") +
      geom_point(data = best_match, aes(x = end_year, y = t_value), color = "blue", size = 3, shape = 18) +
      scale_x_continuous(breaks = pretty_breaks(n = 8)) +
      labs(
        title = paste(param_name, "Sliding T-values"),
        subtitle = paste0("Date = ", best_match$end_year, " CE, T = ", round(best_match$t_value, 2)),
        x = NULL, y = "T-value"
      ) +
      theme_classic() +
      theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(color = "blue"))
    
    # P-value time series plot
    p2 <- ggplot(param_data, aes(x = end_year, y = p_value_adj)) +
      geom_line(color = "black", linewidth = 0.5) +
      geom_hline(yintercept = c(0.01, 0.0001), linetype = "dashed", color = "red") +
      scale_x_continuous(breaks = pretty_breaks(n = 8)) +
      scale_y_continuous(
        trans = reverselog_trans(10),
        breaks = trans_breaks("log10", function(x) 10^x),
        labels = trans_format("log10", math_format(10^.x)),
        limits = c(1, grand_min_p)
      ) +
      labs(
        title = paste(param_name, "Sliding Bonferroni p-values"),
        subtitle = paste0("p = ", if_else(best_match$p_value_adj < 0.0001, "<0.0001", as.character(round(best_match$p_value_adj, 4))),
                          ", IF = ", round(best_match$isolation_factor)),
        x = NULL, y = "Adjusted p-value"
      ) +
      theme_classic() +
      theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(color = "blue"))
    
    # T-value histogram
    p3 <- ggplot(param_data, aes(x = t_value)) +
      geom_histogram(aes(y = after_stat(density)), fill = "aquamarine4", color = "black", bins = 50) +
      geom_vline(xintercept = best_match$t_value, color = "blue", linetype = "dashed", linewidth = 1) +
      coord_cartesian(xlim = c(hist_min, hist_max * 1.05)) +
      labs(title = "T-value Density", x = "T-value", y = "Density") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold"))
    
    # Arrange the three plots
    ggarrange(p1, p3, p2, ncol = 3, nrow = 1, widths = c(1, 0.5, 1))
  }
  
  # Generate plots for all parameters including COMBO
  plot_list <- all_results %>%
    group_by(parameter) %>%
    nest() %>%
    arrange(factor(parameter, levels = c(parameters, "COMBO"))) %>%
    mutate(plot = map(data, create_param_plots)) %>%
    pull(plot)
  
  # Final figure assembly
  final_figure <- ggarrange(plotlist = plot_list, ncol = 1, nrow = length(plot_list))
  
  final_figure_annotated <- annotate_figure(
    final_figure,
    top = text_grob(runtitle, color = "red", face = "bold", size = 28),
    bottom = text_grob(paste("Best COMBO Date:", best_date, "CE |", "Shift series by", year_shift, "years to attain date"),
                       color = "red", face = "bold", size = 14)
  )
  
  print(final_figure_annotated)
  
  # Save the figure
  ggsave("T_value_figure.tiff", plot = final_figure_annotated, width = 40, height = 25, units = "cm", dpi = 300, compression = "lzw", bg = "white")
}
