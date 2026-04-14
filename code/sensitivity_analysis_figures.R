# Load libraries
library(tidyverse)
library(GenBinomApps) # to get Clopper-Pearson confidence interval function
library(cols4all)
library(future)
library(doFuture)
library(foreach)
library(progressr) # progress bars for long parallel computations
library(cowplot)
library(latex2exp)
library(fitdistrplus)
library(data.table)
library(multidplyr)

# Labeling functions
appender_R0 <- function(string) TeX(paste("$R_0 = $", string))
appender_sigma <- function(string) TeX(paste("$\\sigma = $", string))
appender_Thv <- function(string) TeX(paste("$\\tau_{hv} = $", string))
appender_var <- function(var, val) TeX(paste(var, val))

# # Parameters ---- 
b_baseline <- 0.3 # biting rate
Tvh_baseline <- 0.5 # vector transmission-probability susceptible vector being infected after contact with infectious host
Nh_baseline <- 10000 # number of hosts in system
Nv_baseline <- 100000 # number of vectors in system
muv <- 0.1 # vector mortality rate
gammah <- 0.1 # host recovery rate
# 
R0_from_Thv_function <- function(Thv, b, Tvh, Nh, Nv) {
  sqrt(Thv * b^2 * Tvh * Nv / (Nh * gammah * muv))
}

# Sigma values (environmental noise level in 0-0.3) to test
sigmas <- seq(0, 1, by = 0.05)

# Final time point
max_time = 3650

# Load data from Julia output ----
# Sensitivity analysis data
sens_df = read_csv("./data/collect_outputs_sensitivity.csv") 

# Summarize key statistics ----

sens_summary_df <- sens_df |>
  # rename quantities to match plots
  mutate(name = case_when(
    name == "max_value" ~ "max_cases",
    # name == "positive_duration" ~ "duration",
    name == "exceeded_10" ~ "small_outbreak",
    name == "exceeded_100" ~ "big_outbreak",
    name == "positive_at_final" ~ "endemic",
    TRUE ~ name
  )) |> 
  mutate(R0 = round(R0_from_Thv_function(Thv, b, Tvh, Nh, Nv), 3)) |> 
  mutate(across(b:sigma, ~ round(.x, 3))) 


# Process statistics ----
# Calculate confidence intervals
# For both noise types

sens_stats_df <- sens_summary_df |> 
  pivot_wider(names_from = statistic) |> 
  rowwise() |>
  mutate(
    max_val = if_else(name == "max_cases", Nh, max_time / 365),
    lower_ci = max(0, mean - 0.674 * sqrt(variance)),
    upper_ci = if_else(name %in% c("small_outbreak", "big_outbreak", "endemic", "zero_cases"), 
                       min(1, mean + 0.674 * sqrt(variance)),
                       min(max_val, mean + 0.674 * sqrt(variance))
    )
  ) |>
  dplyr::select(-max_val) |>
  ungroup()

sens_stats_df$R0_factor = factor(round(sens_stats_df$R0,3), levels = (unique(round(sens_stats_df$R0,3))))
write_rds(sens_stats_df, "./data/sens_stats.rds")

# Figure S2: Sensitivity analysis - curves --------------------------------
R0_colors = rev(c4a("kovesi.bu_bk_rd", 7))

# For each varied parameter, create one line plot for each value of that parameter
# varied parameters = b, Tvh, Nv, Nh

# Write the code for b then modify
b_line_plots_df <- sens_stats_df |> 
  filter(
    Tvh == Tvh_baseline,
    Nv == Nv_baseline,
    Nh == Nh_baseline
  ) |> unique()

sens_lineplot_func <- function(varied_parameter, chosen_output) {
  
  top_label <- switch(
    varied_parameter,
    "b"   = bquote("Biting rate (" * b * ")"),
    "Tvh" = bquote("Host susceptibility (" * tau[vh] * ")"),
    "Nh"  = bquote("Host population size (" * N[h] * ")"),
    "Nv"  = bquote("Vector population size (" * N[v] * ")"),
    expression("")
  )
  
  appender_top = function(string) {TeX(appender_var(top_label, string))}
  
  
  name_label = case_when(
    chosen_output == "big_outbreak" ~ "$Probability of an outbreak of over 100 cases$",
    chosen_output == "max_cases" ~ "Peak number of cases",
    chosen_output == "duration" ~ "Duration (final time with at least one infection)"
  )
  
  get_varying_df <- function(df, varied_parameter) {
    # Define base values
    bases <- list(b = b_baseline, Tvh = Tvh_baseline, Nh = Nh_baseline, Nv = Nv_baseline)
    
    # Remove the parameter we want to vary
    bases[varied_parameter] <- NULL
    
    # Create filter conditions dynamically
    conditions <- map2(names(bases), bases, 
                       ~ expr(!!sym(.x) == !!.y))
    
    # Apply all conditions
    df %>%
      dplyr::filter(!!!conditions)
  }
  
  get_label_func <- function(output_type) {
    switch(
      output_type,
      big_outbreak = scales::percent_format(accuracy = 1),
      scales::label_number(accuracy = 2)  # default fallback
    )
  }
  
  # Filter down to only the varied parameter
  line_plot_df <- sens_stats_df |> 
    mutate(round_sigma = round(sigma, 3)) |>
    filter(
      name %in% c("big_outbreak", "max_cases", "duration"),
      # R0_factor %in% c(0.95, 1.05, 1, 2, 3, 4, 5),
      # round_sigma %in% round(seq(0.0, 2.0, by = 0.05), 3)
    ) |> 
    get_varying_df(varied_parameter) |> 
    unique()
  
  var_breaks = unique(dplyr::select(line_plot_df, !!sym(varied_parameter))) |> pull()
  
  # Probability of a large outbreak
  line_plot <- line_plot_df |> 
    arrange(R0) |> 
    filter(name == chosen_output) |> 
    ggplot(aes(x = sigma)) +
    geom_line(aes(y = mean, color = factor(!!sym(varied_parameter)), group = factor(!!sym(varied_parameter))),
              lwd = 1) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = factor(!!sym(varied_parameter)), group = factor(!!sym(varied_parameter))), alpha = 0.25) +
    scale_linetype_manual( 
      values = c(1, 3)
    ) +
    scale_color_manual(
      values = (c4a("seaborn.flare", 5)),
      breaks = var_breaks
    ) +
    scale_fill_manual(
      values = (c4a("seaborn.flare", 5)), # alt: "seaborn.crest"
      breaks = var_breaks
    )  +
    theme_cowplot() +
    ggtitle("") +
    scale_x_continuous(
      unname(TeX("Environmental noise $[\\sigma]$")),
      expand = c(0,0)
    ) +
    scale_y_continuous(
      name = name_label,
      labels = get_label_func(chosen_output),
      expand = c(0,0.005)
    ) +
    facet_grid(~ Thv,
               labeller = labeller(.cols = as_labeller(appender_Thv, default = label_parsed))) +
    theme_minimal_grid(10) +  
    guides(
      color = guide_legend(
        title = top_label,
        title.position = "top",
        label.position = "bottom",
        direction = "horizontal",
        nrow = 1,
        override.aes = list(linewidth = 3),
      ),
      fill = guide_none(),
      linetype = guide_none()
    ) +
    theme(
      plot.title = element_text(size = 10),
      legend.text = element_text(size = 7, margin = margin(t = -1)),
      legend.spacing = unit(0, "pt"),
      legend.title = element_text(size = 8, margin = margin(b = -2)),
      legend.position = "top",
      legend.justification = "left",
      legend.direction = "horizontal",  
      legend.spacing.y = unit(0, "pt"),
      legend.key.width = unit(14, "pt"),
      plot.title.position = "panel",
      panel.spacing = unit(5, "mm"),
      legend.box.margin = margin(t = -18, b = -10)
    )
  line_plot
  
}

(sens_lineplot_func("b", "big_outbreak"))


sens_boxwhisker_func <- function(varied_parameter, chosen_output) {
  
  top_label <- switch(
    varied_parameter,
    "b"   = bquote("Biting rate (" * b * ")"),
    "Tvh" = bquote("Host susceptibility (" * tau[vh] * ")"),
    "Nh"  = bquote("Host population size (" * N[h] * ")"),
    "Nv"  = bquote("Vector population size (" * N[v] * ")"),
    expression("")
  )

  appender_top = function(string) {TeX(appender_var(top_label, string))}
  
  name_label = case_when(
    chosen_output == "big_outbreak" ~ "Probability of an outbreak of over 100 cases",
    chosen_output == "max_cases" ~ "Peak number of cases",
    chosen_output == "duration" ~ "Duration (final time with at least one infection)"
  )
  
  get_varying_df <- function(df, varied_parameter) {
    # Define base values
    bases <- list(b = b_baseline, Tvh = Tvh_baseline, Nh = Nh_baseline, Nv = Nv_baseline)
    
    # Remove the parameter we want to vary
    bases[varied_parameter] <- NULL
    
    # Create filter conditions dynamically
    conditions <- map2(names(bases), bases, 
                       ~ expr(!!sym(.x) == !!.y))
    
    # Apply all conditions
    df %>%
      dplyr::filter(!!!conditions)
  }
  
  get_label_func <- function(output_type) {
    switch(
      output_type,
      big_outbreak = scales::percent_format(accuracy = 1),
      scales::label_number(accuracy = 2)  # default fallback
    )
  }
  
  # Filter down to only the varied parameter
  line_plot_df <- sens_stats_df |> 
    mutate(round_sigma = round(sigma, 3)) |>
    filter(
      name %in% c("big_outbreak", "max_cases", "duration"),
      # R0_factor %in% c(0.95, 1.05, 1, 2, 3, 4, 5),
      # round_sigma %in% round(seq(0.0, 2.0, by = 0.05), 3)
    ) |> 
    get_varying_df(varied_parameter) |> 
    mutate(
      varied_value = factor(
        .data[[varied_parameter]],
        levels = sort(unique(.data[[varied_parameter]]))
      ),
      varied_row = factor(
        paste0(
          switch(
            varied_parameter,
            "b"   = "b",
            "Tvh" = "tau[vh]",
            "Nh"  = "N[h]",
            "Nv"  = "N[v]"
          ),
          "==",
          as.character(varied_value)
        ),
        levels = paste0(
          switch(
            varied_parameter,
            "b"   = "b",
            "Tvh" = "tau[vh]",
            "Nh"  = "N[h]",
            "Nv"  = "N[v]"
          ),
          "==",
          levels(varied_value)
        )
      )
    ) |> 
    unique()
  
  var_breaks = unique(dplyr::select(line_plot_df, !!sym(varied_parameter))) |> pull()
  
  # Probability of a large outbreak
  line_plot <- line_plot_df |> 
    arrange(R0) |> 
    filter(name == chosen_output) |> 
    mutate(
      mean     = if_else(rep(chosen_output == "max_cases", n()), pmin(mean, Nh), mean),
      lower_ci = if_else(rep(chosen_output == "max_cases", n()), pmin(lower_ci, Nh), lower_ci),
      upper_ci = if_else(rep(chosen_output == "max_cases", n()), pmin(upper_ci, Nh), upper_ci),
      min      = if_else(rep(chosen_output == "max_cases", n()), pmin(min, Nh), min),
      max      = if_else(rep(chosen_output == "max_cases", n()), pmin(max, Nh), max)
    ) |>
    ggplot(aes(x = sigma)) +
    geom_boxplot(
      aes(
        ymin = lower_ci,
        lower = lower_ci,
        middle = mean,
        upper = upper_ci,
        ymax = upper_ci,
        fill  = varied_value,
        group = interaction(sigma, varied_value, drop = TRUE)
      ),
      stat = "identity"
    ) +
    geom_line(aes(y = mean, group = factor(!!sym(varied_parameter))),
              lwd = 1) +
    scale_color_manual(
      values = (c4a("seaborn.flare", 5)),
      breaks = var_breaks
    ) +
    scale_fill_manual(
      values = (c4a("seaborn.flare", 5)), # alt: "seaborn.crest"
      breaks = var_breaks
    )  +
    theme_cowplot() +
    ggtitle("") +
    scale_x_continuous(
      unname(TeX("Environmental noise $[\\sigma]$")),
      expand = c(0,0)
    ) +
    scale_y_continuous(
      name = "",
      labels = get_label_func(chosen_output),
      expand = c(0, 0.005),
      limits = c(0,NA)
    ) +
    facet_grid(varied_row ~ Thv,
               switch = "y",
               scales = "free_y",
               labeller = labeller(
                 varied_row = label_parsed,
                 .cols = as_labeller(appender_Thv, default = label_parsed))) +
    theme_minimal_grid(8) +  
    guides(
      color = guide_legend(
        title = top_label,
        title.position = "top",
        label.position = "bottom",
        direction = "horizontal",
        nrow = 1,
        override.aes = list(linewidth = 3),
      ),
      fill = guide_none(),
      linetype = guide_none()
    ) +
    theme(
      plot.title = element_text(size = 10),
      legend.text = element_text(size = 7, margin = margin(t = -1)),
      legend.spacing = unit(0, "pt"),
      legend.title = element_text(size = 8, margin = margin(b = -2)),
      legend.position = "top",
      legend.justification = "left",
      legend.direction = "horizontal",  
      legend.spacing.y = unit(0, "pt"),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 90),
      legend.key.width = unit(14, "pt"),
      plot.title.position = "panel",
      panel.spacing = unit(5, "mm"),
      legend.box.margin = margin(t = -18, b = -10)
    ) +
    ggtitle(name_label)
  line_plot
  
}

p <- sens_boxwhisker_func("Tvh", "duration")

ggsave(paste0("./figures/Boxwhisker_test.png"),
       p, width = 9, height = 4, units = "in", dpi = 600)


varied_parameters <- c("b", "Tvh", "Nh", "Nv")
chosen_outputs <- c("big_outbreak", "max_cases", "duration")

for (vp in varied_parameters) {
  for (co in chosen_outputs) {
    
    p <- sens_boxwhisker_func(vp, co)
    
    file_stub <- paste0("./figures/Boxwhisker_", vp, "_", co)
    
    # ggsave(paste0(file_stub, ".png"),
    #        p, width = 9, height = 4, units = "in", dpi = 1200)
    
    ggsave(paste0(file_stub, ".eps"),
           p, width = 9, height = 4, units = "in", dpi = 600)
    
    # ggsave(paste0(file_stub, ".tif"),
    #        p, width = 9, height = 4, units = "in", dpi = 600)
  }
}