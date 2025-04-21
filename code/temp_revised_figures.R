# Functions ----

# Generic plot function
generic_heat_function <- function(output_name, type_name) {
  
  # Set heatmap palette to match unit of output
  palette = case_when(
    output_name %in% c("small_outbreak", "big_outbreak", "endemic") ~ "viridis",
    output_name %in% c("duration", "peak_time", "duration_dieout", "duration_10", "duration_100") ~ "plasma",
    output_name %in% c("max_cases") ~ "turbo",
  )
  
  # Label legend appropriately
  leg_label = case_when(
    output_name %in% c("small_outbreak", "big_outbreak", "endemic") ~ "Probability",
    output_name %in% c("duration", "peak_time", "duration_dieout", "duration_10", "duration_100") ~ "Time\n[years]",
    output_name %in% c("max_cases") ~ "Number\nof cases",
  )
  
  # Define zero proxy values
  zero_proxy = case_when(
    output_name %in% c("small_outbreak", "big_outbreak", "endemic") ~ 1E-4,
    output_name %in% c("duration", "peak_time", "duration_dieout", "duration_10", "duration_100") ~ 1,
    output_name %in% c("max_cases") ~ 1,
  )
  
  
  
  value_labels <- ifelse(
    leg_label == "Probability",
    label_percent(),
    function (x) {x}
  )
  
  sigma_width = unique(diff(full_stats_df$sigma))[2]
  R0_height = unique(diff(full_stats_df$R0))[2]
  
  fixed_df = full_stats_df %>%
    filter(name == output_name) %>% 
    dplyr::select(sigma, R0, mean, type) %>% 
    filter(if_all(c(sigma, R0, mean), ~ is.finite(.x) & !is.na(.x)))
  
  contour_df = smooth_zero_contour(fixed_df, output_name, type_name, zero_proxy)
  
  out = full_stats_df %>% 
    ungroup() %>% 
    filter(name == output_name, type == type_name) %>%
    ggplot(aes(y = R0 - R0_height/2)) +
    # Plot values
    geom_tile(aes(
      x = sigma - sigma_width/2,
      fill = mean,
      width = sigma_width,
      height = R0_height
    )) +
    # Add smoothed contour for zero proxy values
    geom_hline(yintercept = 1, color = "red", lwd = 0.5) +
    geom_path(data = contour_df,
              aes(
                x = sigma,
                y = R0,
                group = group),
              color = "white", lwd = 0.5) +
    # Annotate x-axis for demographic noise
    scale_x_continuous(
      TeX("Environmental noise strength [$\\sigma$]"), 
      limits = c(0, NA),
      expand = c(0,0),
      breaks = seq(0, 2.0, by = 0.25),
      labels = seq(0, 2.0, by = 0.25)
    ) +
    scale_y_continuous(
      name = ifelse(type_name == "enviro",
                    TeX("                                            Basic reproduction number  [$R_0$]"),
                    "          "),
      limits = c(0, 5),
      expand = c(0,0)
    ) +
    # color:
    scale_fill_viridis_c(
      name = output_name,
      option = palette,
      limits = case_when(
        output_name %in% c("duration", "peak_time", "duration_dieout", "duration_10", "duration_100") ~ c(0,10),
        output_name %in% c("small_outbreak", "big_outbreak", "endemic") ~ c(0,1),
        output_name %in% c("max_cases") ~ c(0,10000)
      ),
      breaks = case_when(
        output_name %in% c("duration", "peak_time", "duration_dieout", "duration_10", "duration_100") ~ c(0, 5, 10),
        output_name %in% c("small_outbreak", "big_outbreak", "endemic") ~ c(0,0.5, 1),
        output_name %in% c("max_cases") ~ c(0,5000, 10000)
      ),
      labels = value_labels
    ) +
    # legend:
    guides(
      fill = guide_colourbar(
        title = leg_label,
        position = "right",
        direction = "vertical",
        title.position = "top",
        title.hjust = 0,
        title.vjust = 1,
        barheight = 8,
        show.limits = TRUE,
        draw.ulim = TRUE,
        draw.llim = TRUE,
      )
    ) +
    theme_half_open(10) +
    theme(
      legend.position = "right",
      legend.justification = "center",
      legend.box.just = "right"
    )
  
  return(out)
}


line_plots_df <- all_stats_df %>% 
  mutate(round_sigma = round(sigma, 3)) %>%
  filter(
    name %in% c("big_outbreak", "max_cases", "duration"),
    R0_factor %in% c(0.95, 1.05, 1, 2, 3, 4, 5),
    round_sigma %in% round(seq(0.0, 2.0, by = 0.05), 3)
  ) %>% 
  mutate(label = case_when(name == "big_outbreak" ~ "A. Probability",
                           name == "max_cases" ~ "B. Intensity",
                           name == "duration" ~ "C. Duration"))

# New Figure 3: Probability plots ----

## Subplot A) Line plot
#  = full-size, on left
Big_outbreak_mean <- line_plots_df %>% 
  arrange(R0) %>% 
  filter(name == "big_outbreak") %>% 
  ggplot(aes(x = sigma)) +
  geom_line(aes(y = mean, color = R0_factor, group = R0_factor),
            lwd = 1) +
  scale_color_manual(values = R0_colors) +
  scale_fill_manual(values = R0_colors)  +
  theme_cowplot() +
  labs(color = unname(TeX("$R_0$")), fill = unname(TeX("$R_0$"))) +
  ggtitle("A.") +
  scale_x_continuous(
    unname(TeX("Environmental noise strength $[\\sigma]$")),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    "Probability",
    labels = scales::percent_format(accuracy = 1),
    expand = c(0,0.005)
  ) +
  theme_minimal_grid(10) +
  guides(
    color = guide_legend(
      title = TeX("Basic reproduction number $[R_0]$:"),
      title.position = "top",
      label.position = "bottom",
      reverse = TRUE,
      direction = "horizontal",
      nrow = 1
    )
  ) +
  theme(
    plot.title = element_text(size = 10),
    legend.text = element_text(size = 7),
    legend.spacing = unit(0, "pt"),
    legend.title = element_text(size = 8),
    legend.position = "top",
    legend.justification = "right",
    legend.direction = "horizontal",  
    legend.spacing.y = unit(0, "pt"),
    legend.key.width = unit(14, "pt"),
    legend.key.height = unit(1, "pt"),
    plot.title.position = "panel",
    legend.margin = margin(t = -17, b = -5)
  )

## Subplot B) Heatmap for full model
#  = half-size, on top-right
#  * label clearly with subtitle
#  * center colorbar with below figure
Big_outbreak_heat_temp <- generic_heat_function("big_outbreak", "all")

legend_temp = get_legend(Big_outbreak_heat_temp)

Big_outbreak_heat = Big_outbreak_heat_temp +
  labs(
    title = "B. Full model"
  ) +
  theme(
    plot.title = element_text(size = 10),
    axis.title.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),  # Remove x-axis line
    axis.title.y = element_blank()
  ) +
  guides(fill = "none") #+
  # theme(legend.margin=margin(t=-0.5,l=0.0,b=-0.275,r=0.25, unit='cm'))

## Subplot C) Heatmap for environmental noise sub-model
#  = half-size, on bottom-right
#  * label clearly with subtitle
Big_outbreak_heat_enviro <- generic_heat_function("big_outbreak", "enviro") +
  labs(
    title = "C. Environmental noise submodel"
  ) +
  theme(
    plot.title = element_text(size = 10),
    # axis.title.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    # axis.title.y = element_blank()
  ) +
  guides(fill = "none") #+
  # theme(legend.margin=margin(t=-0.5,l=0.0,b=-0.275,r=0.25, unit='cm'))

Figure3 = ggpubr::ggarrange(Big_outbreak_mean,
                            egg::ggarrange(
                              Big_outbreak_heat, 
                              Big_outbreak_heat_enviro, 
                              ncol = 1,
                              heights = c(0.475,0.525)
                            ),
                            legend_temp,
                            ncol = 3,
                            widths = c(0.4, 0.525, 0.075)
)

ggsave("./figures/New_Figure3.png", Figure3, width = 9, height = 4, units = "in", dpi = 1200)

# New Figure 4: Intensity plots ----

## Subplot A) Line plot
#  = full-size, on left
Peak_cases_mean <- line_plots_df %>% 
  arrange(R0) %>% 
  filter(name == "max_cases") %>% 
  ggplot(aes(x = sigma)) +
  geom_line(aes(y = mean, color = R0_factor, group = R0_factor),
            lwd = 1) +
  scale_color_manual(values = R0_colors) +
  scale_fill_manual(values = R0_colors)  +
  theme_cowplot() +
  labs(color = unname(TeX("$R_0$")), fill = unname(TeX("$R_0$"))) +
  ggtitle("A.") +
  scale_x_continuous(
    unname(TeX("Environmental noise strength $[\\sigma]$")),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    "Number of cases",
    breaks = c(2500, 5000, 7500, 10000),
    labels = c(2500, 5000, 7500, 10000),
    expand = expansion(add = c(10,400))
  ) +
  theme_minimal_grid(10) +
  guides(
    color = guide_legend(
      title = TeX("Basic reproduction number $[R_0]$:"),
      title.position = "top",
      label.position = "bottom",
      reverse = TRUE,
      direction = "horizontal",
      nrow = 1
    )
  ) +
  theme(
    plot.title = element_text(size = 10),
    legend.text = element_text(size = 7),
    legend.spacing = unit(0, "pt"),
    legend.title = element_text(size = 8),
    legend.position = "top",
    legend.justification = "right",
    legend.direction = "horizontal",  
    legend.spacing.y = unit(0, "pt"),
    legend.key.width = unit(14, "pt"),
    legend.key.height = unit(1, "pt"),
    plot.title.position = "panel",
    legend.margin = margin(t = -17, b = -5)
  )

## Subplot B) Heatmap for full model
#  = half-size, on top-right
#  * label clearly with subtitle
#  * center colorbar with below figure
Peak_cases_heat_temp <- generic_heat_function("max_cases", "all")

legend_temp = get_legend(Peak_cases_heat_temp)

Peak_cases_heat = Peak_cases_heat_temp +
  labs(
    title = "B. Full model"
  ) +
  theme(
    plot.title = element_text(size = 10),
    axis.title.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),  # Remove x-axis line
    axis.title.y = element_blank()
  ) +
  guides(fill = "none") #+
# theme(legend.margin=margin(t=-0.5,l=0.0,b=-0.275,r=0.25, unit='cm'))

## Subplot C) Heatmap for environmental noise sub-model
#  = half-size, on bottom-right
#  * label clearly with subtitle
Peak_cases_heat_enviro <- generic_heat_function("max_cases", "enviro") +
  labs(
    title = "C. Environmental noise submodel"
  ) +
  theme(
    plot.title = element_text(size = 10),
    # axis.title.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    # axis.title.y = element_blank()
  ) +
  guides(fill = "none") #+
# theme(legend.margin=margin(t=-0.5,l=0.0,b=-0.275,r=0.25, unit='cm'))

Figure4 = ggpubr::ggarrange(Peak_cases_mean,
                            egg::ggarrange(
                              Peak_cases_heat, 
                              Peak_cases_heat_enviro, 
                              ncol = 1,
                              heights = c(0.475,0.525)
                            ),
                            legend_temp,
                            ncol = 3,
                            widths = c(0.4, 0.525, 0.075)
)

ggsave("./figures/New_Figure4.png", Figure4, width = 9, height = 4, units = "in", dpi = 1200)

# New Figure 5: Intensity plots ----

## Subplot A) Line plot
#  = full-size, on left
Duration_mean <- line_plots_df %>% 
  arrange(R0) %>% 
  filter(name == "duration") %>% 
  ggplot(aes(x = sigma)) +
  geom_line(aes(y = mean, color = R0_factor, group = R0_factor),
            lwd = 1) +
  scale_color_manual(values = R0_colors) +
  scale_fill_manual(values = R0_colors)  +
  theme_cowplot() +
  labs(color = unname(TeX("$R_0$")), fill = unname(TeX("$R_0$"))) +
  ggtitle("A.") +
  scale_x_continuous(
    unname(TeX("Environmental noise strength $[\\sigma]$")),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    "Time [years]",
    breaks = c(2.5, 5.0, 7.5, 10.0),
    expand = expansion(add = c(0,0.25))
  ) +
  theme_minimal_grid(10) +
  guides(
    color = guide_legend(
      title = TeX("Basic reproduction number $[R_0]$:"),
      title.position = "top",
      label.position = "bottom",
      reverse = TRUE,
      direction = "horizontal",
      nrow = 1
    )
  ) +
  theme(
    plot.title = element_text(size = 10),
    legend.text = element_text(size = 7),
    legend.spacing = unit(0, "pt"),
    legend.title = element_text(size = 8),
    legend.position = "top",
    legend.justification = "right",
    legend.direction = "horizontal",  
    legend.spacing.y = unit(0, "pt"),
    legend.key.width = unit(14, "pt"),
    legend.key.height = unit(1, "pt"),
    plot.title.position = "panel",
    legend.margin = margin(t = -17, b = -5)
  )

## Subplot B) Heatmap for full model
#  = half-size, on top-right
#  * label clearly with subtitle
#  * center colorbar with below figure
Duration_heat_temp <- generic_heat_function("duration", "all")

legend_temp = get_legend(Duration_heat_temp)

Duration_heat = Duration_heat_temp +
  labs(
    title = "B. Full model"
  ) +
  theme(
    plot.title = element_text(size = 10),
    axis.title.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),  # Remove x-axis line
    axis.title.y = element_blank()
  ) +
  guides(fill = "none") #+
# theme(legend.margin=margin(t=-0.5,l=0.0,b=-0.275,r=0.25, unit='cm'))

## Subplot C) Heatmap for environmental noise sub-model
#  = half-size, on bottom-right
#  * label clearly with subtitle
Duration_heat_enviro <- generic_heat_function("duration", "enviro") +
  labs(
    title = "C. Environmental noise submodel"
  ) +
  theme(
    plot.title = element_text(size = 10),
    # axis.title.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    # axis.title.y = element_blank()
  ) +
  guides(fill = "none") #+
# theme(legend.margin=margin(t=-0.5,l=0.0,b=-0.275,r=0.25, unit='cm'))

Figure5 = ggpubr::ggarrange(Duration_mean,
                            egg::ggarrange(
                              Duration_heat, 
                              Duration_heat_enviro, 
                              ncol = 1,
                              heights = c(0.475,0.525)
                            ),
                            legend_temp,
                            ncol = 3,
                            widths = c(0.4, 0.525, 0.075)
)

ggsave("./figures/New_Figure5.png", Figure5, width = 9, height = 4, units = "in", dpi = 1200)

















