

###### First run MixSIAR model comparison (MixSIAR_shiner_models.R), then select best model (here using model 3), then run this code with best model's jags #######


library(tidyverse)
library(MixSIAR)
library(ggplot2)
library(patchwork)


# Extract posterior draws for individual-level proportions
p_array <- jags.3$BUGSoutput$sims.list$p.ind  # [iterations x individuals x sources]

# Get dimensions
n_iter <- dim(p_array)[1]
n_ind  <- dim(p_array)[2]
n_src  <- dim(p_array)[3]

# Get source names
source_names <- source$source_names

posterior_df <- map_dfr(1:n_src, function(i_src) {
  draws <- p_array[, , i_src]
  dimnames(draws) <- list(NULL, paste0("id_", 1:n_ind))
  as_tibble(draws) %>%
    pivot_longer(everything(), names_to = "id", values_to = "prop") %>%
    mutate(source = source_names[i_src])
})

# Get fork length data for each individual
tl_vec <- mix3$data$tl
id_vec <- paste0("id_", 1:length(tl_vec))
tl_df  <- tibble(id = id_vec, tl = tl_vec)

# Join fork length to posterior draws
posterior_df <- left_join(posterior_df, tl_df, by = "id")

# Summarize posterior draws per individual/source
summary_df <- posterior_df %>%
  group_by(source, id, tl) %>%
  summarise(
    median = median(prop),
    lower  = quantile(prop, 0.025),
    upper  = quantile(prop, 0.975),
    .groups = "drop"
  )

# Ensure correct capitalization of source levels
summary_df <- summary_df %>%
  mutate(source = case_when(
    source == "creek" ~ "Creek",
    source == "lake"  ~ "Lake",
    TRUE ~ source
  ))

tl_plot <- ggplot(summary_df, aes(x = tl, y = median, color = source, fill = source)) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = TRUE,
    alpha = 0.45
  ) +
  scale_color_manual(values = c("Lake" = "#D1775E", "Creek" = "#3A74B4"), name = "Sources") +
  scale_fill_manual(values = c("Lake" = "#D1775E", "Creek" = "#3A74B4"), name = "Sources") +
  labs(x = "Total Length (mm)", y = "Diet Proportion") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "right",
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.key.size = unit(0.8, "cm"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_blank(),
    plot.margin = unit(c(0.2, 0.2, 1.2, 0.2), "cm"),
    axis.title.x = element_text(vjust = -1.5),
    panel.spacing = unit(1.5, "lines")
  ) +
  ylim(0, 1.0)

tl_plot
graphics.off()

### season + location plots
df.stats <- output_stats(jags.3, mix3, source, output_options)


df_tidy <- df.stats %>%
  as_tibble(rownames = "param") %>%
  filter(str_detect(param, "^p\\.")) %>%
  separate(param, into = c("p", "month", "location", "source"), sep = "\\.") %>%
  mutate(
    month    = str_to_title(month),
    location = str_to_title(location),
    source   = str_to_title(source),
    month    = dplyr::recode(as.character(month),
                             "May" = "Spring",
                             "August" = "Summer",
                             "October" = "Fall"),
    month    = factor(month, levels = c("Spring", "Summer", "Fall"))
  ) %>%
  mutate(
    location = dplyr::recode(location,
                             "Creek" = "Creek Captured",
                             "Lake"  = "Lake Captured"),
    location = factor(location, levels = c("Creek Captured", "Lake Captured"))
  )

# Create the plot
posteriors <- ggplot(df_tidy, aes(x = month, fill = source, color = source)) +
  geom_boxplot(
    aes(
      ymin = `2.5%`, lower = `25%`, middle = `50%`,
      upper = `75%`, ymax = `97.5%`,
      group = interaction(month, source)
    ),
    stat = "identity",
    position = position_dodge(width = 0.75),
    width = 0.6,
    alpha = 0.45,
    size = 0.5
  ) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dashed", color = "grey50", linewidth = 0.5) +
  facet_wrap(~location, ncol = 1) +
  scale_fill_manual(values = c("Creek" = "#3A74B4", "Lake" = "#D1775E"), name = "Sources") +
  scale_color_manual(values = c("Creek" = "#3A74B4", "Lake" = "#D1775E"), name = "Sources") +
  coord_cartesian(ylim = c(0, 0.85)) +
  theme_bw() +
  xlab("Season") +
  ylab("Diet Proportion") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.key.size = unit(0.8, "cm"),
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 14),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_blank(),
    plot.margin = unit(c(0.2, 0.2, 1.2, 0.2), "cm"),
    axis.title.x = element_text(vjust = -1.5),
    panel.spacing = unit(0.8, "lines")
  ) +
  ylim(0, 1)

print(posteriors)

posteriors_clean <- posteriors + theme(legend.position = "none")

combined_plot <- posteriors_clean + plot_spacer() + tl_plot +
  plot_layout(ncol = 3, widths = c(1.6, 0.03, 1.1), guides = "collect") +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 14))
print(combined_plot)






