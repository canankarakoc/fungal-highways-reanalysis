# Reanalysis of HGT frequency along mycelial networks
# Updated analysis for manuscript
# Original data from Tom Berthold experiments
# Canan Karakoc 30 December 2025

# Load required packages
library(ggplot2)
library(plyr)  # load before dplyr so dplyr's summarize is used
library(dplyr)
library(emmeans)
library(multcomp)
library(multcompView)
library(car)
library(cowplot)

# run from the repository folder; the data files go in data/ (see data/README.md)
dir.create("figures", showWarnings = FALSE)

# Set theme for publication-quality figures
# Theme 
# panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
mytheme <- theme_bw() +
  theme(axis.ticks.length = unit(.25, "cm")) +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14, color = "black"), axis.title = element_text(size = 16)) +
  theme(panel.border = element_rect(
    fill = NA, colour = "black",
    size = 1
  )) +
  theme(strip.text.x = element_text(size = 14), strip.background = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(panel.border = element_rect(
    fill = NA, colour = "black",
    linewidth = 1
  )) +
  theme(
    axis.text.x.top = element_blank(), axis.title.x.top = element_blank(),
    axis.text.y.right = element_blank(), axis.title.y.right = element_blank()
  ) +
  theme(
    axis.title.x = element_text(margin = margin(10, 0, 0)),
    axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
    axis.text.x = element_text(margin = margin(10, 0, 0, 0)),
    axis.text.y = element_text(margin = margin(0, 10, 0, 0))
  )

# Read data
data <- read.csv("data/Testdatensatz_dreiExperimente_noZero.csv", 
                 sep = ";", dec = ",", header = TRUE, na.strings = c("NA", ""))


# Clean up vector names for display
data$Vector_display <- factor(data$Vector, 
                              levels = c("negativecontrol", "botrytis", "cunning", "pythium", "glassfibre"),
                              labels = c("Neg. Control", "B. cinerea", "C. elegans", "P. ultimum", "Glass fibre"))

# Set inoculum as ordered factor
data$Inoculum <- factor(data$Inoculum, levels = c("low", "medium", "high"))

# Log-transform cell counts for analysis
data$log_Cellcount <- log10(data$Cellcount + 1)

# Remove NA values
data_clean <- data %>% filter(!is.na(Cellcount))

cat("=== Data Summary ===\n")
print(summary(data_clean))

cat("\n=== Sample sizes per group ===\n")
print(table(data_clean$Vector_display, data_clean$Inoculum))

# Statistical Analysis
cat("\n=== Statistical Analysis ===\n")

# Test assumptions
cat("\n--- Normality test (Shapiro-Wilk on log-transformed data) ---\n")
shapiro_result <- shapiro.test(data_clean$log_Cellcount)
print(shapiro_result)

cat("\n--- Levene's Test for homogeneity of variances ---\n")
levene_vector <- leveneTest(log_Cellcount ~ Vector_display, data = data_clean, center = mean)
print(levene_vector)

# Linear model - Main effects only (interaction did not improve model)
cat("\n--- Linear Model (Vector + Inoculum) ---\n")
model_main <- lm(log_Cellcount ~ Vector + Inoculum, data = data_clean)
print(anova(model_main))

# Model with interaction for comparison
model_full <- lm(log_Cellcount ~ Vector * Inoculum, data = data_clean)

cat("\n--- AIC Comparison (main effects vs. interaction) ---\n")
cat("Main effects model AIC:", AIC(model_main), "\n")
cat("Full model (with interaction) AIC:", AIC(model_full), "\n")

# Type II ANOVA
cat("\n--- Type II ANOVA ---\n")
anova_result <- Anova(model_main, type = "II")
print(anova_result)

# Post-hoc comparisons for Vector using emmeans
cat("\n--- Post-hoc comparisons (Tukey) ---\n")
emm <- emmeans(model_main, ~ Vector)
pairs_result <- pairs(emm, adjust = "tukey")
print(pairs_result)

# Get compact letter display - METHOD 1: Using multcompView
# Extract p-values from pairwise comparisons
pairwise_pvals <- as.data.frame(pairs_result)

# Create a named vector of p-values for multcompLetters
# Format: "group1-group2" = p.value
pval_vector <- setNames(pairwise_pvals$p.value, 
                        paste(pairwise_pvals$contrast))

# Alternative: use Tukey HSD directly
tukey_result <- TukeyHSD(aov(log_Cellcount ~ Vector, data = data_clean))
tukey_pvals <- tukey_result$Vector[, "p adj"]

# Get letters using multcompLetters
letters_result <- multcompLetters(tukey_pvals, threshold = 0.05)
cat("\n--- Compact Letter Display ---\n")
print(letters_result)

# Create data frame with letters for plotting
sig_letters <- data.frame(
  Vector = names(letters_result$Letters),
  letter = as.character(letters_result$Letters)
)

# Map to display names
vector_map <- c("negativecontrol" = "Neg. Control", 
                "botrytis" = "B. cinerea", 
                "cunning" = "C. elegans", 
                "pythium" = "P. ultimum", 
                "glassfibre" = "Glass fibre")

sig_letters$Vector_display <- factor(vector_map[sig_letters$Vector],
                                     levels = c("Neg. Control", "B. cinerea", "C. elegans", 
                                                "P. ultimum", "Glass fibre"))

# Calculate positions for significance letters
letter_positions <- data_clean %>%
  group_by(Vector_display) %>%
  summarize(ypos = max(log_Cellcount, na.rm = TRUE) + 0.25, .groups = "drop") %>%
  left_join(sig_letters %>% dplyr::select(Vector_display, letter), by = "Vector_display")

print(letter_positions)

# Color palette (grayscale for publication)
colors <- c("low" = "#BEBEBE", "medium" = "#787878", "high" = "#2E2E2E")

# Create the figure
p <- ggplot(data_clean, aes(x = Vector_display, y = log_Cellcount)) +
  # Boxplot
  geom_boxplot(aes(fill = Inoculum), outlier.shape = NA, alpha = 0.8, 
               position = position_dodge(width = 0.75), width = 0.65,
               color = "black", linewidth = 0.4) +
  # Individual data points
  geom_point(aes(fill = Inoculum), 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
             shape = 21, size = 1.8, alpha = 0.5, color = "black", stroke = 0.3) +
  # Add significance letters
  #geom_text(data = letter_positions, 
  #          aes(x = Vector_display, y = ypos, label = letter),
  #          size = 5, fontface = "bold", inherit.aes = FALSE) +
  # Labels
  labs(x = "Transport vector",
       y = expression(log[10](Transconjugants + 1)),
       fill = "Inoculum") +
  # Color scale
  scale_fill_manual(values = colors,
                    labels = c("Low", "Medium", "High")) +
  mytheme+
  # Theme adjustments
  theme(axis.text.x = element_text(face = "italic"),
        legend.position = c(0.12, 0.85),
        legend.background = element_rect(fill = "white", color = "grey50", linewidth = 0.3),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  # Adjust y-axis
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  coord_cartesian(ylim = c(-0.1, 4.3))

# Display plot
print(p)

# Save the figure
ggsave("figures/Figure_HGT_frequency.pdf", p, width = 9, height = 5.5, dpi = 300)
ggsave("figures/Figure_HGT_frequency.png", p, width = 9, height = 5.5, dpi = 300)


# Print model summary for manuscript
cat("Vector effect: F =", round(anova_result$`F value`[1], 2), 
    ", df =", anova_result$Df[1], "/", anova_result$Df[3],
    ", p <", ifelse(anova_result$`Pr(>F)`[1] < 0.001, "0.001", round(anova_result$`Pr(>F)`[1], 3)), "\n")
cat("Inoculum effect: F =", round(anova_result$`F value`[2], 2), 
    ", df =", anova_result$Df[2], "/", anova_result$Df[3],
    ", p =", round(anova_result$`Pr(>F)`[2], 3), "\n")

# Mean transconjugants by vector (back-transformed)
means <- data_clean %>%
  group_by(Vector_display) %>%
  summarize(
    mean_log = mean(log_Cellcount),
    mean_count = 10^mean(log_Cellcount) - 1,
    .groups = "drop"
  )
print(means)


# ============================================================
# PATCH: Figure 5 – Simulated transconjugant numbers
# ============================================================

# Read simulation data
fig5 <- read.csv("data/fig5_data.csv", header = TRUE)

# Rename groups to match manuscript species labels and set display order
fig5$Vector <- factor(fig5$Fungi,
                      levels = c("B", "C", "P", "G"),
                      labels = c("B. cinerea", "C. elegans", "P. ultimum", "Glass fibre")
)

# Compute per-group means for the mean-line overlay
fig5_means <- ddply(fig5, .(Vector), summarise,
                    mean_val = mean(TransPerPop)
)

# Build Figure 5 in the style of Figure 2
fig5_plot <- ggplot(fig5, aes(x = Vector, y = TransPerPop)) +
  
  # Grey interquartile boxes (no fill colour since there is no Inoculum factor)
  geom_boxplot(
    aes(group = Vector),
    fill          = "grey80",
    colour        = "grey40",
    outlier.shape = NA,   # hide outlier dots; all points shown via geom_jitter
    # coef default = 1.5 → whiskers extend to 1.5*IQR, matching Fig 2 convention
    width         = 0.55,
    lwd           = 0.4
  ) +
  
  # Individual replicate points (jittered)
  geom_jitter(
    width  = 0.15,
    size   = 1.2,
    alpha  = 0.5,
    colour = "grey30"
  ) +
  
  # Black horizontal segment for the arithmetic mean (mirrors Fig 2 style)
  stat_summary(
    fun      = mean,
    geom     = "crossbar",
    width    = 0.55,
    fatten   = 1.5,
    colour   = "black",
    lwd      = 0.6
  ) +
  
  # Axis labels
  labs(
    x = "Transport vector",
    y = expression("Transconjugants per total cells")
  ) +
  
  # Cowplot clean theme to match Fig 2
  theme_cowplot(font_size = 11) +
  theme(
    axis.text.x  = element_text(face = "italic", size = 10),
    axis.title   = element_text(size = 11),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3)
  )

# Save (same dimensions used for Fig 2 in the manuscript)
save_plot("figures/Fig5_simulated_transconjugants.pdf",
          fig5_plot, base_width = 5.5, base_height = 4.5)
save_plot("figures/Fig5_simulated_transconjugants.png",
          fig5_plot, base_width = 5.5, base_height = 4.5, dpi = 300)

print(fig5_plot)

#------------------------------------------------------
# Exploratory analysis: Hyphal traits vs HGT frequency
#------------------------------------------------------

# Read transconjugant data
data <- read.csv("data/Testdatensatz_dreiExperimente_noZero.csv", 
                 sep = ";", dec = ",", header = TRUE, na.strings = c("NA", ""))

# Calculate means per vector - BOTH log and linear scale
vector_means <- data %>%
  filter(!is.na(Cellcount)) %>%
  group_by(Vector) %>%
  summarize(
    mean_log_TC = mean(log10(Cellcount + 1)),
    sd_log_TC = sd(log10(Cellcount + 1)),
    mean_TC = mean(Cellcount),
    sd_TC = sd(Cellcount),
    n = n(),
    .groups = "drop"
  )

# Hyphal traits from Table 1 of the manuscript
hyphal_traits <- read.csv("data/hyphal_traits.csv", header = TRUE)

# Merge
analysis_data <- merge(vector_means, hyphal_traits, by = "Vector")

# Display labels
analysis_data$Vector_label <- factor(analysis_data$Vector,
                                     levels = c("botrytis", "cunning", "pythium", "glassfibre"),
                                     labels = c("B. cinerea", "C. elegans", "P. ultimum", "Glass fibre"))


traits <- c("contact_angle", "hyphal_diameter", "water_film_width", 
            "vector_coverage", "waterfilm_coverage")
trait_labels <- c("Water contact angle (°)", "Hyphal diameter (µm)", 
                  "Water film width (µm)", "Vector coverage (%)", 
                  "Water film coverage (%)")

for (i in seq_along(traits)) {
  trait <- traits[i]
  x <- analysis_data[[trait]]
  y_log <- analysis_data$mean_log_TC
  y_lin <- analysis_data$mean_TC
  
  r_log <- cor(x, y_log, use = "complete.obs")
  r_lin <- cor(x, y_lin, use = "complete.obs")
  
  cat(sprintf("%-24s %+.3f     %+.3f        %.3f\n", 
              trait, r_log, r_lin, abs(r_log - r_lin)))
}

# Calculate common y-axis range
y_range <- c(
  min(analysis_data$mean_log_TC) - 0.3,
  max(analysis_data$mean_log_TC) + 0.3
)

# Function to create individual correlation plot
create_corr_plot <- function(data, trait, trait_label, y_var = "mean_log_TC", 
                             y_label = expression(Mean~log[10](TC + 1)),
                             y_limits = y_range) {
  
  # Filter NA
  plot_data <- data[!is.na(data[[trait]]), ]
  
  # Calculate correlation
  r_val <- cor(plot_data[[trait]], plot_data[[y_var]])
  
  # Calculate position for r label (middle of data range)
  x_mid <- mean(range(plot_data[[trait]]))
  y_mid <- max(y_limits)  # Use fixed y center for consistency
  
  # Colors (grayscale matching main figure)
  colors <- c("B. cinerea" = "#E69F00", "C. elegans" = "#56B4E9", 
              "P. ultimum" = "#009E73", "Glass fibre" = "#666666")
  
  p <- ggplot(plot_data, aes_string(x = trait, y = y_var)) +
    # Trend line first (behind points)
    geom_smooth(method = "lm", se = FALSE, color = "grey50", 
                linetype = "dashed", linewidth = 0.8) +
    # Error bars
    geom_errorbar(aes(ymin = plot_data[[y_var]] - sd_log_TC/sqrt(n),
                      ymax = plot_data[[y_var]] + sd_log_TC/sqrt(n)),
                  width = diff(range(plot_data[[trait]]))*0.03,
                  color = "grey30", linewidth = 0.4) +
    # Points
    geom_point(aes(fill = Vector_label), shape = 21, size = 4, 
               color = "black", stroke = 0.8) +
    # Correlation coefficient in middle
    annotate("label", x = x_mid, y = y_mid, 
             label = sprintf("italic(r) == %.2f", r_val),
             parse = TRUE, size = 3.5, fill = "white", 
             label.size = 0.3, alpha = 0.9) +
    # Labels
    labs(x = trait_label, y = y_label) +
    # Color scale
    scale_fill_manual(values = colors, name = NULL) +
    # Theme matching main figure
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      axis.line = element_line(color = "black", linewidth = 0.4),
      plot.margin = margin(5, 10, 5, 5)
    ) +
    # Consistent y-axis range across all panels
    coord_cartesian(ylim = y_limits)
  
  return(p)
}

# Create plots for each trait
p1 <- create_corr_plot(analysis_data, "contact_angle", "Water contact angle (°)", y_limits = y_range)
p2 <- create_corr_plot(analysis_data, "hyphal_diameter", "Hyphal diameter (µm)", y_limits = y_range)
p3 <- create_corr_plot(analysis_data, "water_film_width", "Water film width (µm)", y_limits = y_range)
p4 <- create_corr_plot(analysis_data, "vector_coverage", "Vector coverage (%)", y_limits = y_range)
p5 <- create_corr_plot(analysis_data, "waterfilm_coverage", "Water film coverage (%)", y_limits = y_range)
p6 <- create_corr_plot(analysis_data, "lacunarity", "Lacunarity (λ)", y_limits = y_range)

# Create a shared legend
legend_data <- analysis_data[!is.na(analysis_data$lacunarity), ]
colors <- c("B. cinerea" = "#E69F00", "C. elegans" = "#56B4E9", 
            "P. ultimum" = "#009E73", "Glass fibre" = "#666666")

p_legend <- ggplot(analysis_data, aes(x = contact_angle, y = mean_log_TC, fill = Vector_label)) +
  geom_point(shape = 21, size = 4, color = "black", stroke = 0.8) +
  scale_fill_manual(values = colors, name = NULL) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.text = element_text(face = "italic", size = 10),
        legend.key.size = unit(0.8, "cm"))

# Extract legend
legend <- get_legend(p_legend)

# Combine plots
plot_grid_main <- plot_grid(
  p1, p2, p3, p4, p5, p6,
  ncol = 3, nrow = 2,
  labels = c("a", "b", "c", "d", "e", "f"),
  label_size = 12,
  label_fontface = "bold"
)

# Add legend at bottom left
final_plot <- plot_grid(
  plot_grid_main,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

# Add title
title <- ggdraw() + 
  draw_label("Exploratory Analysis: Hyphal Traits vs. HGT Frequency (n = 4 vectors)", 
             fontface = "bold", size = 12, x = 0.5)

final_with_title <- plot_grid(final_plot, ncol = 1, rel_heights = c(0.05, 1))

# Save
ggsave("figures/Figure_hyphal_traits_correlations.pdf", final_with_title, 
       width = 10, height = 7, dpi = 300)
ggsave("figures/Figure_hyphal_traits_correlations.png", final_with_title, 
       width = 10, height = 7, dpi = 300)

cat("Figure saved!\n")

# Also create a version with LINEAR scale for comparison

create_corr_plot_linear <- function(data, trait, trait_label) {
  
  plot_data <- data[!is.na(data[[trait]]), ]
  r_val <- cor(plot_data[[trait]], plot_data$mean_TC)
  x_mid <- mean(range(plot_data[[trait]]))
  y_mid <- max(range(plot_data$mean_TC))
  
  colors <- c("B. cinerea" = "#E69F00", "C. elegans" = "#56B4E9", 
              "P. ultimum" = "#009E73", "Glass fibre" = "#666666")
  
  p <- ggplot(plot_data, aes(x = .data[[trait]], y = mean_TC)) +
    geom_smooth(method = "lm", se = FALSE, color = "grey50", 
                linetype = "dashed", linewidth = 0.8) +
    geom_errorbar(aes(ymin = mean_TC - sd_TC/sqrt(n),
                      ymax = mean_TC + sd_TC/sqrt(n)),
                  width = diff(range(plot_data[[trait]]))*0.03,
                  color = "grey30", linewidth = 0.4) +
    geom_point(aes(fill = Vector_label), shape = 21, size = 4, 
               color = "black", stroke = 0.8) +
    annotate("label", x = x_mid, y = y_mid, 
             label = sprintf("italic(r) == %.2f", r_val),
             parse = TRUE, size = 3.5, fill = "white", 
             label.size = 0.3, alpha = 0.9) +
    labs(x = trait_label, y = "Mean transconjugants") +
    scale_fill_manual(values = colors, name = NULL) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      axis.line = element_line(color = "black", linewidth = 0.4),
      plot.margin = margin(5, 10, 5, 5)
    )
  
  return(p)
}

# Create linear scale plots
p1_lin <- create_corr_plot_linear(analysis_data, "contact_angle", "Water contact angle (°)")
p2_lin <- create_corr_plot_linear(analysis_data, "hyphal_diameter", "Hyphal diameter (µm)")
p3_lin <- create_corr_plot_linear(analysis_data, "water_film_width", "Water film width (µm)")
p4_lin <- create_corr_plot_linear(analysis_data, "vector_coverage", "Vector coverage (%)")
p5_lin <- create_corr_plot_linear(analysis_data, "waterfilm_coverage", "Water film coverage (%)")
p6_lin <- create_corr_plot_linear(analysis_data, "lacunarity", "Lacunarity (λ)")

plot_grid_lin <- plot_grid(
  p1_lin, p2_lin, p3_lin, p4_lin, p5_lin, p6_lin,
  ncol = 3, nrow = 2,
  labels = c("a", "b", "c", "d", "e", "f"),
  label_size = 12,
  label_fontface = "bold"
)

final_plot_lin <- plot_grid(plot_grid_lin, legend, ncol = 1, rel_heights = c(1, 0.08))

title_lin <- ggdraw() + 
  draw_label("Linear Scale: Hyphal Traits vs. HGT Frequency (n = 4 vectors)", 
             fontface = "bold", size = 12, x = 0.5)

final_with_title_lin <- plot_grid(final_plot_lin, ncol = 1, rel_heights = c(0.05, 1))

ggsave("figures/Figure_hyphal_traits_LINEAR.pdf", final_with_title_lin, 
       width = 10, height = 7, dpi = 300)
ggsave("figures/Figure_hyphal_traits_LINEAR.png", final_with_title_lin, 
       width = 10, height = 7, dpi = 300)

cat("Linear scale figure saved!\n")





