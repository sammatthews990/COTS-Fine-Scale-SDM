library(dplyr)
library(ggplot2)
library(tidyr)

files <- c(
  Baseline = "DataProcessing/outputs/model_performance_metrics_clean.rds",
  EcoCentroid = "DataProcessing/outputs/model_performance_metrics_clean_ecoCentroid.rds",
  PseudoReps_N5 = "DataProcessing/outputs/model_performance_metrics_clean_pseudorepsN5.rds"
)

df_overall <- purrr::imap_dfr(files, function(f, name) {
  readRDS(f)$Overall %>% mutate(Strategy = factor(name, levels = c("Baseline", "EcoCentroid", "PseudoReps_N5")))
})

df_oof <- df_overall %>% filter(Split == "Validation (Out-of-Fold)")

cat("=== Summary of Validation Metrics (Out-of-Fold) ===\n")
print(df_oof %>% select(Strategy, roc_auc, accuracy, sens, spec))

# Plot comparative ROC AUC and Sensitivity
df_long <- df_oof %>%
  select(Strategy, ROC_AUC = roc_auc, Accuracy = accuracy, Sensitivity = sens, Specificity = spec) %>%
  pivot_longer(cols = c(ROC_AUC, Accuracy, Sensitivity, Specificity), names_to = "Metric", values_to = "Value")

p <- ggplot(df_long, aes(x = Strategy, y = Value, fill = Strategy)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = round(Value, 3)), vjust = -0.3, size = 3.5, fontface = "bold") +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_fill_manual(values = c("Baseline" = "#95a5a6", "EcoCentroid" = "#3498db", "PseudoReps_N5" = "#2ecc71")) +
  theme_bw(base_size = 12) +
  labs(
    title = "COTS SDM Model Performance Comparison across Spatial Point Strategies",
    subtitle = "Grouped 5-Fold Spatial Cross-Validation by CullSiteName",
    x = "Spatial Sampling Strategy",
    y = "Metric Value"
  ) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

ggsave("DataProcessing/outputs/comparative_strategy_metrics.png", p, width = 9, height = 6, dpi = 300)
cat("Plot saved to DataProcessing/outputs/comparative_strategy_metrics.png\n")
