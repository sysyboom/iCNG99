# ======================================
# 🌿 Essential Gene Prediction Evaluation (Refined Natural Green Edition)
# ======================================

library(readxl)
library(dplyr)
library(caret)
library(ggplot2)
library(VennDiagram)

# ========== 1️⃣ Read data ==========
df <- read_excel("/data/new/models/paper/gene_essentiality_comparison.xlsx")
df$model_predicted <- factor(df$model_predicted, levels = c(0,1))
df$experimental_essential <- factor(df$experimental_essential, levels = c(0,1))

# ========== 2️⃣ Confusion matrix ==========
cm <- confusionMatrix(df$model_predicted, df$experimental_essential, positive = "1")

TP <- as.numeric(cm$table[2,2])
TN <- as.numeric(cm$table[1,1])
FP <- as.numeric(cm$table[2,1])
FN <- as.numeric(cm$table[1,2])

# Build data for heatmap
cm_df <- as.data.frame(cm$table)
colnames(cm_df) <- c("Predicted", "Actual", "Count")
cm_df$Percent <- round(cm_df$Count / sum(cm_df$Count) * 100, 1)

# assign correct labels based on positions
cm_df$Label <- c("TN", "FN", "FP", "TP")[c(1,3,2,4)]
cm_df$Display <- paste0(cm_df$Label, "\n", cm_df$Count, " (", cm_df$Percent, "%)")

# ========== 3️⃣ Confusion Matrix Plot (Elegant Green) ==========
p1 <- ggplot(cm_df, aes(x = Actual, y = Predicted, fill = Count)) +
  geom_tile(color = "white", linewidth = 1.1, radius = unit(4, "pt"),  alpha = 0.8 ) +  # rounded edges
  geom_text(aes(label = Display), size = 5.2, fontface = "bold", color = "black", lineheight = 0.9) +
  scale_fill_gradient(
    low = "#99E2B4",   # very light green
    high = "#036666",  # dark green
    name = "Count"
  ) +
  coord_equal() +
  theme_minimal(base_size = 15) +
  theme(
    axis.text = element_text(face = "bold", color = "black", size = 13),
    axis.title = element_text(face = "bold", color = "black", size = 14),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Confusion Matrix for Essential Gene Prediction (YPAD)",
    x = "Experimental (Actual)",
    y = "Model Prediction"
  )
print(p1)

# ========== 5️⃣ Performance Metrics ==========
TPR <- TP / (TP + FN)
FPR <- FP / (FP + TN)
Precision <- TP / (TP + FP)
F1 <- (2 * Precision * TPR) / (Precision + TPR)
Accuracy <- (TP + TN) / (TP + TN + FP + FN)
MCC <- (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

metrics <- data.frame(
  Metric = c("Accuracy", "Precision", "Recall", "F1", "MCC"),
  Value = round(c(Accuracy, Precision, TPR, F1, MCC), 3)
)

# 1️⃣ 给 Precision 单独打一个标签
  metrics <- metrics %>%
  mutate(Group = ifelse(Metric == "Precision", "Precision", "Other"))

# 2️⃣ 统一颜色，用透明度高亮 Precision
bar_fill <- "#4D908E"  # 单一蓝灰色
  
p2 <- ggplot(metrics, aes(x = Metric, y = Value)) +
  geom_bar(
    aes(alpha = Group),
    stat  = "identity",
    width = 0.7,
    color = NA,
    fill  = bar_fill
  ) +
  geom_text(aes(label = sprintf("%.2f", Value)),
            vjust = -0.5, size = 4.2,
            fontface = "bold", color = "black") +
  scale_alpha_manual(
    values = c("Other" = 0.5, "Precision" = 1.0),
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0.02))  # 顶部空白很小
  ) +
  labs(
    title = "Model Essentiality Prediction Performance (YPAD)",
    x = NULL,
    y = "Value"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid       = element_blank(),
    panel.background = element_blank(),
    plot.background  = element_blank(),
    axis.line   = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.text.x = element_text(face = "bold", size = 11, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(face = "bold", size = 11, color = "black",
                                margin = margin(r = 5)),
    plot.title  = element_text(face = "bold", hjust = 0.5, size = 15,
                               margin = margin(b = 8)),
    plot.margin = margin(15, 20, 10, 20),
    legend.position = "none"
  )

ggsave("performance_metrics_precision_alpha.png", p2,
       width = 10, height = 3.2, dpi = 600, bg = "white")

print(p2)
