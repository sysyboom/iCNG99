library(readxl)
library(tidyverse)
library(emmeans)
library(RColorBrewer)

df_raw <- read_excel("curve.xlsx")


df_long <- df_raw %>%
  pivot_longer(cols = -Strain,
               names_to = c("Time", "Rep"),
               names_sep = "_",
               values_to = "OD") %>%
  mutate(
    Time_n = as.numeric(gsub("h", "", Time)),
    Time_f = factor(Time, levels = unique(Time[order(as.numeric(gsub("h", "", Time)))]))
  )

print(head(df_long))

last_time <- max(df_long$Time_n, na.rm = TRUE)


df_last <- df_long %>% filter(Time_n == last_time)
model_last <- lm(OD ~ Strain, data = df_last)

posthoc_last <- emmeans(model_last, trt.vs.ctrl ~ Strain, ref = "WT")$contrasts %>%
  as.data.frame() %>%
  rename(p.adj = p.value) %>%
  mutate(significance = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*",
    p.adj < 0.1  ~ sprintf("(p.adj=%.3f)", p.adj),
    TRUE ~ "ns"
  )) %>%
  mutate(adjust_method = "Dunnett")  

print(posthoc_last)

legend_map <- posthoc_last %>%
  mutate(Strain_orig = sub(" - WT", "", contrast)) %>%
  mutate(new_name = paste0(Strain_orig, " ", significance)) %>%
  select(Strain_orig, new_name)


df_summary <- df_long %>%
  group_by(Strain, Time_n) %>%
  summarise(mean_OD = mean(OD), sd_OD = sd(OD), .groups = "drop") %>%
  left_join(legend_map, by = c("Strain" = "Strain_orig")) %>%
  mutate(display_name = ifelse(Strain == "WT", "WT (Control)", new_name))


df_summary$display_name <- factor(df_summary$display_name,
                                  levels = c("WT (Control)",
                                             sort(unique(df_summary$display_name[df_summary$display_name != "WT (Control)"]))))

fixed_colors <- c("black", "#F3722C", "#F8961E", "#F9C74F", "#90BE6D", "#43AA8B", "#577590")
strain_names <- levels(df_summary$display_name)
color_map <- setNames(fixed_colors[1:length(strain_names)], strain_names)

print(color_map)

p <- ggplot(df_summary, aes(x = Time_n, y = mean_OD, color = display_name, group = Strain)) +
  geom_errorbar(aes(ymin = mean_OD - sd_OD, ymax = mean_OD + sd_OD),
                width = 1.0, alpha = 0.6, linewidth = 0.5) +
  geom_line(aes(linewidth = (Strain == "WT"))) +
  geom_point(aes(size = (Strain == "WT"))) +
  scale_linewidth_manual(values = c("TRUE" = 1.2, "FALSE" = 0.7), guide = "none") +
  scale_size_manual(values = c("TRUE" = 2.2, "FALSE" = 1.5), guide = "none") +
  scale_color_manual(values = color_map) +
  scale_x_continuous(breaks = unique(df_summary$Time_n),
                     limits = c(min(df_summary$Time_n) - 1, last_time + 3)) +
  theme_classic() +
  labs(
    title = "Growth Curves under the Drug-Tolerant State",
    subtitle = paste0("Statistical significance at ", last_time,
                      "h determined by Dunnett’s post-hoc test (adjusted P-values) vs wild type"),
    x = "Time (hours)",
    y = "OD600",
    color = NULL
  ) +
  theme(
    legend.position = c(0.8, 0.22),  
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.text = element_text(size = 8, face = "bold"),
    legend.key.height = unit(0.35, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    axis.line = element_line(linewidth = 0.6),
    plot.title = element_text(face = "bold", size = 12)
  )

print(p)

