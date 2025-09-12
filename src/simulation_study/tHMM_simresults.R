suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

SAVE = TRUE

allmetrics = readRDS("output/tHMM/simstudy_simplified/allmetrics.rds")
NMI = allmetrics$NMI
dimnames(NMI) = list("Replicate" = 1:nrow(NMI), "Year" = 1:ncol(NMI))
prec = allmetrics$prec_est
dimnames(prec) = list("Replicate" = 1:nrow(prec), "Year" = 1:ncol(prec))


DODGE = 0.75
COLORS = c("#387cbc", "#ff7c04")


df = NMI %>% 
  reshape2::melt(value.name = "NMI") %>%
  full_join(
    reshape2::melt(prec, value.name = "Precision"),
    by = c("Replicate", "Year")
  ) %>% 
  pivot_longer(cols = c("NMI", "Precision"), names_to = "Metric", values_to = "value")

df_interval = df %>% 
  mutate(Year = factor(Year, levels = sort(unique(Year)), labels = as.character(1:8))) %>% 
  group_by(Year, Metric) %>% 
  summarize(valmean = mean(value),
            valmed = median(value),
            valmin = min(value),
            valmax = max(value),
            valq05 = quantile(value, 0.05),
            valq25 = quantile(value, 0.25),
            valq75 = quantile(value, 0.75),
            valq95 = quantile(value, 0.95)) %>% 
  as.data.frame()


df %>% 
  mutate(Year = as.factor(Year)) %>% 
  ggplot() +
  geom_boxplot(aes(x = Year, y = value, fill = Metric), 
               outlier.size = 0.5) + 
  scale_fill_manual(values = COLORS, labels = c("Partition", "Cluster sequence")) +
  scale_color_manual(values = COLORS, labels = c("Partition", "Cluster sequence")) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.025))) +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 14),
        strip.background = element_rect(colour = "white"),
        legend.position = "bottom", legend.title = element_blank()) +
  labs(x = "Year", y = "Metric value")

plt_interval = df_interval %>% 
  mutate(Year = as.numeric(Year)) %>%
  ggplot(aes(x = Year, y = valmed, fill = Metric, col = Metric)) +
  geom_errorbar(aes(ymin = valmin, ymax = valmax), linewidth = 0.5, width = 0, position = position_dodge(width = DODGE)) +
  geom_errorbar(aes(ymin = valq25, ymax = valq75), linewidth = 1.75, width = 0, position = position_dodge(width = DODGE)) +
  geom_point(size = 1.5, col = "black", position = position_dodge(width = DODGE)) +
  scale_color_manual(values = COLORS, labels = c("Partition", "Cluster sequence")) +
  scale_fill_manual(values = COLORS, labels = c("Partition", "Cluster sequence")) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.025))) +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 14),
        strip.background = element_rect(colour = "white"),
        legend.position = "bottom", legend.title = element_blank())

plt_interval

if (SAVE) {
  zoom = 0.75
  width = 9
  height = 5
  filename = "interval_s2s.pdf"
  ggsave(filename = filename, plot = plt_interval, path = "img/tHMM/simstudy_simplified", 
         width = zoom * width, height = zoom * height)
}



# Table with range, interquartile range and median for Precision and NMI
table_NMI = df_interval %>% 
  filter(Metric == "NMI") %>% 
  select(Year, valmin, valq25, valmed, valq75, valmax) %>% 
  rename("Year" = Year,
         "Min" = valmin,
         "Q25" = valq25,
         "Median" = valmed,
         "Q75" = valq75,
         "Max" = valmax) %>% 
  column_to_rownames("Year") %>% t()

table_Precision = df_interval %>%
  filter(Metric == "Precision") %>% 
  select(Year, valmin, valq25, valmed, valq75, valmax) %>% 
  rename("Year" = Year,
         "Min" = valmin,
         "Q25" = valq25,
         "Median" = valmed,
         "Q75" = valq75,
         "Max" = valmax) %>% 
  column_to_rownames("Year") %>% t()

# Convert both into LaTex tables
table_latex = knitr::kable(list(table_NMI, table_Precision),
                           v.align = "t",
                           format = "latex", digits = 2,
                           booktabs = TRUE,
                           row.names = TRUE)
