suppressPackageStartupMessages(library(tidyverse))
source("src/tHMM_utils.R")

# HAMMING DISTANCE FUNCTION ----
# It's already defined in utils.R. Useful to have it here to share this code
# with Raffaella.
hamming_matrix <- function(X) {
  
  if (sum(is.na(X)) > 0) stop("Currently not implemented with NA")
  
  # Ensure data.frame of factors, column-wise
  XD <- as.data.frame(X, stringsAsFactors = FALSE)
  XD[] <- lapply(XD, function(col) {
    f <- factor(col)
    f
  })
  
  # One-hot encode all columns; no intercept
  Z <- model.matrix(~ . - 1, data = XD,
                    contrasts = lapply(XD, contrasts, contrasts = FALSE))
  
  # Hamming = 0.5 * L1 distance on the one-hot representation
  D <- as.matrix(dist(Z, method = "manhattan")) / 2
  
  # Ensure exact zeros on diagonal (numerical nicety)
  diag(D) <- 0
  as.dist(D)
}

# ggplot theme ----
theme_set(theme_bw() + theme(text = element_text(size = 12),
                             legend.key.size = unit(0.4, "cm"), 
                             strip.text = element_text(size = 10, margin = margin(2, 2, 2, 2))))

# Load estimates of centers
mu_est_F = readRDS("output/tHMM/5yrs/PT/female/muestdf_20019.RDS")
mu_est_M = readRDS("output/tHMM/5yrs/PT/male/muestdf_20146.RDS")

# Load datasets
ghe_F = readRDS("data/rds/GHEdf_female_5yrs.RDS")
ghe_M = readRDS("data/rds/GHEdf_male_5yrs.RDS")

## Find best order for cluster labels ----
# Create a single matrix with all estimates
allmu = rbind(mu_est_F %>% mutate(Sex = "F"),
              mu_est_M %>% mutate(Sex = "M")) %>% 
  mutate(Year = case_when(Year == 1 ~ 2000,
                          Year == 2 ~ 2005,
                          Year == 3 ~ 2010,
                          Year == 4 ~ 2015,
                          Year == 5 ~ 2020)) %>% 
  select(-CauseI) %>% mutate(CauseI = as.integer(CauseS)) %>% select(-CauseS) %>% 
  unite("ID", Sex, Year, Cluster) %>% 
  pivot_wider(names_from = Age, values_from = CauseI) %>% 
  column_to_rownames("ID") %>% as.matrix()
colnames(allmu) = levels(ghe_F$Age)

# Hamming distance computed on allmu
HDM_allmu = hamming_matrix(allmu)

# Find order of rows of allmu
set.seed(1)
ph_allmu = pheatmap::pheatmap(allmu, cluster_rows = T, cluster_cols = F, silent = T,
                              clustering_method = "ward.D2", clustering_distance_rows = HDM_allmu)

order_both = order(ph_allmu$tree_row$order)

# Split results in male and female
are_F = (substr(ph_allmu$tree_row$labels, 1, 1) == "F")

ph_allmu_F = list()
ph_allmu_F$order = order_both[are_F]
ph_allmu_F$labels = substr(ph_allmu$tree_row$labels[are_F], 3, 100)

ph_allmu_M = list()
ph_allmu_M$order = order_both[!are_F]
ph_allmu_M$labels = substr(ph_allmu$tree_row$labels[!are_F], 3, 100)


# Save objects
saveRDS(ph_allmu_F, "output/tHMM/5yrs/PT/female/ph_allmu_20019.RDS")
saveRDS(ph_allmu_M, "output/tHMM/5yrs/PT/male/ph_allmu_20146.RDS")
