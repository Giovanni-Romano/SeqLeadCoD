# Packages ----
options(warning = 1)
suppressPackageStartupMessages(
  {
    library(tidyverse)
    library(reshape2) # melt()
    library(abind)
    library(TraMineR)
  }
)

# Load objects ----
## Data ----
ghe = readRDS("data/rds/dataGHE.RDS") %>% unite("ID", c(CountryN, Sex), remove = FALSE) %>% 
  unite("IDshort", c(Country, Sex), remove = FALSE)
ghe.lab = unique(ghe %>% select(c("Age", "CauseS"))) %>% 
  group_by(Age) %>% 
  mutate(CauseI = as.integer(factor(CauseS)))
country_info = ghe %>% select(ID, Region) %>% unique()
## MCMC ----
PPE = readRDS("output/tHMM/PPE_Gnedin20019.RDS")
mu_est = readRDS("output/tHMM/muest_Gnedin20019.RDS")
## Define dimensions ----
ages = levels(ghe$Age)
nages = length(ages)
years = 2000:2021
nyears = length(years)
transitions = paste(years[-nyears], years[-1], sep = "-")
IDs = ghe %>% pull(IDshort) %>% unique %>% sort
nIDs = length(IDs)
ncauses = length(levels(ghe$CauseS))
ntransitions = length(transitions)
niter = 30*1e3
nburnin = 15*1e3


# Palette causes ----
palette = c(pals::alphabet(26)[-(4:5)], pals::alphabet2(26), pals::glasbey(32)[-c(4, 18)])[1:ncauses]
names(palette) =  sort(table(ghe$CauseS %>% as.character()), decreasing = T) %>% names


# Create a single matrix with all estimates
allmu = do.call(rbind, mu_est)

# Create a vector with all cluster sizes
allsizes = apply(PPE, 2, function(x) {
  table(factor(x, levels = 1:max(x, na.rm = TRUE)))
}) %>% melt(varnames = c("Cl"), value.name = "Size") %>%
  rename(Year = L1) %>%
  mutate(ID = paste0("Y", Year, "_cl", Cl)) %>%
  select(ID, Size)

# allmu has integer coding for causes (for pheatmap) --> transform into char
{
  allmu_char= matrix(NA_character_, nrow = nrow(allmu), ncol = ncol(allmu),
                     dimnames = dimnames(allmu))
  for (a in colnames(allmu)) {
    df_sub <- ghe.lab[ghe.lab$Age == a, ]
    map_vals <- setNames(df_sub$CauseS, df_sub$CauseI)
    allmu_char[, a] <- as.character(map_vals[as.character(allmu[, a])])
  }
  All.Seq.Tr<-seqdef(allmu_char)
  All.costs<- seqsubm(All.Seq.Tr, method="TRATE")
  Diss.om <- as.matrix(as.dist(seqdist(All.Seq.Tr,method="OM",indel=1, sm=All.costs,norm=FALSE)))
}

# Find best order
ph = pheatmap::pheatmap(allmu, cluster_rows = T, cluster_cols = F, silent = T,
                        clustering_method = "ward.D2", clustering_distance_rows = as.dist(Diss.om))
# ser = seriation::seriate(as.dist(Diss.om))[[1]]

# All possible clusters to complete data.frame so that all years have same number of rows
all_clusters = paste0("cl", unique(c(PPE)))

center_df = cbind(as.data.frame(allmu_char), 
                  order = order(ph$tree_row$order)) %>%
  rownames_to_column(var = "ID") %>%
  left_join(allsizes, by = "ID") %>%
  separate(ID, into = c("Year", "Cl"), sep = "_") %>%
  mutate(Year = as.integer(gsub("Y", "", Year))) %>%
  pivot_longer(-c(Year, Cl, order, Size), 
               names_to = "Age", values_to = "Cause") %>% 
  mutate(Age = factor(Age, levels = ages))

# Complete imputing "fake" clusters
center_df.complete = center_df %>% 
  mutate(order = as.integer(as.character(order))) %>% 
  complete(Year, Age, Cl = as.character(all_clusters)) %>%
  # Impute order also for "fake" clusters
  group_by(Year, Age) %>% 
  mutate(order = ifelse(is.na(order),
                        max(order, na.rm = TRUE) + dense_rank(Cl[is.na(order)]),
                        order)) %>%
  # Move the really small clusters and the empty ones to the bottom
  mutate(order = if_else(is.na(Size) | Size <= 2, order + 1000, order)) %>%
  ungroup() %>% 
  mutate(order = fct_rev(factor(order)))

saveRDS(center_df, "output/tHMM/center_df_Gnedin20019.RDS")
saveRDS(center_df.complete, "output/tHMM/center_df_complete_Gnedin20019.RDS")

plt_centers = center_df.complete %>% 
  ggplot() +
  geom_tile(aes(x = as.integer(Age), y = order, fill = Cause)) +
  scale_fill_manual(values = palette, na.value = "transparent",
                    breaks = sort(names(palette))) +
  guides(fill = guide_legend(ncol = 2, override.aes = list(color = "grey30"))) +
  scale_x_continuous(breaks = seq(1, nages, by = 3), 
                     labels = ages[seq(1, nages, by = 3)]) +
  facet_wrap(~ Year, ncol = 5, scales = "free_y") +
  theme(
    legend.key.height = unit(0.025, "npc"), legend.key.width = unit(0.015, "npc"),
    legend.position = "right",
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    text = element_text(size = 12))

ggsave(plt_centers, 
       filename = "img/tHMM/Gnedin20019/centers.pdf",
       width = 12, height = 8)


# df from old clusters number to new clusters number
df_clusters = center_df %>% 
  select(Year, Cl, order) %>% rename(Cl_old = Cl) %>%
  distinct() %>%
  group_by(Year) %>% arrange(order) %>%
  mutate(Cl_new = 1:n()) %>% 
  ungroup() %>% 
  arrange(Year, Cl_new) %>% 
  mutate(Cl_old = as.integer(gsub("cl", "", Cl_old)),
         Cl_new = as.integer(Cl_new))


# df with data and clusters info
ppe_df = PPE %>% melt() %>% rename(Year = year, Cl_old = value) %>% 
  left_join(df_clusters, by = c("Year", "Cl_old"))

ghe_ppe = ghe %>% 
  left_join(ppe_df, by = c("ID", "Year"))


saveRDS(ghe_ppe, "output/tHMM/ghe_ppe_Gnedin20019.RDS")
saveRDS(df_clusters, "output/tHMM/df_clusters_Gnedin20019.RDS")
