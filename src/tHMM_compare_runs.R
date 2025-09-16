library(tidyverse)

sex = "female"
OUTFOLDER = paste0("output/tHMM/", sex)
IMGFOLDER = paste0("img/", sex, "/compare_runs")
runs = list.files(OUTFOLDER, pattern = "PPE", full.names = F) %>% str_sub(5, 9)
nruns = length(runs)
years = 2000:2021
nyears = length(years)

# Load PPE ----
PPE = abind::abind(lapply(list.files(OUTFOLDER, pattern = "PPE", full.names = T), readRDS), along = 3)
dimnames(PPE)[[3]] = runs

# NMI PPE across runs ----
NMI = array(NA, dim = c(nruns, nruns, nyears))
dimnames(NMI) = list(Run1 = runs,
                     Run2 = runs,
                     Year = years)

for (y in 1:22){
  for (i in 1:5){
    for (j in 1:i){
      NMI[i, j, y] = NMI[j, i, y] = aricode::NMI(PPE[ , y, i], PPE[ , y, j])
    }
  }
}

plt_NMI = NMI %>% reshape2::melt() %>% 
  mutate(Run1 = factor(Run1), Run2 = factor(Run2)) %>% 
  ggplot(aes(x = Run1, y = Run2)) + 
  geom_raster(aes(fill = value)) +
  geom_text(aes(label = round(value, 2)), size = 2.75) + 
  facet_wrap(~Year, nrow = 3) + 
  scale_fill_viridis_c("NMI", limits = c(0, 1)) + 
  labs(x = "Seed run #1", y = "Seed run #2",
       title = "NMI between partition point estimates across 5 runs") + 
  coord_fixed(ratio = 1) +
  theme(legend.position = "NULL", 
        axis.text.x = element_text(angle = 90, vjust = 0.5))

width = 12; height = 6; zoom = 1 
ggsave(filename = "NMI.pdf", path = IMGFOLDER,
       plot = plt_NMI,
       width = width*zoom, height = height*zoom)

# Number of clusters ----
ncl = apply(PPE, c(2, 3), max)
plt_ncl = ncl %>% as.data.frame() %>% 
  rownames_to_column("Year") %>% 
  pivot_longer(-Year, names_to = "Run", values_to = "NCl") %>%
  mutate(Run = factor(Run), Year = as.integer(Year)) %>% 
  ggplot(aes(x = Year, y = NCl, col = Run)) +
  geom_line() + geom_point() + 
  theme_bw()

width = 12; height = 6; zoom = 1 
ggsave(filename = "ncl.pdf", path = IMGFOLDER,
       plot = plt_ncl,
       width = width*zoom, height = height*zoom)


# Posterior similarity matrices ----
PSM = lapply(list.files(OUTFOLDER, pattern = "PSM", full.names = T), 
             function(x) readRDS(x) %>% abind::abind(along = 3)) %>% 
  abind::abind(along = 4)
dimnames(PSM) = list(Country1 = rownames(PSM), Country2 = colnames(PSM),
                     Year = 2000:2021, Run = runs)


plt_PSM = list()
for (y in seq_along(years)){
  yyy = years[y] %>% as.character()
  
  PSM.y = PSM[ , , yyy, ]
  # order = data.frame(Country = rownames(PSM), order = hclust(as.dist(1 - PSM.y[ , , 1]), method = "ward.D")$order %>% order())
  order = hclust(as.dist(1 - PSM.y[ , , 1]), method = "ward.D")$order
  plt_PSM[[y]] = PSM.y %>% reshape2::melt() %>%
    mutate(Country1 = factor(Country1, levels = levels(Country1)[order]),
           Country2 = factor(Country2, levels = levels(Country2)[order])) %>% 
    ggplot(aes(x = Country1, y = Country2, fill = value)) + 
    geom_raster() + 
    facet_grid(~Run) + 
    labs(title = paste0("Co-clustering probabilities - ", yyy)) + 
    scale_fill_viridis_c("", limits = c(0, 1)) + 
    coord_fixed(ratio = 1) +
    theme_bw() + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                       axis.title = element_blank())
}

width=10; height = 3; zoom = 1
ggsave(filename = "PSM.pdf", path = IMGFOLDER,
       plot = gridExtra::marrangeGrob(plt_PSM, nrow = 1, ncol = 1, top = NULL), 
       height = height*zoom, width = width*zoom)

