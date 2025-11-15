# Packages ----
options(warning = 1)
suppressPackageStartupMessages(
  {
    library(tidyverse)
    library(reshape2) # melt()
    library(mcclust) # comp.psm()
    library(HDInterval) # HDInterval()
    library(patchwork) # plot_layout() and "+" between plots
    library(latex2exp)
    library(abind)
    library(TraMineR) # For Optimal Matching functions
    library(ggalluvial)
    # pkgs for plotting maps
    library(sf)
    library(rnaturalearth)
    library(rnaturalearthdata)
    library(countrycode)
  }
)
source("src/tHMM_utils.R")

# Load objects ----
SEX = "male"
SEED = "20146"
OUTFOLDER = paste0("output/tHMM/5yrs/PT/", SEX, "/")

## Model output ----
mu_est = readRDS(paste0(OUTFOLDER, "muestdf_", SEED, ".RDS"))
PPE = readRDS(paste0(OUTFOLDER, "PPE_all.RDS"))[[as.character(SEED)]]
res = readRDS(paste0(OUTFOLDER, "res_Gnedin", SEED, ".RDS"))

## Data ----
if (SEX == "female"){
  ghe.init = readRDS("data/rds/GHEdf_female_5yrs.RDS")
  population = readxl::read_excel("data/raw/population_WPP2024_FEMALE.xlsx", 
                                       col_types = c(rep("text", 2), rep("numeric", 2))) %>% 
    filter(Year %in% seq(2000, 2020, by = 5))
  ghe = ghe.init %>% 
    left_join(population %>% select(-CountryN), by = c("Country", "Year")) %>% 
    rename(Population = Total)
} else if (SEX == "male"){
  ghe.init = readRDS("data/rds/GHEdf_male_5yrs.RDS")
  population = readxl::read_excel("data/raw/population_WPP2024_MALE.xlsx", 
                                       col_types = c(rep("text", 2), rep("numeric", 2))) %>% 
    filter(Year %in% seq(2000, 2020, by = 5))
  ghe = ghe.init %>% 
    left_join(population %>% select(-CountryN), by = c("Country", "Year")) %>% 
    rename(Population = Total)
} else {
  stop("Wrong sex value")
}


ghe.lab = unique(ghe %>% select(c("Age", "CauseS"))) %>% 
  group_by(Age) %>% 
  mutate(CauseI = as.integer(factor(CauseS)))
m = ghe.lab %>% group_by(Age) %>% summarise(m = n()) %>% pull(m)

ghe_PPE.init = left_join(ghe, 
                         PPE %>% as.data.frame() %>% rownames_to_column("CountryN") %>% 
                           pivot_longer(-CountryN, names_to = "Year", values_to = "Cl") %>% 
                           mutate(Year = as.numeric(Year)),
                         by = c("CountryN", "Year")
) %>% 
  # Shorten some country name that are way too long
  mutate(CountryN = case_when(
    CountryN == "Bolivia (Plurinational State of)" ~ "Bolivia",
    CountryN == "Democratic People's Republic of Korea" ~ "North Korea",
    CountryN == "Democratic Republic of the Congo" ~ "Dem. Rep. Congo",
    CountryN == "Iran (Islamic Republic of)" ~ "Iran",
    CountryN == "Lao People's Democratic Republic" ~ "Laos",
    CountryN == "Micronesia (Federated States of)" ~ "Micronesia",
    CountryN == "Netherlands (Kingdom of the)" ~ "Netherlands",
    CountryN == "Saint Vincent and the Grenadines" ~ "St. Vinc. and the Gren.",
    CountryN == "United Kingdom of Great Britain and Northern Ireland" ~ "United Kingdom",
    CountryN == "Bolivarian Republic of Venezuela" ~ "Venezuela",
    TRUE ~ CountryN
  ),
  Cl = as.character(Cl))


## Define dimensions ----
ages = levels(ghe$Age)
nages = length(ages)
years = seq(2000, 2020, by = 5)
nyears = length(years)
transitions = paste(years[-nyears], years[-1], sep = "-")
ntransitions = length(transitions)
countries = ghe %>% pull(CountryN) %>% unique %>% sort
ncountries = length(countries)
causes = readr::read_csv2("data/raw/causes_5yrs.csv") %>% select(CauseS) %>% unique() %>% unlist() %>% unname() %>% sort() 
ncauses = length(causes)


# Show/Save flags ----
SHOW = TRUE
SAVE = TRUE
if (SAVE){
  IMGFOLDER = paste0("img/tHMM/5yrs/PT/", SEX, "/", SEED)
  if (!dir.exists(IMGFOLDER)){
    dir.create(IMGFOLDER, recursive = TRUE)
  }
}

# Palette ----
palette_causes = c(pals::alphabet(26)[-(4:5)], pals::alphabet2(26), pals::glasbey(16)[-4])[1:ncauses]
names(palette_causes) =  causes
palette_regions = c(
  "Western Offshoots" = "#EF476F",
  "Latin America & Caribbean" = "#F78C6B",
  "Post-Sovietic" = "#6a994e",
  "Europe (Non P-S)" = "#a7c957",
  "Middle East & North Africa" = "#a68a64",
  "Sub-Saharan Africa" = "#d4a017",
  "South Asia" = "#8ecae6",
  "East Asia & Pacific" = "#118AB2"
)



# Centers ----
## Compute cluster sizes ----
cluster_sizes = apply(PPE, 2, function(x) {
  table(factor(x, levels = 1:max(x, na.rm = TRUE)))
}) %>% 
  # I name Cl.OG the original labels of the clusters, as I will name "Cl" the 
  # one I'll find with pheatmap
  melt(varnames = c("Cl.OG"), value.name = "Size") %>%
  rename(Year = L1) %>% 
  mutate(Year = as.numeric(Year),
         Cl.OG = as.character(Cl.OG))

## Compute cluster Population sizes ----
cluster_pop = ghe_PPE.init %>% 
  select(Year, Country, Population, Cl) %>% 
  distinct() %>% 
  group_by(Year, Cl) %>% 
  summarise(Cl_Pop = sum(Population)) %>% 
  rename(Cl.OG = Cl)


## Find best order for cluster labels ----
# Create a single matrix with all estimates
allmu = mu_est %>% 
  mutate(Year = case_when(Year == 1 ~ 2000,
                          Year == 2 ~ 2005,
                          Year == 3 ~ 2010,
                          Year == 4 ~ 2015,
                          Year == 5 ~ 2020)) %>% 
  select(-CauseS) %>% 
  unite("ID", Year, Cluster) %>% 
  pivot_wider(names_from = Age, values_from = CauseI) %>% 
  column_to_rownames("ID") %>% as.matrix()
colnames(allmu) = ages

# Load order of rows of allmu
ph_allmu = readRDS(paste0(OUTFOLDER, "ph_allmu_", SEED, ".RDS"))


## Create df with some cluster info ----
#   1) conversion from old Cl lab to new one
#   2) cluster sizes
#   3) y-coordinates for tiles plots   
max_nclusters = max(c(PPE))
cluster_info = 
  # dataframe with Year_Cl.OG and position in pheatmap output
  data.frame(ID = ph_allmu$labels, order = ph_allmu$order) %>% 
  separate(ID, into = c("Year", "Cl.OG"), sep = "_") %>% 
  mutate(Year = as.integer(Year)) %>% 
  left_join(cluster_sizes , by = c("Year", "Cl.OG")) %>% 
  # Create variable with relative size of cluster with respect to total population (n = 183 countries)
  mutate(Sizeprop = Size / ncountries) %>% 
  # Move the really small clusters and the empty ones to the bottom
  mutate(order = if_else(Size <= 2, order + 1000, order)) %>%
  # For each year, I arrange the df by "order" so that every year the cluster
  # with the smallest "order" value will have label 1
  group_by(Year) %>% arrange(order) %>%
  mutate(Cl = factor(row_number(), levels = 1:max_nclusters)) %>% 
  # Add cluster_pop df to define height of blocks
  left_join(cluster_pop, by = c("Year", "Cl.OG")) %>%
  # Create variables for y-coordinates of tiles plots
  mutate(bottomy = cumsum(lag((Cl_Pop), default = 0)),
         topy = cumsum(Cl_Pop)) %>%
  ungroup()

# saveRDS(cluster_info,
#         paste0(OUTFOLDER, "cluster_info_", SEED, ".RDS"))

# Add new Cl lab to ghe_PPE 
ghe_PPE = ghe_PPE.init %>% rename(Cl.OG = Cl) %>% 
  left_join(cluster_info %>% select(Year, Cl.OG, Cl, Size, Cl_Pop), 
            by = c("Year", "Cl.OG"))

# Create dataframe with centers and new cluster lab
mu_df = mu_est %>% 
  rename(Cl.OG = Cluster) %>% 
  mutate(Age = factor(ages[Age], levels = ages),
         Year = case_when(Year == 1 ~ 2000,
                          Year == 2 ~ 2005,
                          Year == 3 ~ 2010,
                          Year == 4 ~ 2015,
                          Year == 5 ~ 2020),
         Cl.OG = as.character(Cl.OG)) %>% 
  full_join(cluster_info %>% 
              select(Year, Cl.OG, Cl), 
            by = c("Year", "Cl.OG")) %>% 
  select(Year, Age, Cl.OG, Cl, CauseI, CauseS)

## Regional composition of each cluster ----
propregion_df = ghe_PPE %>% 
  select(Year, CountryN, RegionN, Cl, Population, Cl_Pop) %>% 
  distinct() %>% 
  group_by(Year, Cl, RegionN, Cl_Pop) %>% 
  summarise(Pop_regionclust = sum(Population)) %>% 
  mutate(Prop_regionclust = Pop_regionclust / Cl_Pop) %>% 
  ungroup() %>%
  select(Year:RegionN, Prop_regionclust)

## Strength of the centroid ----
# How strong is the centroid for each (Year, Cl, Age)?
# We define strength for each age (for each year) as the % of countries in the 
# cluster which has Cause equal to the mode
strength_mu = ghe_PPE %>% 
  select(Year, CountryN, Age, CauseS, Cl) %>% 
  left_join(mu_df %>% select(Year, Age, CauseS, Cl) %>% rename(mu = CauseS),
            by = c("Year", "Cl", "Age")) %>% 
  group_by(Year, Cl, Age) %>% 
  summarise(strength_mode = mean(CauseS == mu)) %>% 
  ungroup()

list_plot_centers = list()

## Plot centers ----
for (y in seq_along(years)){
  
  yyyy = years[y]
  
  # Create tmp df for part of the plot regarding regional composition
  propregion_df.tmp = propregion_df %>% 
    filter(Year == yyyy) %>%
    # Need to join propregion_df with cluster_info to get the top and 
    # bottom y-coordinate for each cluster
    left_join(cluster_info %>% select(Cl, Year, bottomy, topy),
              by = c("Year", "Cl")) %>%
    group_by(Year) %>%
    arrange(Cl, Prop_regionclust) %>%
    mutate(topy2 = cumsum(Prop_regionclust*(topy - bottomy)),
           bottomy2 = lag(topy2, default = 0)) %>%
    ungroup() %>% 
    select(Year, RegionN, bottomy2, topy2, Cl) %>%
    # Create "fake" variable to use facets. labels = "" to avoid printing anything
    mutate(Type = factor("Region", labels = ""))
  
  plt_centers = mu_df %>% 
    filter(Year == yyyy) %>%
    # Add y-coordinates from cluster_info
    left_join(cluster_info %>% select(Year, Cl, bottomy, topy),
              by = c("Year", "Cl")) %>% 
    # Add centroid strength for transparency
    left_join(strength_mu, by = c("Year", "Age", "Cl")) %>% 
    # Define numeric x-coordinate for geom_rect
    mutate(leftx = as.integer(Age) - 0.5,
           rightx = as.integer(Age) + 0.5,
           Type = "Center") %>% 
    ggplot() +
    # Tiles for centroids
    geom_rect(aes(xmin = leftx, xmax = rightx, ymin = -bottomy, ymax = -topy, 
                  # Color by cause and transparency by strength of the centroid
                  fill = CauseS, alpha = strength_mode), 
              col = "black") +
    # Name of the cause within each tile
    geom_text(aes(x = (rightx+leftx)/2, y = -(bottomy+topy)/2, label = CauseS),
              angle = 90, size = 2) +
    # Scale for transparency
    #   - "range" determines the min and max used transparency (min = 0.1 to avoid 
    #       completely transparent tiles)
    #   - "limits" determines the min and max possible values assumed by the variable
    #       which determines the transparency (strength_mode)
    scale_alpha_continuous(range = c(0.1, 1), limits = c(0, 1), guide = "none") +
    # Scale for causes of death
    scale_fill_manual("Cause", values = palette_causes) +
    # Need new scale for fill for regions
    ggnewscale::new_scale_fill() +
    # Tile for regions
    geom_rect(data = propregion_df.tmp,
              aes(xmin = -1, xmax = 0, ymin = -bottomy2, ymax = -topy2, 
                  fill = RegionN),
              # Slightly lighter color of tiles' borders wrt causes of death
              col = "grey30") +
    # Scale fill for regions
    scale_fill_manual("Region", values = palette_regions) +
    # Scale for x-axis to print actual labels and not numeric coordinates
    scale_x_continuous(breaks = c(-0.5, 1:nages), labels = c("Region", ages)) +
    # Use facets to split Centroids/Regions ( ~Type) and add small gaps
    # between clusters (Cl ~).
    #   - scales = "free" to avoid different x-axis between Region and Centroid
    #     facets and different y-axis between different clusters (each has its
    #     own range of bottomy-topy)
    facet_grid(Cl ~ Type, scales = "free", space = "free") +
    labs(title = paste0("Year ", yyyy)) + 
    theme(
      # Reduce a bit title sizes
      title = element_text(size = 9),
      # Fix legend position and key sizes
      legend.position = "right",
      legend.key.height = unit(0.025, "npc"), legend.key.width = unit(0.025, "npc"),
      # Vertical gaps sizes between clusters
      panel.spacing.y = unit(2, "pt"),
      # Remove strips of clusters (background and then text)
      strip.background.y = element_blank(),
      strip.text.y = element_blank(),
      # Remove useless axis text/ticks/titles
      axis.title = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    # Put expand = FALSE to avoid useless spacing beyond limits of axis scales
    coord_cartesian(expand = FALSE)
  
  if (SAVE){
    if (!dir.exists(paste0(IMGFOLDER, "/heatmap_centers_Population"))){dir.create(paste0(IMGFOLDER, "/heatmap_centers_Population"))}
    
    ggsave(plot = plt_centers,
           filename = paste0("heatmap_centers_Population/heatmap_Pop_", yyyy, ".pdf"),
           path = IMGFOLDER,
           height = 6, width = 9)
  }
  
  list_plot_centers[[y]] = plt_centers
}

# Create common legend for all years. To do so create not-printed plot with
# all years and clusters, to gett al possible causes
lgd = ggpubr::get_legend(mu_df %>%
                           left_join(cluster_info %>% select(Year, Cl, bottomy, topy),
                                     by = c("Year", "Cl")) %>% 
                           mutate(leftx = as.integer(Age) - 0.5,
                                  rightx = as.integer(Age) + 0.5,
                                  Type = "Center") %>% 
                           ggplot() +
                           geom_rect(aes(xmin = leftx, xmax = rightx, ymin = -bottomy, ymax = -topy, 
                                         fill = CauseS), col = "black") +
                           scale_fill_manual("Cause", values = palette_causes) + 
                           guides(fill = guide_legend(nrow = 5))+
                           facet_grid(Year ~ Type, scales = "free", space = "free") + 
                           ggnewscale::new_scale_fill() +
                           geom_rect(data = propregion_df.tmp,
                                     aes(xmin = -1, xmax = 0, ymin = -bottomy2, ymax = -topy2, fill = RegionN),
                                     col = "grey30") +
                           scale_fill_manual("Region", values = palette_regions) +
                           guides(fill = guide_legend(nrow = 4))+
                           scale_x_continuous(breaks = c(-0.5, 1:nages), labels = c("Region", ages)) + 
                           theme(legend.position = "top",
                                 legend.key.spacing.x = unit(1, "cm")))

# Put all year-specific plots together
plt_centers_together = ggpubr::ggarrange(plotlist = list_plot_centers, 
                                         common.legend = T, legend.grob = lgd,
                                         ncol = 3, nrow = 2)

ggsave(plot = plt_centers_together,
       filename = "estimated_centers_Population.pdf",
       path = IMGFOLDER,
       width = 15, height = 12)


# Flow coloured by Region ----
df_flow = ghe_PPE %>% 
  select(Country, RegionN, Cl, Year, Population) %>% distinct() %>% 
  mutate(Year = as.integer(factor(Year))) %>% 
  group_by(Year) %>% 
  mutate(freq = Population/sum(Population)) %>% 
  ungroup()

df_flow_summary <- df_flow %>% 
  mutate(Cl = case_when(Cl == 1 ~ "01",
                        Cl == 2 ~ "02",
                        Cl == 3 ~ "03",
                        Cl == 4 ~ "04",
                        Cl == 5 ~ "05",
                        Cl == 6 ~ "06",
                        Cl == 7 ~ "07",
                        Cl == 8 ~ "08",
                        Cl == 9 ~ "09",
                        TRUE ~ Cl)) %>% 
  mutate(Year = as.integer(factor(Year))) %>%
  unite("ID", Cl:RegionN) %>% 
  mutate(ID = factor(ID,
                     levels = as.vector(t(outer(sprintf("%02d", 1:20),
                                                levels(df_flow$RegionN),
                                                paste, sep = "_")))))

IDs = df_flow_summary$ID %>% unique() %>% sort()
IDs_regions = substr(IDs, start = 4, stop = 100)
palette_regions_cluster = palette_regions[IDs_regions]
names(palette_regions_cluster) = IDs

plt_flow_regions = df_flow %>%
  ggplot() +
  geom_flow(map = aes(x = Year, y = freq, stratum = Cl, 
                      alluvium = Country, fill = RegionN),
            alpha = 0.65) + 
  geom_bar(data = df_flow_summary,
           mapping = aes(x = Year, y = freq, fill = ID),
           color = "transparent",
           linewidth = 0.01,
           width = 0.335,
           stat = "identity", position = "fill") +
  geom_stratum(map = aes(x = Year, y = freq, stratum = Cl,
                         alluvium = Country), color = "black", linewidth = 0.55, alpha = 0) +
  scale_x_continuous(expand = c(0.025, 0.025), labels = c("2000", "2005", "2010", "2015", "2020")) +
  scale_y_continuous(expand = c(0.01, 0.015)) +
  scale_fill_manual(values = c(palette_regions, palette_regions_cluster)) + 
  theme_bw() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

if (SHOW) plt_flow_regions

if (SAVE){
  width = 12; height = 8; zoom = 0.75
  filename = "flow_regions_Population.pdf"
  ggsave(filename = filename,
         plot = plt_flow_regions,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}



# Network plot ----
library(igraph)
Regions = ghe %>% select(CountryN, RegionN) %>% unique %>% pull(RegionN)
names(Regions) = ghe %>% select(CountryN, RegionN) %>% unique %>% pull(CountryN)

source("src/utils.R")
palette_cluster = ggsci::pal_d3("category20")(20)[-8]
names(palette_cluster) = 1:length(palette_cluster)

if (SAVE) {if (!dir.exists(paste0(IMGFOLDER, "/networks_Population"))){dir.create(paste0(IMGFOLDER, "/networks_Population"))}}

for (t in 1:nyears){
  
  Z = ghe_PPE %>% filter(Year  == years[t]) %>% select(Country, Cl) %>% distinct() %>% pull(Cl) %>% as.integer()
  names(Z) = ghe_PPE %>% filter(Year  == years[t]) %>% select(Country, Cl) %>% distinct() %>%  pull(Country) %>% as.character()
  Z = Z[rownames(res$Y)]
  
  HSM = as.matrix((19-hamming_matrix(res$Y[ , , t]))/19)
  HSM_clusters = matrix(NA, nrow = max(Z), ncol = max(Z))
  for (j1 in 1:max(Z)){
    for (j2 in 1:max(Z)){
      tmp = HSM[Z == j1, Z == j2]
      HSM_clusters[j1, j2] = if_else(is.matrix(tmp), mean(tmp[upper.tri(tmp, diag = F)]), mean(tmp))
    }
  }
  
  HSM_clusters[HSM_clusters < 2/19] = 0
  HSM_clusters = HSM_clusters
  
  t_standard <- df_flow %>% filter(Year == t) %>% group_by(Cl, RegionN) %>% 
    summarise(t = sum(freq)) %>% 
    pivot_wider(id_cols = Cl, names_from = RegionN, values_from = t) %>% 
    column_to_rownames("Cl") %>% as.matrix()
  t_standard[is.na(t_standard)] = 0
  
  # define the colors in the pie-charts
  values <- list()
  for (h in 1:max(Z)){values[[h]] <- t_standard[h , sort(colnames(t_standard))]}
  pie_colors <- list()
  for (h in 1:max(Z)) {pie_colors[[h]] <- palette_regions[sort(names(palette_regions))]}
  
  mark_colors = adjustcolor(palette_cluster, alpha.f = 0.7)
  
  # transform the block probability matrix into an igraph object
  net_Y <- graph.adjacency(HSM_clusters, mode=c("undirected"), weighted=TRUE, diag=FALSE)
  
  # node sizes are proportional to cluster cardinality
  V(net_Y)$size <- rowSums(t_standard)*250
  
  # edge sizes are proportional to the estimated block probabilities
  # Note: for graphical purposes, the block probabilities below 0.1 are not displayed
  E(net_Y)$width <- E(net_Y)$weight^2*50
  
  # additional graphical settings
  V(net_Y)$label <- NA
  V(net_Y)$frame.color <- "black"
  E(net_Y)$color <- "grey"
  
  # node positions are obtained via force–directed placement
  seed = 2 #c(1, 1, 3, 1, 1)[t]
  set.seed(seed)
  l <- layout_with_fr(net_Y, 
                      minx = rep(-1.15, max(Z)), maxx = rep(1.15, max(Z)),
                      miny = rep(-0.65, max(Z)), maxy = rep(0.65, max(Z)))
  # l <- norm_coords(l, ymin=-0.75, ymax=0.75, xmin=-1.25, xmax=1.25)
  
  mark.g = list(); for(i in 1:max(Z)){mark.g[[i]] = i}
  
  pdf(file = paste0(IMGFOLDER, "/networks_Population/network_Pop_", years[t], ".pdf"),
      width = 9, height = 6)
  par(mar=c(0,0,0,0)+0.1)
  plot(net_Y,
       rescale=F, layout=l, edge.curved=0.2,vertex.shape="pie",
       vertex.pie=values, vertex.pie.color=pie_colors, 
       mark.groups = mark.g,
       mark.col = mark_colors, 
       mark.border = NA,
       mark.expand = rowSums(t_standard)*125 + 10,
       mark.shape = 0.5)
  dev.off()
}