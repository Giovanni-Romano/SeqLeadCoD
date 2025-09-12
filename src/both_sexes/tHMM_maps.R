# Packages ----
options(warning = 1)
suppressPackageStartupMessages(
  {
    library(tidyverse)
    library(reshape2) # melt()
    library(patchwork) # plot_layout() and "+" between plots
    library(abind)
    
    # pkgs for plotting maps
    library(sf)
    library(rnaturalearth)
    library(rnaturalearthdata)
    library(countrycode)
    
  }
)
source("src/utils.R")
cmndargs = c("Gnedin", "20019")#commandArgs(trailingOnly = TRUE)

# Load objects ----
TYPE = cmndargs[1]
SEED = cmndargs[2]
NAME = paste0(TYPE, SEED)
OUTFOLDER = "output/tHMM/"
## Data ----
ghe = readRDS("data/rds/dataGHE.RDS") %>% unite("ID", c(CountryN, Sex), remove = FALSE) %>% 
  unite("IDshort", c(Country, Sex), remove = FALSE)
ghe.lab = unique(ghe %>% select(c("Age", "CauseS"))) %>% 
  group_by(Age) %>% 
  mutate(CauseI = as.integer(factor(CauseS)))
m = ghe.lab %>% group_by(Age) %>% summarise(m = n()) %>% pull(m)
ghe.w = ghe %>% 
  left_join(ghe.lab, by = c("Age", "CauseS")) %>% 
  select(c("ID", "Age", "Year", "CauseI"))  %>% 
  arrange(ID, Age, Year)   %>% 
  pivot_wider(names_from = Age, values_from = CauseI)
ghe.array = abind(map(split(ghe.w, ghe.w$Year), \(x) x %>% select(-Year) %>% column_to_rownames("ID")), along = 3)
country_info = ghe %>% select(ID, Region) %>% unique()

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
niter = (RDS$input$ctr$ctr_mcmc$nchain + RDS$input$ctr$ctr_mcmc$nburnin)
nburnin = RDS$input$ctr$ctr_mcmc$nburnin

# Show/Save flags ----
SHOW = TRUE
SAVE = TRUE
if (SAVE){
  IMGFOLDER = paste0("img/tHMM/", NAME, "/")
  if (!dir.exists(IMGFOLDER)){
    dir.create(IMGFOLDER, recursive = TRUE)
  }
}


# Load world map
world <- bind_rows(cbind(ne_countries(scale = "medium", returnclass = "sf"), "Sex" = "F"),
                   cbind(ne_countries(scale = "medium", returnclass = "sf"), "Sex" = "M"))

df = ghe_PPE %>% 
  select(Year, Country, Sex, Cl.ss) %>% 
  unique() %>% 
  rename(Cluster = Cl.ss, iso_a3 = Country)

# Merge map with cluster assignments
# Function to create maps per year/sex
plot_clusters <- function(year_input) {
  df_year <- df %>% filter(Year == year_input)
  
  clusters = df_year %>% select(Cluster) %>% unique %>% pull %>% sort
  palette = pals::glasbey(length(clusters))
  
  map_data <- world %>%
    left_join(df_year, by = c("iso_a3_eh" = "iso_a3", "Sex")) %>%  # Match on ISO3 codes
    mutate(Sex = if_else(Sex == "F", "Female", "Male"))
  
  robin_crs <- "+proj=robin"
  
  ggplot(map_data) +
    geom_sf(aes(fill = factor(Cluster)), color = "gray60", size = 0.1) +
    facet_wrap(~Sex) +
    scale_fill_manual(values = palette, na.value = "white") +
    # theme_minimal() +
    coord_sf(crs = robin_crs) +
    labs(title = paste("Cluster Assignment -", year_input),
         fill = "Cluster") +
    theme(legend.position = "none")
}



unique_years <- unique(df$Year)

plots <- lapply(unique_years, function(y) plot_clusters(y))

# If you want to save each plot:
for (i in seq_along(unique_years)) {
  ggsave(paste0("img/tHMM/maps/cluster_map_", unique_years[i], ".pdf"), plots[[i]], width = 12, height = 4)
}
