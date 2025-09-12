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
  }
)
source("src/utils.R")
source("src/tHMM_utils.R")
cmndargs = c("Gnedin", "20019")#commandArgs(trailingOnly = TRUE)

# Set some theme aspects for ggplot
theme_set(theme(text = element_text(size = 12),
                legend.key.size = unit(0.4, "cm"), 
                strip.text = element_text(size = 10, margin = margin(2, 2, 2, 2))))

# Load objects ----
TYPE = cmndargs[1]
SEED = cmndargs[2]
NAME = paste0(TYPE, SEED)
OUTFOLDER = "output/tHMM/"
## Model output ----
RDS = readRDS(paste0(OUTFOLDER, "res_", NAME, ".RDS"))
C = RDS$output$C
alpha = RDS$output$alpha[-1, ]
mu_est = readRDS(paste0(OUTFOLDER, "muest_", NAME, ".RDS"))
PSM = readRDS(paste0(OUTFOLDER, "PSM_", NAME, ".RDS"))
PPE = readRDS(paste0(OUTFOLDER, "PPE_", NAME, ".RDS"))
df_clusters = readRDS(paste0(OUTFOLDER, "df_clusters_", NAME, ".RDS"))
center_df = readRDS(paste0(OUTFOLDER, "center_df_", NAME, ".RDS")) %>% 
  mutate(Cl = as.integer(gsub("cl", "", Cl)),
         order = factor(order, levels = sort(unique(order), decreasing = T))) %>% 
  left_join(df_clusters %>% select(Year, Cl_old, Cl_new), 
            by = c("Year", "Cl" = "Cl_old"))

center_df.complete = readRDS(paste0(OUTFOLDER, "center_df_complete_", NAME, ".RDS")) %>% 
  mutate(Cl = as.integer(gsub("cl", "", Cl)),
         order = factor(order, levels = sort(unique(order), decreasing = T))) %>% 
  left_join(df_clusters %>% select(Year, Cl_old, Cl_new), 
            by = c("Year", "Cl" = "Cl_old"))

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
ghe.array = abind::abind(map(split(ghe.w, ghe.w$Year), \(x) x %>% select(-Year) %>% column_to_rownames("ID")), along = 3)

# UN Regions (DO NOT EXECUTE, USELESS)
{
  un_regions = read_delim("data/raw/Regions_UN.csv",
                          delim = ";", escape_double = FALSE, col_types = cols(`M49 Code` = col_skip(),
                                                                               `ISO-alpha2 Code` = col_skip()),
                          trim_ws = TRUE) %>%
    select(`ISO-alpha3 Code`, `Sub-region Name`, `Intermediate Region Name`) %>%
    mutate(Region = if_else(!is.na(`Intermediate Region Name`),
                            `Intermediate Region Name`, `Sub-region Name`)) %>%
    rename(Country = `ISO-alpha3 Code`) %>% select(Country, Region) %>%
    mutate(Region = case_when(Region %in% c("Melanesia", "Micronesia",
                                            "Polynesia", "Australia and New Zealand") ~ "Oceania",
                              
                              Region == "Southern Asia" ~ "South Asia",
                              Region == "Eastern Asia" ~ "East Asia",
                              Region == "Western Asia" ~ "West Asia",
                              Region == "South-eastern Asia" ~ "South East Asia",
                              Region == "Central Asia" ~ "Centr. Asia",
                              
                              Region == "Northern Africa" ~ "North Africa",
                              Region == "Eastern Africa" ~ "East Africa",
                              Region == "Middle Africa" ~ "Middle Africa",
                              Region == "Western Africa" ~ "West Africa",
                              Region == "Southern Africa" ~ "South Africa",
                              
                              Region %in% c("Caribbean", "Central America", "South America")  ~ "Latin America",
                              Region == "Northern America" ~ "North America",
                              
                              Region == "Eastern Europe" ~ "East Europe",
                              Region == "Northern Europe" ~ "North Europe",
                              Region == "Southern Europe" ~ "South Europe",
                              Region == "Western Europe" ~ "West Europe",
                              
                              TRUE ~ as.character(Region))) %>%
    # make Region a factor putting levels of same area together (Africa, America, Asia, Europe, Oceania)
    mutate(Region = factor(Region,
                           levels = c("North Africa", "East Africa", "Middle Africa", "West Africa", "South Africa",
                                      "Centr. America", "South America", "North America",
                                      "East Europe", "North Europe", "South Europe", "West Europe",
                                      "South Asia", "East Asia", "West Asia", "South East Asia", "Centr. Asia",
                                      "Oceania")))
}

# My regions
{
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  postsoviet = c("ARM", "AZE", "BLR", "EST", "GEO", "KAZ", "KGZ", "LVA", 
                 "LTU", "MDA", "RUS", "TJK", "TKM", "UKR", "UZB") 
  myregions = world %>% 
    as.data.frame() %>% 
    select(iso_a3_eh, region_wb, subregion) %>% 
    mutate(postsoviet = if_else(iso_a3_eh %in% postsoviet, "Yes", "No")) %>%
    mutate(myregions = case_when(postsoviet == "Yes" ~ "Post-Sovietic",
                                 iso_a3_eh == "TUR" ~ "Middle East & North Africa",
                                 region_wb == "Europe & Central Asia" ~ "Europe (Non P-S)",
                                 iso_a3_eh %in% c("USA", "CAN", "AUS", "NZL") ~ "Western Offshoots",
                                 TRUE ~ region_wb)) %>% 
    mutate(myregions = factor(myregions,
                              levels = c("Western Offshoots", "Latin America & Caribbean",
                                         "Sub-Saharan Africa", "Middle East & North Africa", 
                                         "Europe (Non P-S)", "Post-Sovietic",
                                         "South Asia", "East Asia & Pacific"))) %>%
    select(iso_a3_eh, myregions) %>%
    unique() 
}

alternative_regions = unique(myregions) %>% rename(Country = iso_a3_eh,
                                                   Region = myregions)

country_info = ghe %>% select(ID, Country, Region) %>% unique() %>% 
  rename(Region_WHO = Region) %>% 
  left_join(alternative_regions, by = "Country")

ghe_PPE = left_join(readRDS(paste0(OUTFOLDER, "ghe_ppe_", NAME, ".RDS")) %>% rename(Region_WHO = Region), 
                    country_info %>% select(ID, Region),
                    by = c("ID")) %>% 
  # Shorten same country names that are way too long
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
  )) %>%
  # Redefine ID accordingly
  rename(ID.long = ID) %>% unite("ID", CountryN, Sex, sep = "_", remove = F) %>%
  # Shorten regions
  mutate(Region.short = factor(case_when(
    Region == "Western Offshoots" ~ "WO",
    Region == "Latin America & Caribbean" ~ "LAC",
    Region == "Sub-Saharan Africa" ~ "SSA",
    Region == "Middle East & North Africa" ~ "MENA",
    Region == "Europe (Non P-S)" ~ "EU",
    Region == "Post-Sovietic" ~ "P-S",
    Region == "South Asia" ~ "SA",
    Region == "East Asia & Pacific" ~ "EAP",
    TRUE ~ Region
  ), levels = c("WO", "LAC", "SSA", "MENA", "EU", "P-S", "SA", "EAP")))


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
SAVE = FALSE
if (SAVE){
  IMGFOLDER = paste0("img/tHMM/", NAME, "/")
  if (!dir.exists(IMGFOLDER)){
    dir.create(IMGFOLDER, recursive = TRUE)
  }
}

# Palette causes ----
palette = c(pals::alphabet(26)[-(4:5)], pals::alphabet2(26), pals::glasbey(32)[-c(4, 18)])[1:ncauses]
names(palette) =  sort(table(ghe$CauseS %>% as.character()), decreasing = T) %>% names
# # 0) Read dataframe with groups
# df_groups = read.csv2("data/raw/causes.csv") %>%
#   select(CauseS, GioGroup) %>% rename(Cause = CauseS, Group = GioGroup)
# 
# # 1) Define one base color per GioGroup:
# base_hues <- c(
#   Birth = "#7FFFD4", # dark blue
#   Cancers = "blue", # orange
#   Cardio = "red", # red
#   Digestive = "#a68a64", # brown
#   DisturbBehav = "#9932cc",
#   Female = "pink",
#   InfectPara = "green", # light green
#   InfectResp = "orange", # yellow
#   Injury = "#999", # light yellow
#   Male = "lightblue",
#   Other = "yellow", # light green
#   SuddenEvents = "#BC8F8F" # light yellow
# )
# 
# # 2) For each group, generate n shades from light → base → dark:
# df_colors <- df_groups %>%
#   group_by(Group) %>%
#   mutate(n_in_group = n(),
#          fill_hex = {
#            n_col = n_in_group[1]
#            pal_fun <- scales::pal_seq_gradient(
#              scales::col_lighter(base_hues[cur_group_id()], -12*n_col/5),
#              scales::col_lighter(base_hues[cur_group_id()], +12*n_col/5),
#              "Lab"
#            )
#         pal_fun(seq(0, 1, length.out = n_col))
#          }) %>%
#   ungroup()
# 
# # 3) Build your named vector for scale_fill_manual:
# palette <- setNames(df_colors$fill_hex, df_colors$Cause)

# # 4) And finally in ggplot:
# ggplot(df_colors, aes(x=1, y=Cause, fill=Cause)) +
#   geom_tile() +
#   scale_fill_manual(values = palette) +
#   theme_minimal() +
#   facet_wrap(~ Group, ncol = 3, scales = "free")

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


# Traceplots ----
## Traceplots number of clusters ----
ncl_trace = apply(C, c(2, 3), function(x) length(unique(x)))
dimnames(ncl_trace) = list("year" = years, "iteration" = 1:niter)
plt_trcncl = ncl_trace %>% melt(value.name = "ncl") %>% 
  filter(iteration > 100) %>%
  ggplot(aes(x = iteration, y = ncl)) +
  geom_line() +
  geom_vline(xintercept = nburnin, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = c(101, seq(5000, niter, by = 5000))) +
  facet_wrap(~ year) +
  labs(x = "Iteration", y = "Number of clusters")

if (SHOW) plt_trcncl

if (SAVE){
  width = 12; height = 8; zoom = 1
  filename = "trcncl.pdf"
  ggsave(filename = filename,
         plot = plt_trcncl,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}

## Traceplot alpha ----
# Create a vector pasting years[x] with years[x+1] with a dash in between
dimnames(alpha) = list("transition" = transitions, 
                       "iteration" = 1:niter)
plt_trcalpha = alpha %>% melt(value.name = "alpha") %>% 
  ggplot(aes(x = iteration, y = alpha)) +
  geom_line() +
  geom_vline(xintercept = nburnin, linetype = "dashed", color = "red") +
  facet_wrap(~ transition, scales = "fixed") +
  labs(x = "Iteration", y = "Alpha")

if (SHOW) plt_trcalpha

if (SAVE){
  width = 12; height = 8; zoom = 1
  filename = "trcalpha.pdf"
  ggsave(filename = filename,
         plot = plt_trcalpha,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}

# Summary of alphas ----
summary_alpha = apply(alpha[ , -(1:nburnin)], 1, function(x) {
  HDI = hdi(x, credMass = 0.95) %>% unname()
  c("mean" = mean(x), "sd" = sd(x), HDI_l = HDI[1], HDI_u = HDI[2])
}) %>% t() %>% as.data.frame()

if (SHOW) View(summary_alpha %>% round(2))

if (SAVE){
  totxt <- knitr::kable(summary_alpha %>% round(2), format = "simple")  # Or "markdown"
  writeLines(totxt, 
             con = paste0(IMGFOLDER, "summary_alpha.txt"))
}

# Plot summary with mean as points and HDI (black) and sd (gold) as error bars
plt_summaryalpha = summary_alpha %>% 
  rownames_to_column("transition") %>% 
  mutate(transition = factor(transition, levels = transitions)) %>% 
  ggplot(aes(y = transition, x = mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = HDI_l, xmax = HDI_u), width = 0.2, lty = 1) +
  # geom_errorbar(aes(xmin = mean-sd, xmax = mean+sd), width = 0.2, col = "gold2", lty = 1) +
  coord_flip() +
  # labs(y = "Transition", x = "Alpha", title = "Posterior mean with 95% HDI (black) and sd (gold)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        text = element_text(size = 10))

if (SHOW) plt_summaryalpha
if (SAVE) {
  width = 12; height = 8; zoom = 0.35
  filename = "summary_alpha.pdf"
  ggsave(filename = filename,
         plot = plt_summaryalpha,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}


# Compute NMI between PPE in different years
NMI_PPE = matrix(NA, nrow = nyears, ncol = nyears)
dimnames(NMI_PPE) = list(Year1 = years, Year2 = years)
for (t in 1:nyears){
  for (s in t:nyears){
    NMI_PPE[t, s] = NMI_PPE[s, t] = aricode::NMI(PPE[ , t], PPE[ , s])
  }
}

NMI_PPE.consecutiveyears = sapply(2:nyears, function(j) NMI_PPE[j-1, j])

plt_summaryalpha = cbind(summary_alpha, nmi = NMI_PPE.consecutiveyears) %>% 
  rownames_to_column("transition") %>% 
  mutate(transition = factor(transition, levels = transitions)) %>% 
  ggplot(aes(y = transition)) +
  geom_point(aes(x = mean), size = 3) +
  geom_point(aes(x = nmi), size = 3, pch = 4, col = "red") +
  geom_errorbar(aes(xmin = HDI_l, xmax = HDI_u), width = 0.2, lty = 1) +
  # geom_errorbar(aes(xmin = mean-sd, xmax = mean+sd), width = 0.2, col = "gold2", lty = 1) +
  coord_flip() +
  # labs(y = "Transition", x = "Alpha", title = "Posterior mean with 95% HDI (black) and sd (gold)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        text = element_text(size = 10))

# Prob. CC and PPE ----
PPE.ncl = apply(PPE, 2, max)

plt_trcncl2 = plt_trcncl + 
  geom_hline(data = data.frame(year = years, PPE.ncl = PPE.ncl),
             mapping = aes(yintercept = PPE.ncl), col = "gold2", lty = 2)

if (SHOW) plt_trcncl2

if (SAVE){
  width = 12; height = 8; zoom = 1
  filename = "trcncl2.pdf"
  ggsave(filename = filename,
         plot = plt_trcncl2,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}

# Plot CENTERS ----
plt_centers = center_df.complete %>% 
  # filter(Year == 2020) %>% 
  ggplot() +
  geom_tile(aes(x = as.integer(Age), y = fct_rev(order), fill = Cause)) +
  # geom_text(aes(x = as.integer(Age), y = fct_rev(order), label = Cause),
  #           size = 2, color = "black", check_overlap = TRUE) +
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

if(SHOW) plt_centers

if (SAVE){
  width = 12; height = 8; zoom = 1
  ggsave(plt_centers, 
         filename = "img/tHMM/Gnedin20019/centers.pdf",
         width = width, height = height)
  ggsave(plt_centers + 
           theme(legend.position = "none"), 
         filename = "img/tHMM/Gnedin20019/centers_nolgd.pdf",
         width = width, height = height)
  
}

for (t in 2000:2021){
  plt_centers_1Y = center_df.complete %>% 
    filter(Year == t) %>% 
    ggplot() +
    geom_tile(aes(x = as.integer(Age), y = fct_rev(order), fill = Cause)) +
    # geom_text(aes(x = as.integer(Age), y = fct_rev(order), label = Cause),
    #           size = 2, color = "black", check_overlap = TRUE) +
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
      text = element_text(size = 20),
      strip.text = element_text(size = 20))
  
  
  if (SAVE){
    width = 12; height = 8; zoom = 0.75

    name = paste0("img/tHMM/Gnedin20019/centers/centers_", t, ".pdf")
    ggsave(plt_centers_1Y + 
             theme(legend.position = "none"), 
           filename = name,
           width = width*zoom, height = height*zoom)
    
  }
}


# Try a bigger plot with legend at the bottom
plt_centers2 = center_df.complete %>% 
  ggplot() +
  geom_tile(aes(x = as.integer(Age), y = fct_rev(order), fill = Cause)) +
  scale_fill_manual(values = palette, na.value = "transparent",
                    breaks = sort(names(palette))) +
  guides(fill = guide_legend(nrow = 6, override.aes = list(color = "grey30"))) +
  scale_x_continuous(breaks = seq(1, nages, by = 3), 
                     labels = ages[seq(1, nages, by = 3)]) +
  facet_wrap(~ Year, ncol = 5, scales = "free_y") +
  theme(
    legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"),
    legend.position = "bottom", legend.direction = "horizontal",
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    text = element_text(size = 12))

if (SHOW) plt_centers2

if (SAVE){
  width = 9; height = 12; zoom = 1
  ggsave(plt_centers2, 
         filename = "img/tHMM/Gnedin20019/centers2.pdf",
         width = width, height = height)
}

# plot just to visualize colors-causes association
palette %>% as.data.frame() %>% 
  rownames_to_column("Cause") %>% 
  ggplot(aes(x = 1, y = Cause, fill = Cause)) +
  geom_tile() +
  geom_text(aes(label = Cause), size = 3, color = "yellow2") +
  scale_fill_manual(values = palette) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        text = element_text(size = 12))


# Plot data by PPE ----
plt_dataclust = list()
for (yyyy in years){
  cat(yyyy, "\t")
  plt = ghe_PPE %>% 
    filter(Year == yyyy) %>%
    arrange(Cl_new, Region, Sex) %>%
    mutate(ID = factor(ID, levels = unique(ID))) %>% 
    ggplot(aes(x = Age, y = ID, fill = CauseS)) +
    geom_raster() + 
    facet_grid(rows = vars(Cl_new), scales = "free_y", space = "free_y") + 
    geom_text(aes(label = CauseS), size = 2) +
    scale_fill_manual(values = palette) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  data_aux = unique(ghe_PPE %>% filter(Year == yyyy) %>%
                      select(ID, Cl_new, Sex, Region)) %>%
    arrange(Cl_new, Region, Sex) %>%
    mutate(ID = factor(ID, levels = unique(ID))) %>% 
    pivot_longer(cols = -c(ID, Cl_new), , names_to = "x", values_to = "fill")
  
  colors_aux = c("M"= "lightblue", "F" = "pink", palette_regions)
  
  plt_aux = data_aux %>% 
    ggplot() + 
    geom_tile(aes(x = x, y = ID, fill = fill)) +
    facet_grid(rows = vars(Cl_new), scales = "free_y", space = "free_y") + 
    scale_fill_manual("Region", values = colors_aux) +
    # geom_text(aes(x = x, y = ID, label = fill), size = 2) +
    scale_color_manual(values = c("black", "white"), breaks = c("black", "white")) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())
  
  plt_out = plt_aux + plt + plot_layout(widths = c(1, 20))
  plt_dataclust[[yyyy+1-2000]] = plt_out
  
  if (SAVE){
    width = 16; height = 40; zoom = 1
    if (!dir.exists(paste0(IMGFOLDER, "heatmap/"))){
      dir.create(paste0(IMGFOLDER, "heatmap/"), recursive = TRUE)
    }
    
    ggsave(plot = plt_out,
           filename = paste0("heatmap/heatmap_", yyyy, ".pdf"),
           path = IMGFOLDER,
           height = height*zoom, width = width*zoom)
  }
}

# Entropy of region and sex within cluster ----
# entropy_df = ghe_PPE %>% 
#   group_by(Year, Cl_new) %>%
#   mutate(H_region = DescTools::Entropy(table(Region)) / log(6, base = 2),
#          H_sex = DescTools::Entropy(table(Sex)) / log(2, base = 2))
# 
# plt_entropy = entropy_df %>% 
#   ggplot(aes(x = H_sex, y = H_region, color = factor(Cl_new))) +
#   # geom_point() +
#   geom_text(aes(label = Cl_new)) +
#   facet_wrap(~ Year) + 
#   labs(x = "Sex", y = "Region")
# 
# if (SHOW) plt_entropy
# 
# if (SAVE) {
#   width = 12; height = 8; zoom = 1
#   filename = "entropy.pdf"
#   ggsave(filename = filename,
#          plot = plt_entropy,
#          path = IMGFOLDER,
#          width = width*zoom, height = height*zoom)
# }


# Test independence PPE vs Region and Sex ----
pvalue_region = lapply(years, function(t) {tmp = ghe_PPE %>% filter(Year == t); chisq.test(tmp$Cl_new, tmp$Region, simulate.p.value = T, B = 10^4)})
pvalue_sex = lapply(years, function(t) {tmp = ghe_PPE %>% filter(Year == t); chisq.test(tmp$Cl_new, tmp$Sex, simulate.p.value = T, B = 10^4)})

# Plot composition of clusters wrt sex and region----
type = "stack" # "fill" for relative proportions, "stack" for counts
## Sex ----
prop_df = ghe_PPE %>%
  group_by(Year, Cl_new) %>%
  summarise(p_f = mean(Sex == "F"), .groups = "drop") %>% 
  # arrange(Year, p_f) %>%
  # mutate(Cl_new_ord = 1:n()) %>%
  # mutate(Cl_new_lab = paste(Year, Cl_new, sep = " - ")) %>%
  # select(Year, Cl_new_lab, Cl_new_ord)
  select(Year, Cl_new)

plt_propSex = ghe_PPE %>%
  select(Year, Cl_new, Country, Sex) %>% 
  unique() %>% 
  # mutate(Cl_new_lab = paste(Year, Cl_new, sep = " - ")) %>%
  left_join(prop_df %>% select(Year, Cl_new), by = c("Year", "Cl_new")) %>%
  # left_join(prop_df %>% select(Year, Cl_new_lab, Cl_new_ord), by = c("Year", "Cl_new_lab")) %>%
  ggplot(aes(x = Cl_new, fill = Sex)) +
  geom_bar(position = type, col = "grey30", linewidth = 0.25) +
  facet_wrap(~ Year, scales = "free_x", ncol = 6) +
  scale_fill_manual(values = c("F" = "pink", "M" = "lightblue"), 
                    labels = c("Female", "Male")) +
  guides(fill = guide_legend(nrow = 2)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(5/6, 1/8),
        legend.direction = "horizontal")

if (SHOW) plt_propSex

if (SAVE){
  width = 12; height = 6; zoom = 1
  ggsave(filename = paste0("clusters_propSex_", type, ".pdf"),
         plot = plt_propSex,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}

## Region ----
prop_df <- ghe_PPE %>% 
  group_by(Year, Cl_new) %>% 
  summarise(
    p_SSA = mean(Region == "Sub-Saharan Africa"),
    p_MENA = mean(Region == "Middle East & North Africa"),
    p_EAP = mean(Region == "East Asia & Pacific"),
    p_LAC = mean(Region == "Latin America & Caribbean"),
    p_EU = mean(Region == "Europe (Non P-S)"),
    p_PS = mean(Region == "Post-Sovietic"),
    p_NA = mean(Region == "North America"),
    p_SA = mean(Region == "South Asia"),
  ) %>% 
  select(Year, Cl_new)


plt_propRegion = ghe_PPE %>%
  left_join(prop_df, by = c("Year", "Cl_new")) %>%
  select(Year, ID, Region, Cl_new) %>% 
  unique() %>% 
  ggplot(aes(x = Cl_new, fill = Region)) +
  geom_bar(position = type, col = "grey30", linewidth = 0.25) +
  scale_fill_manual(values = palette_regions) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
  facet_wrap(~ Year, scales = "free_x", ncol = 6) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(5/6, 1/8),
        legend.direction = "horizontal")

if (SHOW) plt_propRegion

if (SAVE){
  width = 12; height = 6; zoom = 1
  ggsave(filename = paste0("clusters_propRegion_", type, ".pdf"),
         plot = plt_propRegion,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}

## Are male populations more homogenous? ----
# Number of clusters by sex
plt_nclsex = ghe_PPE %>%
  select(-c(Age, CauseC:CauseT)) %>%
  distinct() %>% 
  group_by(Year, Sex, Cl_new) %>% 
  summarise(n = n()) %>% 
  group_by(Year, Sex) %>% 
  summarise(n_cl = n()) %>% 
  ggplot() +
  geom_line(aes(x = Year, y = n_cl, colour = Sex)) +
  scale_y_continuous("Number of clusters", breaks = c(10, 15, 20, 25))


if (SHOW) plt_nclsex

if (SAVE){
  width = 12; height = 6; zoom = 0.75
  ggsave(filename = "nclsex.pdf",
         plot = plt_nclsex,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}

# helper functions
gini_index <- function(p) {
  # p: proportions (must sum to 1)
  # Gini = 1 - sum(p^2) is the "Gini-Simpson index" (diversity form)
  1 - sum(p^2)
}

shannon_entropy <- function(p) {
  # Shannon entropy in bits
  -sum(ifelse(p > 0, p * log2(p), 0))
}

# aggregate counts
plot_data = ghe_PPE %>%
  select(-c(Age, CauseC:CauseT)) %>% distinct %>% 
  count(Year, Sex, Region.short, Region, Cl_new) 

ncl_byyear = plot_data %>%
  group_by(Year) %>% 
  summarise(ncl_year = length(unique(Cl_new)))
  

stats_by_group = plot_data %>% 
  left_join(ncl_byyear, by = "Year") %>% 
  group_by(Year, Sex, Region.short) %>%
  mutate(p = n / sum(n)) %>%
  summarise(
    Region = unique(Region),
    ncl_year = unique(ncl_year),
    ncl = length(unique(Cl_new)),
    ncl_norm = ncl/ncl_year,
    H       = shannon_entropy(p),                 # entropy (bits)
    H_norm  = if (ncl_year > 1) H / log2(ncl_year) else 0,  # normalized to [0,1]
    Gini = gini_index(p),
    Gini_norm = if (ncl_year > 1) Gini * ncl_year / (ncl_year) else 0,
    .groups = "keep"
  )

# stacked bar plot: each bar = Region, filled by clusters
plt_clustregionsex.bar = ggplot(plot_data, 
       aes(x = Sex, y = n, 
           fill = tidytext::reorder_within(factor(Cl_new), n, interaction(Year, Sex, Region.short)))) +
  geom_bar(stat = "identity", position = "stack", col = "grey30", show.legend = F) +
  facet_grid(Region.short ~ Year, scale = "free_y") +
  labs(x = "Region", y = "Clusters' split",
       fill = "Cluster (Cl_new)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

if (SHOW){plt_clustregionsex.bar}

if (SAVE){
  width = 9; height = 6; zoom = 2
  ggsave(filename = "clustregionsex_bar.pdf",
         plot = plt_clustregionsex.bar,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}


plt_clustregionsex.indices = stats_by_group %>%
  pivot_longer(c("ncl", "Gini_norm", "H_norm"), names_to = "stat") %>%
  mutate(stat = factor(stat, levels = c("H_norm", "Gini_norm", "ncl"),
                       labels = c("Entropy (norm.)", "Gini (norm.)", "Nmbr. of clust."))) %>%
  ggplot() +
  geom_line(aes(x = Year, y = value, col = Sex)) +
  ggh4x::facet_nested(stat ~ Region.short, scales = "free_y") +
  ggh4x::facetted_pos_scales(y = list(scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0, 1)),
                                      scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0, 1)),
                                      scale_y_continuous(breaks = seq(5, 22, by = 5), limits = c(0, 15)))) +
  scale_x_continuous(minor_breaks = 2000:2022) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.major.x = element_line(linewidth = 1.25),
        axis.title.y = element_blank())


if (SHOW){plt_clustregionsex.indices}

if (SAVE){
  width = 12; height = 6; zoom = 1
  ggsave(filename = "clustregionsex_indices.pdf",
         plot = plt_clustregionsex.indices,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}

# Waffle Regions ----
{
  round_preserve_sum <- function(x, total = 100) {
    x_scaled <- x / sum(x) * total
    floored <- floor(x_scaled)
    remainder <- total - sum(floored)
    if (remainder > 0) {
      indices <- order(x_scaled - floored, decreasing = TRUE)[1:remainder]
      floored[indices] <- floored[indices] + 1
    }
    return(floored)
  }
  
  waffle_data <- ghe_PPE %>%
    select(Year, Country, ID, Region, Cl_new) %>% unique() %>% # Remove duplicates given by 19 ages
    left_join(prop_df, by = c("Year", "Cl_new")) %>%
    select(Year, Cl_new, Region) %>%
    count(Year, Cl_new, Region) %>% 
    group_by(Year, Cl_new) %>%
    mutate(pct = round_preserve_sum(n, total = 100)) %>%
    ungroup()
  
  
  # Step 1: Detect truly missing clusters
  observed_clusters <- waffle_data %>%
    distinct(Year, Cl_new)
  
  # Step 2: Build complete cluster set
  full_clusters <- expand.grid(
    Year = unique(ghe_PPE$Year),
    Cl_new = 1:22
  )
  
  # Step 3: Identify missing ones
  missing_clusters <- anti_join(full_clusters, observed_clusters, by = c("Year", "Cl_new"))
  
  # Step 4: Add them to your data with Region = "Empty"
  waffle_data_complete <- waffle_data %>%
    bind_rows(
      missing_clusters %>%
        mutate(Region = "Empty", n = 81, pct = 100)
    )
  
  # Optional: Add gray color
  palette_regions[["Empty"]] <- "grey90"
  
  library(waffle)
  
  plt_waffle = ggplot(waffle_data_complete, aes(fill = Region, values = pct)) +
    geom_waffle(n_rows = 10, size = 0, color = "transparent") +
    facet_grid(Year ~ Cl_new) +
    scale_fill_manual(values = palette_regions) +
    coord_equal() +
    theme(
      panel.spacing.x = unit(0.1, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )
  
  plt_waffle_counts = ggplot(waffle_data_complete, aes(fill = Region, values = n)) +
    geom_waffle(n_rows = 9, size = 0, color = "transparent") +
    facet_wrap(~Cl_new)+
    facet_grid(Year ~ Cl_new) +
    scale_fill_manual(values = palette_regions) +
    coord_equal() +
    theme(
      panel.spacing.x = unit(0.1, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )
  
  ggsave(filename = "waffle_regions.pdf",
         plot = plt_waffle,
         path = IMGFOLDER,
         width = 8, height = 8)
  
  ggsave(filename = "waffle_regions_counts.pdf",
         plot = plt_waffle_counts,
         path = IMGFOLDER,
         width = 8, height = 8)
  
}


# Network plot by PPE(coord)+sex(shape)+region(color) ----
source("src/tHMM_networkplots.R")
## Pairwise hamming distances between rows ----
HDM = array(NA, dim = c(dim(ghe.array)[1], dim(ghe.array)[1], dim(ghe.array)[3]))
for (t in 1:dim(ghe.array)[3]) {
  for (i in 1:(dim(ghe.array)[1])) {
    for (j in i:dim(ghe.array)[1]) {
      HDM[i, j, t] = sum(ghe.array[i, , t] != ghe.array[j, , t])
      HDM[j, i, t] = HDM[i, j, t]
    }
  }
}
dimnames(HDM) = list("Country1" = dimnames(ghe.array)[[1]], 
                     "Country2" = dimnames(ghe.array)[[1]],
                     "Year" = dimnames(ghe.array)[[3]])
## Hamming similarity matrix (HSM) ----
HSM = 19 - HDM

## Network plot ----
plt_net = list()
for (t in 1:22){
  yyyy = years[t]
  tmp = plotgraph(HSM[,,t], 
                  thres.plot = 14, thres.pos = 17, 
                  comm = PPE[,t], weight.comm = 10,
                  positions = NULL)
  
  plt_net[[t]] = tmp
  
  if (SAVE){
    if (!dir.exists(paste0(IMGFOLDER, "network/"))){
      dir.create(paste0(IMGFOLDER, "network/"), recursive = TRUE)
    }
    
    width = 7; height = 5.5; zoom = 1
    ggsave(plot = tmp,
           filename = paste0("network/net_", yyyy, ".pdf"),
           path = IMGFOLDER,
           height = height*zoom, width = width*zoom)
  }
  
}


# Flow plot ----
plt_flow = plot_tHMM_flow(ppe = PPE.custord, tosort = F)

if (SHOW) plt_flow

if (SAVE){
  width = 12; height = 8; zoom = 1
  filename = "flow.pdf"
  ggsave(filename = filename,
         plot = plt_flow,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}


# Matrix of NMI between estimated partitions ----
NMI = matrix(NA, nrow = nyears, ncol = nyears)
for (i in 1:nyears) {
  for (j in 1:nyears) {
    NMI[i, j] = aricode::NMI(PPE.seqsort[, i], PPE.seqsort[, j])
  }
}
rownames(NMI) = years
colnames(NMI) = years

plt_NMI = NMI %>% data.frame() %>%
  rename_with(~ sub("^X", "", .), .cols = everything()) %>% 
  rownames_to_column("Year1") %>% 
  pivot_longer(-Year1, names_to = "Year2", values_to = "NMI") %>% 
  ggplot(aes(x = Year1, y = Year2, fill = NMI)) +
  geom_tile() +
  scale_fill_viridis(limits = c(0, 1)) +
  labs(x = "Year1", y = "Year2", fill = "NMI") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.title = element_blank(),
        text = element_text(size = 10),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.box.margin = margin(t = 0, r = -5, b = 0, l = -8))

if (SHOW) plt_NMI + geom_text(aes(label = round(NMI, 2)), size = 2.5, color = "black")

if (SAVE){
  width = 5; height = 4; zoom = 0.75
  filename = "NMI.pdf"
  ggsave(filename = filename,
         plot = plt_NMI,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}


# Estimates of sigma ----
param_sample = readRDS(paste0(OUTFOLDER, "param_", NAME, ".RDS"))
sigma_sample = param_sample$sigma
sigma_summary = lapply(sigma_sample, 
                       function(S) 
                       {tmp = apply(S[ , , -(1:10000)], 1:2, function(x) {
                         HDI = hdi(x, credMass = 0.95) %>% unname()
                         c("mean" = median(x), "sd" = sd(x), HDI_l = HDI[1], HDI_u = HDI[2])
                       })
                       tmp
                       })
omega_summary = lapply(sigma_sample, 
                       function(S) 
                       {tmp = apply(exp(-1/S[ , , -(1:10000)]), 1:2, function(x) {
                         HDI = hdi(x, credMass = 0.95) %>% unname()
                         c("mean" = mean(x), "sd" = sd(x), HDI_l = HDI[1], HDI_u = HDI[2])
                       })
                       tmp
                       })
for (t in 1:22){
  Cl_names = df_clusters %>% filter(Year == years[t]) %>% arrange(Cl_old) %>% pull(Cl_new)
  dimnames(sigma_summary[[t]]) = list("Stat" = c("mean", "sd", "HDI_l", "HDI_u"),
                                      "Cl" = paste0("Cl", Cl_names),
                                      "Age" = ages)
  dimnames(omega_summary[[t]]) = list("Stat" = c("mean", "sd", "HDI_l", "HDI_u"),
                                      "Cl" = paste0("Cl", Cl_names),
                                      "Age" = ages)
}
names(sigma_summary) = names(omega_summary) = years

sigma_df = sigma_summary %>% melt %>% rename(Year = L1) %>% mutate(Year = as.integer(Year)) %>% 
  pivot_wider(names_from = Stat, values_from = value) %>% 
  mutate(AgeInt = as.integer(Age)) %>% 
  left_join(df_clusters %>% select(Year, Cl_new, order) %>% unique() %>% 
              mutate(Cl = paste0("Cl", Cl_new)), 
            by = c("Year", "Cl"))

omega_df = omega_summary %>% melt %>% rename(Year = L1) %>% mutate(Year = as.integer(Year)) %>%
  pivot_wider(names_from = Stat, values_from = value) %>% 
  mutate(AgeInt = as.integer(Age)) %>% 
  left_join(df_clusters %>% select(Year, Cl_new, order) %>% unique() %>% 
              mutate(Cl = paste0("Cl", Cl_new)), 
            by = c("Year", "Cl"))

tmptheme = theme(axis.title = element_blank(),
                 axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                 legend.position = "inside",
                 legend.position.inside = c(5/6, 0.5/16),
                 legend.direction = "horizontal")

sigma_plot = sigma_df %>%
  ggplot(aes(x = AgeInt, y = mean)) +
  geom_point(aes(color = Cl), alpha = 0.4, size = 0.9, show.legend = F) +
  stat_summary(geom = "line", fun = median, col = "black", size = 0.4, show.legend = FALSE) +
  stat_summary(geom = "line", fun = min, col = "black", size = 0.25, lty = 2, show.legend = FALSE) +
  stat_summary(geom = "line", fun = max, col = "black", size = 0.25, lty = 2, show.legend = FALSE) +
  scale_x_continuous(breaks = seq(1, nages, by = 3), labels = ages[seq(1, nages, by = 3)]) +
  facet_wrap(~ Year, scales = "fixed", nrow = 4) +
  scale_shape_manual("", values = 16, labels = TeX(r"($\hat{\sigma}_{kxt}$)")) +
  guides(shape = guide_legend(override.aes = list(size = 3, alpha = 0.5)),
         linetype = guide_legend(title = "", override.aes = list(lwd = 0.8))) +
  # Fake points and lines for the legend
  geom_point(data = data.frame(x = 2, y = mean(sigma_df$mean), type = "Cluster"), aes(x = x, y = y, shape = type), alpha = 0, size = 0) +
  geom_line(data = data.frame(x = rep(2, 6), y = rep(0.2, 6), meas = c("Median", "Minimum/Maximum")), 
            aes(x = x, y = y, linetype = meas), size = 0) +
  tmptheme


if (SHOW) sigma_plot

omega_plot = omega_df %>%
  ggplot(aes(x = AgeInt, y = mean)) +
  geom_point(aes(color = Cl), alpha = 0.5, size = 0.9, show.legend = F) +
  stat_summary(geom = "line", fun = median, col = "black", size = 0.4, show.legend = FALSE) +
  stat_summary(geom = "line", fun = min, col = "black", size = 0.25, lty = 2, show.legend = FALSE) +
  stat_summary(geom = "line", fun = max, col = "black", size = 0.25, lty = 2, show.legend = FALSE) +
  scale_x_continuous(breaks = seq(1, nages, by = 3), labels = ages[seq(1, nages, by = 3)]) +
  facet_wrap(~ Year, scales = "fixed", nrow = 4) +
  # Fake points and lines for the legend
  geom_point(data = data.frame(x = 2, y = mean(omega_df$mean), type = "Cluster"), aes(x = x, y = y, shape = type), alpha = 0, size = 0) +
  geom_line(data = data.frame(x = rep(2, 6), y = rep(0.2, 6), meas = c("Median", "Minimum/Maximum")), 
            aes(x = x, y = y, linetype = meas), size = 0) +
  scale_shape_manual("", values = 16, labels = TeX(r"($\hat{\omega}_{kxt}$)")) +
  guides(shape = guide_legend(override.aes = list(size = 3, alpha = 0.5)),
         linetype = guide_legend(title = "", override.aes = list(lwd = 0.5))) +
  tmptheme

if (SHOW) omega_plot


if (SAVE){
  width = 12; height = 7; zoom = 0.85
  filename = "sigma.pdf"
  ggsave(filename = filename,
         plot = sigma_plot,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
  
  filename = "omega.pdf"
  ggsave(filename = filename,
         plot = omega_plot,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}

## Omega heatmap ----
# Complete imputing "fake" clusters
all_clusters = unique(center_df$Cl)
omega_df.complete = omega_df %>% 
  mutate(order = as.integer(as.character(order))) %>% 
  complete(Year, Age, Cl = as.character(all_clusters)) %>%
  # Impute order also for "fake" clusters
  group_by(Year) %>% 
  mutate(order = ifelse(is.na(order),
                        max(order, na.rm = TRUE) + dense_rank(Cl[is.na(order)]),
                        order)) %>%
  ungroup() %>% 
  mutate(order = factor(order, levels = 484:1))

omega_heatmap = omega_df.complete %>% 
  ggplot(aes(x = AgeInt, y = order)) +
  geom_raster(aes(fill = mean)) +
  scale_fill_viridis_c(TeX(r"($\hat{\omega}_{kxt}$)"),
                       trans = "log", breaks = c(0.002, 0.02, 0.13)) +
  # scale_fill_viridis_c(limits = c(0, max(omega_df$mean))) +
  scale_x_continuous(breaks = seq(1, nages, by = 3), labels = ages[seq(1, nages, by = 3)]) +
  facet_wrap(~ Year, scales = "free_y", nrow = 5) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 22),
        legend.position = "inside",
        legend.position.inside = c(3.5/5, 0.125/5),
        legend.direction = "horizontal",
        legend.text = element_text(angle = 270, size = 20, hjust = 1, vjust = 0.5),
        legend.key.height = unit(0.025, "npc"),  # height of the gradient bar
        legend.key.width = unit(0.05, "npc"))   # width of the gradient bar)


## Sigma heatmap ----
sigma_df.complete = sigma_df %>% 
  mutate(order = as.integer(as.character(order))) %>% 
  complete(Year, Age, 
           Cl = paste0("Cl", all_clusters)) %>%
  # Impute order also for "fake" clusters
  group_by(Year) %>% 
  mutate(order = ifelse(is.na(order),
                        max(order, na.rm = TRUE) + dense_rank(Cl[is.na(order)]),
                        order)) %>%
  ungroup() %>% 
  mutate(order = factor(order, levels = 484:1))

sigma_heatmap = sigma_df.complete %>% 
  ggplot(aes(x = AgeInt, y = order)) +
  geom_raster(aes(fill = mean)) +
  scale_fill_viridis_c(TeX(r"($\hat{\sigma}_{kxt}$)"),
                       trans = "log", 
                       breaks = c(0.15, 0.25, 0.65)) +
  scale_x_continuous(breaks = seq(1, nages, by = 3), labels = ages[seq(1, nages, by = 3)]) +
  facet_wrap(~ Year, scales = "free_y", nrow = 5) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(3.5/5, 0.125/5),
        legend.direction = "horizontal",
        legend.text = element_text(angle = 270, hjust = 1, vjust = 0.5),
  )

if (SHOW) sigma_heatmap

# sigma_bar = sigma_df %>% 
#   mutate(Group = if_else(Age== "0-1", "Perinatal", "Others")) %>% 
#   group_by(Year, Group) %>% 
#   summarise(Min = min(mean),
#             Median = median(mean) #,
#             # Max = max(mean)
#   ) %>%
#   pivot_longer(cols = c(Min, Median),#, max), 
#                names_to = "statistic", 
#                values_to = "value") %>%
#   ggplot() + 
#   geom_bar(aes(x = Year, y = value), stat = "identity") +
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
#   facet_grid(statistic~Group, scales = "free_y") +
#   theme(axis.title = element_blank())

sigma_df.bar = sigma_df %>% 
  mutate(Group = if_else(Age== "0-1", "Perinatal", "Others")) %>% 
  group_by(Year, Group) %>% 
  summarise(Min = min(mean),
            Median = median(mean) #,
            # Max = max(mean)
  ) %>%
  pivot_longer(cols = c(Min, Median),#, max), 
               names_to = "statistic", 
               values_to = "value") %>%
  mutate(Group = factor(Group, levels = c("Perinatal", "Others")))

sigma_bar = ggplot() + 
  geom_bar(data = sigma_df.bar %>% filter(Group == "Perinatal"),
           mapping = aes(x = Year, y = value, fill = "Perinatal"), 
           alpha = 0.5, stat = "identity") +
  geom_bar(data = sigma_df.bar %>% filter(Group == "Others"),
           mapping = aes(x = Year, y = value, fill = "Others"),
           alpha = 0.5, stat = "identity") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  facet_grid(rows = vars(statistic), scales = "free_y") +
  theme(axis.title = element_blank())

if (SHOW) sigma_bar

# Put sigma heatmap and geombar in same plot
if (SAVE){
  plt_out = sigma_heatmap + sigma_bar + plot_layout(ncol = 1, heights = c(4, 1))
  
  width = 12; height = 12; zoom = 1
  ggsave(plot = plt_out,
         filename = "sigmaheatbar.pdf",
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
  width = 9; height = 6; zoom = 0.75
  ggsave(plot = sigma_heatmap,
         filename = "sigmaheat.pdf",
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}


# STROKE EU clusters ----

# Plot ----
# cause_breaks = ghe_PPE %>% 
#   filter(Year == 2000 & Cl_new == 8 | Year == 2019 & Cl_new %in% c(7, 15)) %>% 
#   group_by(CauseS) %>% 
#   summarise(n = n()) %>% 
#   # arrange(desc(n))
#   filter(n >= 10) %>% pull(CauseS) %>% as.character()

plt = ghe_PPE %>% 
  filter(Year == 2000 & Cl_new == 8 | Year == 2019 & Cl_new %in% c(7, 15)) %>%
  arrange(Cl_new, Region, Sex) %>%
  mutate(ID = factor(ID, levels = unique(ID))) %>% 
  ggplot(aes(x = Age, y = ID, fill = CauseS)) +
  geom_tile() +
  ggh4x::facet_nested(
    rows = vars(Year, Cl_new),
    scales = "free_y", space = "free_y"
  ) +
  scale_fill_manual("Cause", values = palette) +
  guides(fill = guide_legend(override.aes = list(color = "grey30"))) +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right")



data_aux = unique(ghe_PPE %>% 
                    filter(Year == 2000 & Cl_new == 8 | Year == 2019 & Cl_new %in% c(7, 15)) %>%
                    select(Year, ID, Cl_new, Sex, Region)) %>%
  arrange(Region, Sex) %>%
  mutate(ID = factor(ID, levels = unique(ID))) %>%
  pivot_longer(cols = -c(Year:Cl_new), , names_to = "x", values_to = "fill") %>% 
  mutate(x = factor(x, levels = c("Region", "Sex")))


# colors_aux = c("lightblue", "pink", "#EF476F", "#F78C6B", "#FFD166", "#06D6A0", "#118AB2", "#073B4C")
# names(colors_aux) = c("M", "F", "AF", "AM", "EM", "EU", "SEA", "WP")
colors_aux = c("M"= "lightblue", "F" = "pink", palette_regions)

plt_aux = ggplot(data = data_aux) +
  geom_tile(aes(x = x, y = ID, fill = fill)) +
  ggh4x::facet_nested(
    rows = vars(Year, Cl_new),
    scales = "free_y", space = "free_y"
  ) +
  scale_fill_manual("Region", values = colors_aux, breaks = unique(ghe_PPE$Region), 
                    labels = unique(ghe_PPE$Region.short)) +
  scale_y_discrete(breaks = c("Georgia_F", "Serbia_F", "Portugal_F", "Montenegro_F", "North Macedonia_F",
                              "Albania_F", "Romania_F", "North Macedonia_M")) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(color = "grey30"))) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        legend.position = "left")

plt_out = (plt_aux + plt + plot_layout(widths = c(1, 20),  
                                       guides = "keep"))
# &
#   theme(legend.position = "bottom",
#         legend.justification.bottom = "center")



if (SHOW) plt_out

if (SAVE){
  width = 13; height = 6.5; zoom = 1
  ggsave(plot = plt_out,
         filename = "strokeEU_clusters.pdf",
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}




## NMI between "stroke" clusters in 2000 and 2019 ----
PPEstroke = full_join(ghe_PPE %>% 
                        filter(Year == 2000, Cl_new %in% strok_cl[[1]]) %>%
                        select(ID, Cl_new) %>% unique(),
                      ghe_PPE %>% 
                        filter(Year == 2019, Cl_new %in% strok_cl[[2]]) %>%
                        select(ID, Cl_new) %>% unique(),
                      by = "ID") %>% 
  mutate(Cl2000 = if_else(is.na(Cl_new.x), 0, Cl_new.x),
         Cl2019 = if_else(is.na(Cl_new.y), 0, Cl_new.y))

table(PPEstroke$Cl2000, PPEstroke$Cl2019)

# Diabetes clusters ----
diabetes_clusters = center_df %>% 
  group_by(Year, Cl_new) %>%
  mutate(Diabetes = if_else(sum(Cause == "Diabetes") > 3, 
                            "Yes", "No")) %>%
  select(Year, Cl_new, Diabetes) %>% unique()

# cause_breaks = ghe_PPE %>%
#   left_join(diabetes_clusters, by = c("Year", "Cl_new")) %>%
#   filter(Diabetes == "Yes") %>%
#   group_by(CauseS) %>% 
#   summarise(n = n()) %>% 
#   filter(n >= 10) %>% pull(CauseS) %>% as.character()

plt = ghe_PPE %>% 
  left_join(diabetes_clusters, by = c("Year", "Cl_new")) %>%
  filter(Diabetes == "Yes") %>% 
  arrange(Region, Sex) %>%
  mutate(ID = factor(ID, levels = unique(ID))) %>%
  filter(Year %in% c(2000, 2009, 2018)) %>% 
  ggplot(aes(x = Age, y = ID, fill = CauseS)) +
  geom_tile() +
  ggh4x::facet_nested(
    rows = vars(Year, Cl_new),
    scales = "free_y", space = "free_y"
  ) +
  scale_fill_manual("Cause", values = palette) +
  guides(fill = guide_legend(override.aes = list(color = "grey30"))) +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right", legend.justification.right = "top")

data_aux = unique(ghe_PPE %>% 
                    left_join(diabetes_clusters, by = c("Year", "Cl_new")) %>%
                    filter(Diabetes == "Yes") %>%
                    select(Year, ID, Cl_new, Sex, Region)) %>%
  arrange(Region, Sex) %>%
  mutate(ID = factor(ID, levels = unique(ID))) %>%
  pivot_longer(cols = -c(Year:Cl_new), , names_to = "x", values_to = "fill")

colors_aux = c("M"= "lightblue", "F" = "pink", palette_regions)

plt_aux = ggplot(data = data_aux %>% filter(Year %in% c(2000, 2009, 2018))) +
  geom_tile(aes(x = x, y = ID, fill = fill)) +
  ggh4x::facet_nested(
    rows = vars(Year, Cl_new),
    scales = "free_y", space = "free_y"
  ) +
  scale_fill_manual("Region", values = colors_aux, breaks = unique(ghe_PPE$Region),
                    labels = unique(ghe_PPE$Region.short)) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(color = "grey30")))+
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    strip.text = element_blank(),
    legend.position = "left",
    legend.justification.right = "center")


plt_out = plt_aux + plt + plot_layout(widths = c(1, 20), guides = "keep")
# & 
#   theme(legend.position = "right",
#         legend.justification.bottom = "center")


if (SHOW) plt_out

if (SAVE){
  width = 12; height = 4; zoom = 1.1
  ggsave(plot = plt_out,
         filename = "diabetes_clusters.pdf",
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}



# SOMALIA (famine) ----
ghe %>% 
  filter(CountryN == "Somalia") %>%
  select(Year, CauseS, Sex, Age) %>%
  mutate(Age = factor(Age, levels = ages, ordered = TRUE)) %>%
  ggplot(aes(x = Age, y = Year, fill = CauseS)) +
  geom_tile() +
  scale_fill_manual(values = palette) +
  facet_grid(~Sex) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 12),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.position.inside = c(0.5, -0.1))

# SELF-HARM ----
## Number of self-harm per year ----

# Total by year
ghe %>% 
  group_by(Year) %>%
  summarise(Selfharm = sum(CauseS == "Self-harm")) %>%
  ungroup()

# Plot male vs female every year
ghe %>% 
  group_by(Year, Sex) %>%
  summarise(Selfharm = sum(CauseS == "Self-harm")) %>%
  ungroup() %>% 
  ggplot() +
  geom_line(aes(x = Year, y = Selfharm, color = Sex))

# Difference between male and female every year
ghe %>% 
  group_by(Year, Sex) %>%
  summarise(Selfharm = sum(CauseS == "Self-harm")) %>%
  ungroup() %>% 
  pivot_wider(names_from = "Sex", values_from = "Selfharm") %>% 
  mutate(diff = M - F)

# Plot male vs female every year and age
ghe %>% 
  group_by(Year, Age, Sex) %>%
  summarise(Selfharm = sum(CauseS == "Self-harm")) %>%
  ungroup() %>% 
  ggplot() +
  geom_line(aes(x = Year, y = Selfharm, color = Sex)) + facet_wrap(~Age, ncol = 6)

# Plot male vs female every year, age and region
ghe_PPE %>% 
  filter(Age %in% ages[4:11]) %>% 
  group_by(Year, Age, Region) %>%
  summarise(Selfharm = sum(CauseS == "Self-harm") / n()) %>%
  ungroup() %>% 
  ggplot() +
  geom_line(aes(x = Year, y = Selfharm, color = Region)) + facet_wrap(~Age) +
  scale_color_manual(values = palette_regions)

# Without Age
ghe_PPE %>% 
  group_by(Year, Region) %>%
  summarise(Selfharm = sum(CauseS == "Self-harm") / n()) %>%
  ungroup() %>% 
  ggplot() +
  geom_line(aes(x = Year, y = Selfharm, color = Region)) +
  scale_color_manual(values = palette_regions)


maxselfharm = ghe %>% 
  group_by(Year, Age) %>%
  summarise(Selfharm = sum(CauseS == "Self-harm")) %>% pull(Selfharm) %>% max

plt_nmbr_selfharm = 
  ghe %>% 
  group_by(Year, Age) %>%
  summarise(Selfharm = sum(CauseS == "Self-harm")) %>%
  ungroup() %>% 
  ggplot(aes(x = as.integer(Age), y = Selfharm)) +
  # Add three vertical shaded areas for three interval of x-axis: [3.5, 5.5], [5.5, 7.5], [7.5, 12.5]
  # geom_rect(data = data.frame(xmin = c(3.5, 5.5, 7.5), xmax = c(5.5, 7.5, 12.5),
  #                             phase = c("Teenagers", "Young Adult", "Adult")),
  #           aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = maxselfharm, fill = phase), alpha = 0.025) +
  geom_rect(aes(xmin = 3.5, xmax = 5.5, ymin = 0, ymax = maxselfharm,
                fill = "Teen"), alpha = 0.025) +
  geom_rect(aes(xmin = 5.5, xmax = 7.5, ymin = 0, ymax = maxselfharm, 
                fill = "YAd"), alpha = 0.025) +
  geom_rect(aes(xmin = 7.5, xmax = 10.5, ymin = 0, ymax = maxselfharm,
                fill = "Ad"), alpha = 0.025) +
  geom_line() + 
  geom_point(size = 0.5) + 
  scale_fill_manual(values = c("#60d394", "#ffd97d", "#ff9b85"),
                    name = "Age",
                    labels = c("Adolescence", "Early Adulthood", "Adulthood")) +
  scale_x_continuous(breaks = seq(1, nages, by = 3), 
                     labels = ages[seq(1, nages, by = 3)]) +
  scale_y_continuous(limits = c(0, maxselfharm)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3),
                             title.position = "left")) +
  facet_wrap(~Year, ncol = 6) + 
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "inside",
        legend.position.inside = c(5/6, 1/32))

if (SHOW) plt_nmbr_selfharm

if (SAVE){
  width = 12; height = 6; zoom = 1
  ggsave(filename = "selfharm_nmbr.pdf",
         plot = plt_nmbr_selfharm,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
  
  width = 12; height = 4; zoom = 1
  ggsave(filename = "selfharm_nmbr_short.pdf",
         plot = plt_nmbr_selfharm,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}


plt_nmbr_selfharm2 =
  ghe %>% 
  group_by(Year, Age) %>%
  summarise(Selfharm = sum(CauseS == "Self-harm")) %>%
  ungroup() %>% 
  ggplot(aes(x = as.integer(Age), y = Selfharm, color = Year, group = Year)) +
  geom_line() +
  geom_point(size = 0.5) +
  scale_color_viridis_c(direction = -1) + 
  scale_x_continuous(breaks = seq(1, nages, by = 3), 
                     labels = ages[seq(1, nages, by = 3)]) +
  scale_y_continuous(limits = c(0, maxselfharm)) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 12),
        legend.position = "right")


if (SHOW) plt_nmbr_selfharm2

if (SAVE){
  width = 12; height = 6; zoom = 1
  ggsave(filename = "selfharm_nmbr2.pdf",
         plot = plt_nmbr_selfharm2,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}


plt_nmbr_selfharm3 =
  ghe %>% 
  group_by(Year, Age) %>%
  summarise(Selfharm = sum(CauseS == "Self-harm")) %>%
  ungroup() %>% 
  ggplot(aes(x = as.integer(Age), y = Selfharm, color = Year, group = Year)) +
  geom_smooth(span = 0.25,
              se = FALSE,
              linewidth = 0.7) + 
  scale_color_viridis_c(direction = -1) + 
  scale_x_continuous(breaks = seq(1, nages, by = 3), 
                     labels = ages[seq(1, nages, by = 3)]) +
  scale_y_continuous(limits = c(0, maxselfharm)) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 12),
        legend.position = "right")

if (SAVE){
  width = 12; height = 6; zoom = 1
  ggsave(filename = "selfharm_nmbr3.pdf",
         plot = plt_nmbr_selfharm3,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}



listplt_nmbr_selfharm = lapply(years, function(y){
  ghe %>% 
    filter(Year == y) %>% 
    group_by(Year, Age) %>%
    summarise(Selfharm = sum(CauseS == "Self-harm")) %>%
    ungroup() %>% 
    ggplot(aes(x = as.integer(Age), y = Selfharm)) +
    geom_line() + geom_point(size = 0.5) + 
    scale_x_continuous(breaks = seq(1, nages, by = 3), 
                       labels = ages[seq(1, nages, by = 3)]) +
    scale_y_continuous(limits = c(0, maxselfharm)) +
    facet_wrap(~Year, ncol = 6) + 
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          text = element_text(size = 12),
          legend.position = "none")
})

## Plot centers of clusters with self-harm ----
center_df.selfharm = center_df %>% 
  group_by(Year, Cl_new) %>%
  mutate(Selfharm = if_else(sum(Cause == "Self-harm") >= 2,
                            "Yes", "No")) %>%
  filter(Selfharm == "Yes")

center_df %>% 
  group_by(Year, Cl_new) %>%
  mutate(Selfharm = if_else(sum(Cause == "Self-harm") >= 2, 
                            "Yes", "No")) %>%
  filter(Selfharm == "Yes") %>% 
  mutate(order = factor(order,
                        levels = center_df.selfharm %>% select(Year, order) %>% unique %>% arrange(Year, order) %>% pull(order),
                        labels = center_df.selfharm %>% select(Year, order, Cl_new) %>% unique %>% arrange(Year, order) %>% pull(Cl_new))) %>% 
  ggplot() +
  geom_tile(aes(x = as.integer(Age), y = order, fill = Cause), color = "grey30", linewidth = .05) +
  scale_fill_manual(values = palette, na.value = "transparent",
                    breaks = sort(names(palette))) +
  guides(fill = guide_legend(ncol = 2, override.aes = list(color = "grey30"))) +
  scale_x_continuous(breaks = seq(1, nages, by = 3), 
                     labels = ages[seq(1, nages, by = 3)]) +
  facet_wrap(~ Year, ncol = 5, scales = "free") +
  theme(
    legend.key.height = unit(0.025, "npc"), legend.key.width = unit(0.015, "npc"),
    legend.position = "right",
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.text.y = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    text = element_text(size = 12))




## Plot selfharm vs Sex and Region ----
selfharm_clusters = center_df %>% 
  group_by(Year, Cl_new) %>%
  mutate(Selfharm = if_else(sum(Cause == "Self-harm") >= 2, 
                            "Yes", "No")) %>%
  select(Year, Cl_new, Selfharm) %>% unique()

dfplot_selfharm = ghe_PPE %>%
  left_join(selfharm_clusters, by = c("Year", "Cl_new")) %>%
  filter(Selfharm == "Yes") %>%
  select(Year, Cl_new, Country, Sex, Region) %>% 
  unique() %>% 
  left_join(prop_df %>% select(Year, Cl_new), by = c("Year", "Cl_new")) %>% 
  mutate(Cl_new = factor(Cl_new, levels = unique(prop_df$Cl_new))) %>%
  group_by(Year, Cl_new, Sex) %>% mutate(nSex = n()) %>% ungroup() %>%
  group_by(Year, Cl_new, Region) %>% mutate(nRegion = n()) %>% ungroup() %>%
  select(Year, Cl_new, Region, Sex, nSex, nRegion) %>% unique %>%
  pivot_longer(cols = c(nSex, nRegion), names_to = "nType", values_to = "n")


nclmax = length(unique(dfplot_selfharm))
ymax =  dfplot_selfharm %>% filter(nType == "nSex") %>%  select(Year, Cl_new, Sex, n) %>% unique() %>% 
  group_by(Year, Cl_new) %>% summarise(n = sum(n)) %>% pull(n) %>% max

listplt_selfharm = lapply(years, function(yyyy){
  d = dfplot_selfharm %>%
    filter(Year == yyyy)
  
  ncl = length(unique(d$Cl_new))
  add = .15  * ncl/ nclmax
  width = .3 * ncl/ nclmax
  
  tmp = d %>% 
    ggplot(data = d, mapping = aes(y = n)) +
    geom_col(data = filter(d, nType == "nSex") %>% select(-Region) %>% unique(), 
             aes(x = dense_rank(Cl_new) - add, fill = Sex), width = width, col = "gray30", linewidth = 0.25) +
    scale_fill_manual(values = c("F" = "pink", "M" = "lightblue"), labels = c("Female", "Male"), guide = guide_legend(order = 1)) +
    scale_x_continuous(breaks = unique(dense_rank(d$Cl_new)), labels = unique(d$Cl_new)) +
    ggnewscale::new_scale_fill() +
    geom_col(data = filter(d, nType == "nRegion") %>% select(-Sex) %>% unique(), 
             aes(x = dense_rank(Cl_new) + add, fill = Region), width = width, col = "gray30") +
    scale_fill_manual(values = palette_regions, guide = guide_legend(order = 2),
                      labels = levels(ghe_PPE$Region.short)) +
    scale_y_continuous(breaks = seq(0, max(dfplot_selfharm$n), by = 10),
                       limits = c(0, ymax)) +
    facet_wrap(~Year) + 
    theme(text = element_text(size = 12),
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  if (yyyy != 2000){
    tmp = tmp + theme(legend.position = "none")
  } else {
    tmp = tmp + theme(legend.direction = "horizontal",
                      text = element_text(size = 12))
  }
  
  if (!(yyyy %in% c(2000, 2006, 2012, 2018))){
    tmp = tmp + theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank())
  } else {
    tmp = tmp
  }
  
})

lgnd = cowplot::get_legend(listplt_selfharm[[1]])
listplt_selfharm[[1]] = listplt_selfharm[[1]] + theme(legend.position = "none")

plt_out = ggpubr::ggarrange(plotlist = listplt_selfharm, ncol = 6, nrow = 4, legend = "none") +
  annotation_custom(lgnd, xmin = 4/6, xmax = 1, ymin = 0, ymax = 1/4)

if (SAVE){
  width = 12; height = 6; zoom = 1
  ggsave(plot = plt_out,
         filename = "selfharm_SexRegion.pdf",
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}

# Add curves of number of selfharm per age x year
plotlist <- c(listplt_selfharm, listplt_nmbr_selfharm)

# Layout: 4 rows x 12 columns = 48 slots
# We'll place the 44 plots and leave 4 slots as NA (blanks)
layout_matrix <- matrix(c(
  1:6,    23:28,
  7:12,   29:34,
  13:18,  35:40,
  19:22,  0, 0, 41:44, 0, 0
), nrow = 4, ncol = 12, byrow = TRUE)

plt_out = (deeptime::ggarrange2(plots = plotlist,
                                layout = layout_matrix,
                                widths = c(rep(1, 6), rep(1/2, 6)))) %>% ggplotify::as.ggplot() +
  annotation_custom(lgnd, xmin = 2/3-2/6, xmax = 2/3, ymin = 0, ymax = 1/4)

# if (SHOW) plt_out

if (SAVE){
  width = 18; height = 6; zoom = 1
  ggsave(plot = plt_out,
         filename = "selfharm_SexRegion2.pdf",
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}


# Plot of self-harm clusters ----
selfharm_plots = list()

for (yyyy in years){
  plt = ghe_PPE %>% 
    left_join(selfharm_clusters, by = c("Year", "Cl_new")) %>%
    filter(Selfharm == "Yes") %>% 
    arrange(Region, Sex) %>%
    mutate(ID = factor(ID, levels = unique(ID))) %>%
    filter(Year == yyyy) %>% 
    ggplot(aes(x = Age, y = ID, fill = CauseS)) +
    geom_tile() +
    ggh4x::facet_nested(
      rows = vars(Year, Cl_new),
      scales = "free_y", space = "free_y"
    ) +
    scale_fill_manual("Cause", values = palette) +
    guides(fill = guide_legend(override.aes = list(color = "grey30"))) +
    theme(
      legend.key.height = unit(0.25, "cm"), legend.key.width = unit(0.25, "cm"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      text = element_text(size = 12),
      legend.position = "right", legend.justification.right = "top")
  
  data_aux = unique(ghe_PPE %>% 
                      left_join(selfharm_clusters, by = c("Year", "Cl_new")) %>%
                      filter(Selfharm == "Yes") %>%
                      select(Year, ID, Cl_new, Sex, Region)) %>%
    filter(Year == yyyy) %>% 
    arrange(Region, Sex) %>%
    mutate(ID = factor(ID, levels = unique(ID))) %>%
    pivot_longer(cols = -c(Year:Cl_new), , names_to = "x", values_to = "fill")
  
  colors_aux = c("M"= "lightblue", "F" = "pink", palette_regions)
  
  plt_aux = ggplot(data = data_aux) +
    geom_tile(aes(x = x, y = ID, fill = fill)) +
    ggh4x::facet_nested(
      rows = vars(Year, Cl_new),
      scales = "free_y", space = "free_y"
    ) +
    scale_fill_manual("Region", values = colors_aux, breaks = unique(ghe_PPE$Region),
                      labels = unique(ghe_PPE$Region.short)) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(color = "grey30")))+
    theme(
      legend.text = element_text(size = 8), legend.title = element_text(size = 8),
      legend.key.height = unit(0.25, "cm"), legend.key.width = unit(0.25, "cm"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      # axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      strip.text = element_blank(),
      text = element_text(size = 12),
      legend.position = "none")
  
  plt_tmp = plt_aux + plt + plot_layout(widths = c(1, 20), guides = "keep") +
    theme(legend.position = "none")
  selfharm_plots[[yyyy - 2000 + 1]] = plt_tmp
}


# & 
#   theme(legend.position = "right",
#         legend.justification.bottom = "center")
plt_out = selfharm_plots %>% 
  wrap_plots(ncol = 4, guides = "collect", 
             heights = ) +
  plot_annotation(title = "Self-harm clusters across years",
                  theme = theme(legend.position = "none"))

# if (SHOW) plt_out

if (SAVE){
  width = 30; height = 60; zoom = 1
  ggsave(plot = plt_out,
         limitsize = F,
         filename = "selfharm_clusters.pdf",
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}


# Table of CoD ----
# transform in latex code for a table
library(knitr)
library(kableExtra)
ghetable = ghe %>%  select(CauseT, CauseS, Age) %>% unique() %>% 
  mutate(Age = factor(Age, ordered = T)) %>% 
  group_by(CauseT, CauseS) %>%
  # get only the first and last age
  summarise(Agemin = min(Age), Agemax = max(Age),
            Age = paste0(min(Age), " to ", max(Age)), .groups = "drop") %>%
  arrange(Agemin, Agemax, CauseS) %>% select(-(Agemin:Agemax))

ghetable %>% 
  select(CauseS, CauseT, Age) %>% 
  arrange(CauseS) %>% 
  kable(caption = "Causes of Death",
        format = "latex",
        , booktabs = TRUE,
        label = "tab:causes_of_death",
        col.names = c("Cause", "Abbreviation", "Age Range"),
        row.names = FALSE,
        align = "l",
        linesep = "") %>%           # <-- suppress \addlinespace
  kable_styling(latex_options = "hold_position") %>%
  kable_styling() %>%
  row_spec(0, bold = TRUE) %>%           # header row
  row_spec(nrow(ghetable), hline_after = TRUE)


# PCA on PPE ----
M = do.call(cbind, apply(PPE[ , 1:5], 2, vec2mat))
# library(logisticPCA)
pca = logisticPCA(M, k = 15, m = 0, main_effects = F)
U = pca$U
rownames(U) = paste0("Cl", 1:nrow(U))
annotations = data.frame(year = as.character(do.call(c, lapply(1:5, function(t) 
  rep(years[t], max(PPE[,t]))))))
rownames(annotations) = rownames(U)

pheatmap::pheatmap(U, cluster_rows = T, treeheight_row = 0, treeheight_col = 0, 
                   annotation_row = annotations, clustering_method = "ward.D2")
pca$PCs
