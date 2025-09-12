# Want to compare the PPE from all the run also with those obtained by independent AP models

# Packages ----
options(warning = 1)
suppressPackageStartupMessages(
  {
    library(tidyverse)
    library(patchwork)
    library(reshape2)
  }
)
source("src/utils.R")

# Folders and flags
OUTFOLDER.tHMM = "output/tHMM/"
SHOW = TRUE
SAVE = TRUE
if (SAVE){
  IMGFOLDER = "img/tHMM/comparisons/"
  if (!dir.exists(IMGFOLDER)){
    dir.create(IMGFOLDER, recursive = TRUE)
  }
}

# Load ----
## Data ---- 
{
  ghe = readRDS("data/rds/dataGHE.RDS") %>%
    unite("ID", c(CountryN, Sex), remove = FALSE)
  ghe.lab = unique(ghe %>% select(c("Age", "CauseS"))) %>% 
    group_by(Age) %>% 
    mutate(CauseI = as.integer(factor(CauseS)))
  ghe.w = ghe %>% 
    left_join(ghe.lab, by = c("Age", "CauseS")) %>% 
    select(c("ID", "Age", "Year", "CauseI"))  %>% 
    arrange(ID, Age, Year)   %>% 
    pivot_wider(names_from = Age, values_from = CauseI)
  Y = abind::abind(map(split(ghe.w, ghe.w$Year), \(x) x %>% select(-Year) %>% column_to_rownames("ID")), along = 3)
}

## Define dimensions ----
years = dimnames(Y)[[3]]
IDs = dimnames(Y)[[1]]
ncauses = length(unique(ghe$CauseS))

## PPE and PSM ----
{
  PPE_tHMM = lapply(list.files(path = OUTFOLDER.tHMM, pattern = "PPE_.+\\.RDS", full.names = TRUE),
                    readRDS)
  PSM_tHMM = lapply(list.files(path = OUTFOLDER.tHMM, pattern = "PSM_.+\\.RDS", full.names = TRUE),
                    readRDS)
  
  PPE_AP = readRDS("output/AP_indmod/ppe.RDS")
  dimnames_AP = dimnames(PPE_AP)
  PPE_AP = seqsort(PPE_AP)
  PSM_AP = lapply(readRDS("output/AP_indmod/psm.RDS"), 
                  function(x) {
                    dimnames(x) = list("ID1" = IDs, "ID2" = IDs)
                    x})
  names(PSM_AP) = names(PSM_tHMM[[1]])
}

### Join the 6 PPEs and the 6 PSM together ----
{
  PPE = append(list(PPE_AP), PPE_tHMM)
  names(PPE) = c("AP", 
                 list.files(path = OUTFOLDER.tHMM, pattern = "PPE_.+\\.RDS", full.names = FALSE) %>% 
                   gsub(".RDS", "", .) %>% gsub("PPE_", "", .))
  PSM = list()
  PSM[[1]] = PSM_AP
  PSM[2:8] = PSM_tHMM
  names(PSM) = names(PPE)
}

# Sort rows and columns in the same order
for (j in seq_along(PPE)) {
  PPE[[j]] = PPE[[j]][order(rownames(PPE[[j]])), order(colnames(PPE[[j]]))]
  for (t in seq_along(PSM[[j]])) {
    PSM[[j]][[t]] = PSM[[j]][[t]][order(rownames(PSM[[j]][[t]])), order(colnames(PSM[[j]][[t]]))]
  }
}


# Palette causes ----
palette = c(pals::alphabet(26)[-(4:5)], pals::alphabet2(26), pals::glasbey(32))[1:ncauses]
names(palette) =  sort(table(ghe$CauseS %>% as.character()), decreasing = T) %>% names


# Compare PPEs ----
## NMI and ARI ----
NMI = ARI = array(NA, dim = c(length(PPE), length(PPE), dim(PPE[[1]])[2]),
                  dimnames = list("Run1" = names(PPE), 
                                  "Run2" = names(PPE), 
                                  "Year" = colnames(PPE[[1]])))


for (i in seq_along(PPE)) {
  for (j in seq_along(PPE)) {
    for (t in seq_len(dim(PPE[[1]])[2])) {
      NMI[i, j, t] = aricode::NMI(PPE[[i]][, t], PPE[[j]][, t])
    }
  }
}

plt_NMI = NMI %>% reshape2::melt(value.name = "NMI") %>% 
  ggplot() +
  geom_raster(aes(x = Run1, y = Run2, fill = NMI)) +
  scale_fill_viridis_c(limits = c(0, 1)) +
  geom_text(aes(x = Run1, y = Run2, label = round(NMI, 2)), size = 3) +
  facet_wrap(~ Year, ncol = 5) +
  labs(x = "Run 1", y = "Run 2", fill = "NMI") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank()
  )

if (SHOW) plt_NMI

if (SAVE){
  width = 12; height = 10; zoom = 1.25
  filename = "NMI.pdf"
  ggsave(plot = plt_NMI,
         filename = filename,
         path = IMGFOLDER,
         height = height * zoom, width = width * zoom)
}


# Plot of data with PPE from AP and Gnedin20019 ----
ppe.df = PPE$AP %>% reshape2::melt(varnames = c("ID", "Year"), value.name = "ClAP") %>% 
  left_join(PPE$Gnedin20019 %>% reshape2::melt(varnames = c("ID", "Year"), value.name = "ClGnedin"),
            by = c("ID", "Year"))
ghe_ppe = left_join(ghe, ppe.df, by = c("ID", "Year"))

for (yyyy in years){
  cat(yyyy, "\t")
  plt = ghe_ppe %>% 
    filter(Year == yyyy) %>%
    ggplot(aes(x = Age, y = ID, fill = CauseS)) +
    geom_raster() + 
    ggh4x::facet_nested(rows = vars(ClAP, ClGnedin), scales = "free_y", space = "free_y") +
    # facet_grid(rows = vars(Cl.ss), scales = "free_y", space = "free_y") + 
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
  
  data_aux = unique(ghe_ppe %>% filter(Year == yyyy) %>%
                      select(ID, ClAP, ClGnedin, Sex, Region)) %>% 
    pivot_longer(cols = -c(ID, ClAP, ClGnedin), , names_to = "x", values_to = "fill") %>% 
    mutate(x = factor(x, levels = c("Sex", "Region")))
  
  colors_aux = c("lightblue", "pink", "#EF476F", "#F78C6B", "#FFD166", "#06D6A0", "#118AB2", "#073B4C")
  names(colors_aux) = c("M", "F", "AF", "AM", "EM", "EU", "SEA", "WP")
  
  plt_aux = ggplot(data = data_aux) + 
    geom_tile(aes(x = x, y = ID, fill = fill)) +
    ggh4x::facet_nested(rows = vars(ClAP, ClGnedin), scales = "free_y", space = "free_y") +
    # facet_grid(rows = vars(ClAP), scales = "free_y", space = "free_y") + 
    scale_fill_manual("Region", values = colors_aux) +
    geom_text(aes(x = x, y = ID, label = fill), size = 2) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())
  
  plt_out = plt_aux + plt + plot_layout(widths = c(1, 20))
  
  ggsave(plot = plt_out,
         filename = paste0("nested_heatmaps_Gnedin20019/nested_heatmap_", yyyy, ".pdf"),
         path = IMGFOLDER,
         height = 40, width = 16)
}


# Load full RDS files to compare runs as for number of clusters and alpha trajectories ----
{
  listRDS = lapply(list.files(path = OUTFOLDER.tHMM, pattern = "res_.+\\.RDS", 
                              full.names = TRUE),
                   readRDS)
  listC = lapply(listRDS, function(x) x$output$C)
  listAlpha = lapply(listRDS, function(x) x$output$alpha)
  
  nburn = listRDS[[1]]$input$ctr$ctr_mcmc$nburnin
  niter = nburn + listRDS[[1]]$input$ctr$ctr_mcmc$nchain
}


ncl_trace = lapply(listC, function(C) apply(C, c(2, 3), function(x) length(unique(x)))) %>% 
  abind::abind(ncl_trace, along = 3)
dimnames(ncl_trace) = list("Year" = years, "Iteration" = 1:niter, "Run" = names(PPE)[-1])
ncl_PPEAP = apply(PPE_AP, 2, max) %>% as.data.frame() %>% rownames_to_column("Year") %>% 
  rename(nclAP = ".") %>% mutate(Year = as.integer(Year), Run = "AP_PPE")
ncl_PPEGnedin20019 = apply(PPE$Gnedin20019, 2, max) %>% as.data.frame() %>% 
  rownames_to_column("Year") %>% rename(nclGnedin = ".") %>% 
  mutate(Year = as.integer(Year), Run = "Gnedin20019_PPE")
plt_trcncl = ncl_trace %>% melt(value.name = "ncl") %>% 
  filter(Iteration > 100) %>%
  ggplot(aes(x = Iteration, y = ncl, col = Run)) +
  geom_line() +
  geom_hline(data = ncl_PPEAP, aes(yintercept = nclAP, col = Run), 
             linetype = "dashed") +
  geom_hline(data = ncl_PPEGnedin20019, aes(yintercept = nclGnedin, col = Run), 
             linetype = "dashed") +
  geom_vline(data = data.frame(nburn = nburn, Year = as.integer(years), Run = "burn-in"),
             aes(xintercept = nburn, col = Run), linetype = 3) +
  scale_color_manual(values = c(scales::hue_pal()(6), "black", "blue", "red"),
                     name = "Run") +
  scale_x_continuous(breaks = c(101, seq(5000, niter, by = 5000))) +
  facet_wrap(~ Year, ncol = 8) +
  labs(x = "Iteration", y = "Number of clusters") +
  theme(legend.position = "bottom")


if (SHOW) plt_trcncl

if (SAVE){
  width = 18; height = 6; zoom = 1.25
  filename = "ncl_trace.pdf"
  ggsave(plot = plt_trcncl,
         filename = filename,
         path = IMGFOLDER,
         height = height * zoom, width = width * zoom)
}  
