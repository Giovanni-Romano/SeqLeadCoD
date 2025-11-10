suppressPackageStartupMessages(library(tidyverse))
library(patchwork)
library(DescTools) # for Entropy
library(mcclust)
# library(mcclust.ext)
library(salso)
source('src/ArgientoPaci/code/complement_functions.R')
source("src/utils.R")

listRDS = lapply(list.files(path = "output/AP_indmod", pattern = "APout_.*\\.rds", full.names = T),
                 readRDS)[order(list.files(path = "output/AP_indmod", pattern = "APout_.*\\.rds", full.names = F) %>% 
                                  gsub(".rds", "", .) %>% str_split_i("_", 2) %>% as.numeric())]
ghe = readRDS("data/rds/dataGHE.RDS") %>% unite("ID", c(CountryN, Sex), remove = FALSE)
ages = levels(ghe$Age)
years = 2000:2021
IDs = ghe %>% pull(ID) %>% unique %>% sort
nIDs = length(IDs)
nyears = length(years)
ncauses = length(levels(ghe$CauseS))
nages = length(ages)

SHOW = TRUE
SAVE = FALSE

# Palette causes ----
palette = c(pals::alphabet(26)[-(4:5)], pals::alphabet2(26), pals::glasbey(32))[1:ncauses]
names(palette) =  sort(table(ghe$CauseS %>% as.character()), decreasing = T) %>% names


# Traceplots number of clusters ----
ncl_trace = lapply(listRDS, function(RDS) apply(RDS$sim_ghe$C, 1, function(x) length(unique(x))) %>% 
                     as.data.frame() %>% rownames_to_column("Iter") %>% rename(ncl = "."))
# names(ncl_trace) = years
# plt_ncl = ncl_trace %>% bind_rows(.id = "Year") %>% 
#   mutate(Year = as.numeric(Year), Iter = as.numeric(Iter)) %>% 
#   ggplot(aes(x = Iter, y = ncl)) +
#   geom_line() +
#   facet_wrap(~ Year) +
#   labs(title = "Traceplot number of clusters",
#        x = "Iteration", y = "Number of clusters") +
#   theme(text = element_text(size = 7))

if (SAVE){
  pdf(file = "img/AP_results/traceplots_ncl.pdf", height = 8, width = 12)
  par(mfrow = c(4, 6))
  for (y in 1:nyears){
    plot(ncl_trace[[y]]$Iter[-1], ncl_trace[[y]]$ncl[-1], type = "l", xlab = "Iteration", ylab = "Number of clusters",
         main = paste0("Traceplot - year=", years[y]))
  }
  dev.off()
}

burn = 1:5002

# Posterior distributions nmbr of clusters ----
post_k = list()
prior_k = list()
for (j in seq_along(listRDS)){
  RDS = listRDS[[j]]
  post_k.tmp = prop.table(table(RDS$sim_ghe$k[-burn]))
  prior_k.tmp = AntMAN::AM_prior_K_Pois(n=nIDs, RDS$param$gam, Lambda = RDS$param$Lambda)
  post_k[[j]] = post_k.tmp
  prior_k[[j]] = prior_k.tmp
}

if (SAVE){pdf("img/AP_results/posterior_nmbr_clusters.pdf", width = 9, height = 6)}
par(mfrow = c(4, 6), mar = rep(2, 4))
for (j in seq_along(listRDS)){
  
  post_k.tmp = post_k[[j]]
  prior_k.tmp = prior_k[[j]]
  
  plot(10:40, prior_k.tmp[10:40], ylim = c(0, max(c(post_k.tmp, prior_k.tmp))), type = "h", 
       xlab = "", ylab = "")
  segments(as.numeric(names(post_k.tmp)), rep(0, length(post_k.tmp)), 
           as.numeric(names(post_k.tmp)), post_k.tmp, 
           col = 2)
  
}
if (SAVE){dev.off()}

# Partition point estimates ----
ppe = readRDS("output/AP_indmod/ppe.RDS")
psm = readRDS("output/AP_indmod/psm.RDS")
# psm = replicate(nyears, matrix(NA, nrow = nIDs, ncol = nIDs), simplify = F)
# ppe = matrix(NA, nrow = nIDs, ncol = nyears)
# rownames(ppe) = IDs; colnames(ppe) = years
# for (j in 1:nyears){
#   cat(j, "\t")
#   ord = order(rownames(listRDS[[j]]$data$ghe.chr))
#   sample.tmp = listRDS[[j]]$sim_ghe$C[-burn, ]
#   psm.tmp = comp.psm(sample.tmp)
#   psm[[j]] = psm.tmp[ord, ord]
#   ppe.tmp = salso(sample.tmp)
#   ppe[ , j] = ppe.tmp[ord]
# }
# saveRDS(ppe, "output/AP_indmod/ppe.RDS")
# saveRDS(psm, "output/AP_indmod/psm.RDS")
ppe.seqsort = seqsort(ppe); dimnames(ppe.seqsort) = dimnames(ppe)
ppe.df = ppe.seqsort %>% reshape2::melt(varnames = c("ID", "Year"), value.name = "Cl.ss")
ghe_ppe = left_join(ghe, ppe.df, by = c("ID", "Year"))

## NMI ----
NMI_mat = matrix(NA, nrow = nyears, ncol = nyears)
for (i in 1:nyears){
  for (j in 1:nyears){
    NMI_mat[i, j] = aricode::NMI(ppe[, i], ppe[, j])
  }
}
dimnames(NMI_mat) = list(years, years)
plt_NMI = NMI_mat %>% reshape2::melt(varnames = c("Year1", "Year2"), value.name = "NMI") %>% 
  ggplot(aes(x = Year1, y = Year2)) +
  geom_tile(aes(fill = NMI)) + 
  scale_fill_continuous(type = "viridis", limits = c(0, 1)) +
  labs(title = "NMI independent ArgientoPaci models",
       x = "Year", y = "Year") +
  theme(text = element_text(size = 7))

if (SHOW) plt_NMI

if (SAVE) ggsave(plot = plt_NMI,
                 filename = "NMI_years.pdf",
                 path = "img/AP_results",
                 height = 4, width = 6)

# Entropy Sex and Region within each estimated cluster ----
plt_entropy = ghe_ppe %>% group_by(Year, Cl.ss) %>% 
  summarise(n = n(), Sex = Entropy(table(Sex)), Region = Entropy(table(Region))) %>% 
  pivot_longer(-c(Year, Cl.ss, n), names_to = "Variable", values_to = "Entropy") %>%
  mutate(NormEntropy = Entropy / log(n)) %>%
  ggplot(aes(x = Year, y = NormEntropy, color = as.character(Variable))) +
  scale_color_discrete("") + 
  geom_line() +
  geom_point() +
  facet_wrap(~ Cl.ss)

if (SHOW) plt_entropy

if (SAVE) ggsave(plot = plt_entropy,
                 filename = "entropy_withinclusters.pdf",
                 path = "img/AP_results",
                 height = 6, width = 6)

# marginal posterior mode of c_jk ----
# post.cent = lapply(listRDS, function(RDS) apply(RDS$gibbs_param$Cent, c(2,3), moda))
# for (j in seq_along(post.cent)){
#   colnames(post.cent[[j]]) = ages
#   rownames(post.cent[[j]]) = 1:nrow(post.cent[[j]])
# }
# 
# post.cent_lab = list()
# for (j in seq_along(post.cent)){
#   tmp = matrix(NA, nrow = nrow(post.cent[[j]]), ncol = ncol(post.cent[[j]]))
#   for (i in 1:nrow(post.cent[[j]])) {
#     for (a in 1:nages){
#       ghe.lab = listRDS[[j]]$data$ghe.lab %>% filter(Age == ages[a])
#       tmp[i,a] = as.character(ghe.lab$CauseS[ghe.lab$CauseI == post.cent[[j]][i,a]])
#     }
#   }
#   dimnames(tmp) = dimnames(post.cent[[j]])
#   post.cent_lab[[j]] = tmp
# }
# 
# 
# for (y in 1:nyears){
#   yyyy = 1999+y
#   cat(yyyy, "\t")
#   
#   tab = table(ppe[,y], ppe.seqsort[,y])
#   Cl2Clss = apply(tab, 1, function(x) which(x > 0)) %>% as.data.frame() %>% rownames_to_column("Cl") %>% rename(Cl.ss = ".")
#   
#   plt = post.cent_lab[[y]] %>% as.data.frame() %>% rownames_to_column("Cl") %>%
#     left_join(., Cl2Clss, by = "Cl") %>%
#     pivot_longer(-c(Cl, Cl.ss), names_to = "Age", values_to = "CauseS") %>%
#     mutate(Age = factor(Age, levels = levels(ghe$Age))) %>%
#     ggplot(aes(x = Age, y = Cl.ss)) +
#     geom_tile(aes(fill = CauseS), color = "white") +
#     geom_text(aes(label = CauseS), size = 2) +
#     facet_grid(rows = vars(Cl.ss), scales = "free_y") +
#     scale_fill_manual("", values = palette) +
#     # guides(fill = guide_legend(nrow = 2)) +
#     theme(legend.position = "none",
#           legend.text = element_text(size = 8),
#           panel.grid = element_blank(),
#           plot.background = element_rect(fill = "white"),
#           panel.background = element_rect(fill = "white"),
#           axis.title.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           axis.text.y = element_blank(),
#           axis.text.x = element_text(angle = 45, hjust = 1))
#   
#   ggsave(plot = plt,
#          filename = paste0("heatmap_centers", yyyy, ".pdf"),
#          path = "img/AP_results/heatmaps_centers",
#          height = 6, width = 12)
# }

# Heatmaps data ordered by clusters ----
for (yyyy in years){
  cat(yyyy, "\t")
  plt = ghe_ppe %>% 
    filter(Year == yyyy) %>%
    ggplot(aes(x = Age, y = ID, fill = CauseS)) +
    geom_raster() + 
    facet_grid(rows = vars(Cl.ss), scales = "free_y", space = "free_y") + 
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
                      select(ID, Cl.ss, Sex, Region)) %>% pivot_longer(cols = -c(ID, Cl.ss), , names_to = "x", values_to = "fill") %>% 
    mutate(x = factor(x, levels = c("Sex", "Region")))
  
  colors_aux = c("lightblue", "pink", "#EF476F", "#F78C6B", "#FFD166", "#06D6A0", "#118AB2", "#073B4C")
  names(colors_aux) = c("M", "F", "AF", "AM", "EM", "EU", "SEA", "WP")
  
  plt_aux = ggplot(data = data_aux) + 
    geom_tile(aes(x = x, y = ID, fill = fill)) +
    facet_grid(rows = vars(Cl.ss), scales = "free_y", space = "free_y") + 
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
         filename = paste0("heatmap_", yyyy, ".pdf"),
         path = "img/AP_results/heatmaps_orderppe",
         height = 40, width = 16)
}


# Coclustering matrices ----
for (j in seq_along(psm)){rownames(psm[[j]]) = colnames(psm[[j]]) = IDs}
names(psm) = years

for (y in c(1, 6, 11, 15, 20, 22)){
  
  yyyy = years[y]
  
  psm.df = psm[[y]] %>% reshape2::melt(varnames = c("ID1", "ID2"), value.name = "ProbCC") %>% 
    left_join(ppe.df %>% filter(Year == yyyy) %>% select(-Year) %>% rename(Cl1 = Cl.ss), by = c("ID1" = "ID")) %>% 
    left_join(ppe.df %>% filter(Year == yyyy) %>% select(-Year) %>% rename(Cl2 = Cl.ss), by = c("ID2" = "ID"))
  
  cat(yyyy, "\t")
  plt_out = psm.df %>% 
    ggplot(aes(x = ID1, y = ID2, fill = ProbCC)) +
    geom_raster() + 
    facet_grid(Cl2 ~ Cl1, scales = "free", space = "free") +
    scale_fill_viridis_c(begin = 0, end = 1) + 
    ggtitle(paste0("Posterior prob. of coclust. in ", year)) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
          # axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          # axis.text.y = element_text(angle = 0, hjust = 1, size = 6))
  
  ggsave(plot = plt_out,
         filename = paste0("probcc_", yyyy, ".pdf"),
         path = "img/AP_results/probcc",
         height = 7, width = 8)
}
