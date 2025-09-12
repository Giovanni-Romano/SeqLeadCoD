options(warn = 1)
suppressPackageStartupMessages(library(tidyverse))
library(abind)
library(RcppArmadillo)
library(ggh4x)
source("src/tHMM_utils.R")
source("src/tHMM_updates.R")
source("src/tHMM_gibbs.R")
# Rcpp::sourceCpp('src/cpp_utils.cpp')
# Rcpp::sourceCpp('src/ArgientoPaci/code/gibbs_utility.cpp')
# Rcpp::sourceCpp('src/ArgientoPaci/code/hyperg2.cpp')

ghe = readRDS("data/rds/dataGHE.RDS") %>%
  unite("ID", c(CountryN, Sex), remove = FALSE)
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

# Load ghe_PPE
cmndargs = c("Gnedin", "20019")#commandArgs(trailingOnly = TRUE)
TYPE = cmndargs[1]
SEED = cmndargs[2]
NAME = paste0(TYPE, SEED)
OUTFOLDER = "output/tHMM/"
ghe_PPE = readRDS(paste0(OUTFOLDER, "ghe_ppe_", NAME, ".RDS")) %>% 
  rename(Region_WHO = Region) %>% 
  left_join(alternative_regions, by = "Country") %>% 
  mutate(Region.short = factor(case_when(
    Region == "Western Offshoots" ~ "WO",
    Region == "Latin America & Caribbean" ~ "LAC",
    Region == "Sub-Saharan Africa" ~ "SSA",
    Region == "Middle East & North Africa" ~ "MENA",
    Region == "Europe (Non P-S)" ~ "Eu",
    Region == "Post-Sovietic" ~ "P-S",
    Region == "South Asia" ~ "SA",
    Region == "East Asia & Pacific" ~ "EAP",
    TRUE ~ Region
  ), levels = c("SSA", "MENA", "LAC", "WO", "Eu", "P-S", "SA", "EAP"), ordered = TRUE))

# Pairwise hamming distances between rows ----
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

# Selected years for plots ----
yearssub = c(2000, 2007, 2014, 2021)


HSM.df = HSM %>% reshape2::melt() %>%
  mutate(value = value/19)

ghe_PPE2000 = ghe_PPE %>% filter(Year == 2000) %>% select(ID, order) %>% unique()
countryorder = ghe_PPE2000 %>% arrange(order) %>% pull(ID)

HSM.dfplot = HSM.df %>% 
  filter(Year %in% yearssub) %>% 
  left_join(ghe_PPE %>% select(Year, ID, Region.short, order) %>% unique() %>% rename(Region1 = Region.short, Order1 = order), by = c("Year", "Country1" = "ID")) %>%
  left_join(ghe_PPE %>% select(Year, ID, Region.short, order) %>% unique() %>% rename(Region2 = Region.short, Order2 = order), by = c("Year", "Country2" = "ID")) %>%
  mutate(Country1 = factor(Country1, levels = countryorder, ordered = T),
         Country2 = factor(Country2, levels = countryorder, ordered = T))

listplot_HSM = lapply(yearssub, function(y)
  ggplot(HSM.dfplot %>% filter(Year == y)) +
    geom_raster(aes(x = Country1, y = Country2, fill = value)) +
    scale_fill_gradient(low = "#edf6f9", high = "#264653", limits = c(0, 1)) +
    facet_nested(rows = vars(Region2), cols = vars(Year, Region1), scales = "free", space = "free",
                 strip = strip_nested()) +
    # coord_equal() + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          panel.spacing = unit(0.1, "lines"))
)



{
  width = 15; height = 4; zoom = 1
  ggsave("img/tHMM/Gnedin20019/HSM.pdf",
         gridExtra::grid.arrange(grobs = listplot_HSM, 
                                 nrow = 1),
         width = width * zoom, height = height * zoom)
}

# Selected years for plots for slides ----
yearssub = c(2000, 2010, 2020)


HSM.df = HSM %>% reshape2::melt() %>%
  mutate(value = value/19)

ghe_PPE2000 = ghe_PPE %>% filter(Year == 2000) %>% select(ID, order) %>% unique()
countryorder = ghe_PPE2000 %>% arrange(order) %>% pull(ID)

HSM.dfplot = HSM.df %>% 
  filter(Year %in% yearssub) %>% 
  left_join(ghe_PPE %>% select(Year, ID, Region.short, order) %>% unique() %>% rename(Region1 = Region.short, Order1 = order), by = c("Year", "Country1" = "ID")) %>%
  left_join(ghe_PPE %>% select(Year, ID, Region.short, order) %>% unique() %>% rename(Region2 = Region.short, Order2 = order), by = c("Year", "Country2" = "ID")) %>%
  mutate(Country1 = factor(Country1, levels = countryorder, ordered = T),
         Country2 = factor(Country2, levels = countryorder, ordered = T))

listplot_HSM = lapply(yearssub, function(y)
  ggplot(HSM.dfplot %>% filter(Year == y)) +
    geom_raster(aes(x = Country1, y = Country2, fill = value)) +
    scale_fill_gradient(low = "#edf6f9", high = "#264653", limits = c(0, 1)) +
    facet_nested(rows = vars(Region2), cols = vars(Year, Region1), scales = "free", space = "free",
                 strip = strip_nested()) +
    # coord_equal() + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          panel.spacing = unit(0.1, "lines"))
)



{
  width = 12; height = 4.2; zoom = 1
  ggsave("img/tHMM/Gnedin20019/HSM_slides.pdf",
         gridExtra::grid.arrange(grobs = listplot_HSM, 
                                 nrow = 1),
         width = width * zoom, height = height * zoom)
}




# # HDM plot ----
# library(seriation)
# library(ggh4x)
# regions <- unique(ghe$Region) %>% as.character()
# method = "complete"
# {
#   # Create empty lists to collect row and column orderings
#   row_ordered <- c()
#   col_ordered <- c()
#   
#   # Create an environment to avoid duplicates
#   row_seen <- new.env()
#   col_seen <- new.env()
#   
#   
#   get_seriated_order <- function(mat, method = "ward.D2") {
#     # mat[] = as.numeric(mat)
#     # o = seriate(mat, margin = 1, ...)
#     # 
#     # get_order(o)
#     p = pheatmap::pheatmap(mat, clustering_method = method, silent = TRUE)
#     p$tree_row$order
#   }
#   
#   for (r1 in regions) {
#     # for (r2 in regions) {
#     
#     rows <- country_info %>% filter(Region == r1) %>% pull(ID)
#     cols <- country_info %>% filter(Region == r1) %>% pull(ID)
#     
#     block <- HSM[rows, cols, 1]
#     
#     # Seriate rows
#     row_ord <- rows[get_seriated_order(as.matrix(block), method)]
#     col_ord <- cols[get_seriated_order(as.matrix(t(block)), method)]
#     
#     # Add to master lists in order — avoiding duplicates
#     for (r in row_ord) {
#       if (!exists(r, envir = row_seen)) {
#         row_ordered <- c(row_ordered, r)
#         assign(r, TRUE, envir = row_seen)
#       }
#     }
#     for (c in col_ord) {
#       if (!exists(c, envir = col_seen)) {
#         col_ordered <- c(col_ordered, c)
#         assign(c, TRUE, envir = col_seen)
#       }
#     }
#     # }
#   }
#   
#   all_countries <- country_info$ID
#   row_ordered <- unique(c(row_ordered, all_countries))
#   col_ordered <- unique(c(col_ordered, all_countries))
#   
#   HSM.dfplot = HSM %>% reshape2::melt() %>%
#     mutate(value = value/19) %>% 
#     # mutate(value = factor(value, levels = 0:19, ordered = T)) %>%
#     left_join(ghe %>% select(ID, Region) %>% unique() %>% rename(Region1 = Region), by = c("Country1" = "ID")) %>%
#     left_join(ghe %>% select(ID, Region) %>% unique() %>% rename(Region2 = Region), by = c("Country2" = "ID")) %>%
#     mutate(Country1 = fct_rev(factor(Country1, levels = col_ordered, ordered = T)),
#            Country2 = factor(Country2, levels = row_ordered, ordered = T)) %>% 
#     filter(Year %in% yearssub)
#   
#   listplot_HSM = lapply(yearssub, function(y)
#     HSM.dfplot %>%
#       filter(Year == y) %>% 
#       ggplot(aes(x = Country1, y = Country2, fill = value)) +
#       geom_raster() +
#       # scale_fill_viridis_c(name = "Hamming Distance") +
#       # scale_fill_distiller(palette = "Greys", direction = 1, 
#       #                       limits = c(0, 1)) +
#       scale_fill_gradient(low = "#edf6f9", high = "#264653", limits = c(0, 1)) +
#       facet_nested(rows = vars(Region2), cols = vars(Year, Region1), scales = "free", space = "free",
#                    strip = strip_nested()) +
#       theme(axis.text = element_blank(),
#             axis.ticks = element_blank(),
#             axis.title = element_blank(),
#             legend.position = "none",
#             panel.spacing = unit(0.1, "lines"))
#   )
#   
#   # gridExtra::grid.arrange(grobs = listplot_HDM, ncol = 4)
#   
#   {
#     width = 15; height = 4; zoom = 0.75
#     ggsave(paste0("img/tHMM/Gnedin20019/HSM_", method, ".pdf"),
#            gridExtra::grid.arrange(grobs = listplot_HSM, 
#                                    nrow = 1),
#            width = width * zoom, height = height * zoom)
#   }
# }
# 
# ## MDS to find positions of nodes in networks ----
# set.seed(1234)
# mds_results = apply(HDM+1e-3, 3, function(x) {
#   MASS::isoMDS(x, k = 2)
# }, simplify = F)
# 
# 
# ## Plot MDS ----
# colors_aux = c("#EF476F", "#F78C6B", "#FFD166", "#06D6A0", "#118AB2", "#073B4C")
# names(colors_aux) = c("AF", "AM", "EM", "EU", "SEA", "WP")
# 
# plot_MDS = lapply(mds_results, function(x) x$points %>% as.data.frame %>% rownames_to_column("ID")) %>%
#   bind_rows(.id = "Year") %>%
#   left_join(ghe %>% select(ID, Region) %>% unique(), by = c("ID" = "ID")) %>% 
#   filter(Year %in% yearssub) %>% 
#   ggplot() +
#   geom_point(aes(x = V1, y = V2, color = Region), size = 1.5) +
#   # geom_text(aes(x = V1, y = V2, label = ID), size = 3, hjust = 0.5, vjust = -0.5) +
#   scale_color_manual(values = colors_aux) +
#   facet_wrap(~ Year, nrow = 1) +
#   labs(x = "MDS Dimension 1", y = "MDS Dimension 2")
# 
# {
#   width = 15; height = 3.5; zoom = 1
#   ggsave("img/tHMM/Gnedin20019/MDS.pdf", plot_MDS, 
#          width = width * zoom, height = height * zoom)
#   }
# 
# ## Plot MDS all years ----
# plot_MDS.all = lapply(mds_results, function(x) x$points %>% as.data.frame %>% rownames_to_column("ID")) %>%
#   bind_rows(.id = "Year") %>%
#   left_join(ghe %>% select(ID, Region) %>% unique(), by = c("ID" = "ID")) %>% 
#   ggplot() +
#   geom_point(aes(x = V1, y = V2, color = Region), size = 1.5) +
#   # geom_text(aes(x = V1, y = V2, label = ID), size = 3, hjust = 0.5, vjust = -0.5) +
#   scale_color_manual(values = colors_aux) +
#   facet_wrap(~ Year, nrow = 5) +
#   labs(x = "MDS Dimension 1", y = "MDS Dimension 2")
# 
# {
#   width = 15; height = 15; zoom = 1
#   ggsave("img/tHMM/Gnedin20019/MDS_all.pdf", plot_MDS.all, 
#          width = width * zoom, height = height * zoom)
#   }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Network plots ----
# 
# # binary_adjacency = {
# #   tmp = HSM
# #   tmp[tmp < 10] = 0
# #   tmp[tmp >= 10] = 1
# #   tmp
# # }
# 
# thres_adjacency = {
#   thres = 18
#   tmp = HSM
#   tmp[tmp <= thres] = 0
#   tmp = (tmp - thres) / (19 - thres)
#   tmp[tmp < 0] = 0
#   tmp
# }
# 
# library(igraph)
# library(ggraph)
# library(graphlayouts)
# library(ggforce)
# library(tidygraph)
# 
# listplot_net = list()
# for (j in seq_along(yearssub)) {
#   t = yearssub[j] - 2000 + 1
#   cat(j/length(yearssub), "\t")
#   g = graph_from_adjacency_matrix(thres_adjacency[,,t] %>% unname, mode = "undirected", weighted = TRUE, diag = FALSE)
#   comm = cluster_louvain(g, weights = E(g)$weight)
#   
#   V(g)$ID = rownames(PPE)
#   V(g)$community = membership(comm) %>% as.character()
#   V(g)$size = degree(g)
#   V(g)$Region = (ghe %>% select(ID, Region) %>% unique() %>% pull(Region))[match(V(g)$ID, ghe %>% select(ID, Region) %>% unique() %>% pull(ID))]
#   
#   ref_year = (2000 - 2000) + 1
#   pos = mds_results[[ref_year]]$points
#   listplot_net[[j]] = ggraph(g, layout = "manual", x = pos[, 1], y = pos[, 2]) +
#     geom_edge_link0(aes(color = weight, alpha = weight), show.legend = FALSE) +
#     geom_node_point(aes(fill = Region), shape = 21, size = 5) +
#     scale_edge_color_gradient(low = "grey90", high = "grey20") +
#     theme(legend.position = "none")
# }
# 
# {
#   width = 8; height = 8; zoom = 1
#   for (j in seq_along(yearssub)) {
#     ggsave(paste0("img/tHMM/Gnedin20019/net", j, ".pdf"),
#            listplot_net[[j]],
#            width = width * zoom, height = height * zoom)
#   }
# }
# 
# # graph_tbl = as_tbl_graph(g)
# # 
# # graph_tbl = graph_tbl %>%
# #   mutate(community = as.factor(V(g)$community)) %>% 
# #   left_join(unique(ghe %>% select(ID, Region)), by = "ID")
# # 
# # ggraph(graph_tbl %>% filter(size > 0), layout = "fr") +
# #   geom_edge_link(aes(alpha = weight), edge_colour = "grey66") +
# #   geom_node_point(aes(color = Region), size = 3) +
# #   theme_void()
# 
# 
# # V(net)$grp = as.character(PPE[ , t])
# # bb <- layout_as_backbone(net, keep = 0.4)
# # E(net)$col <- FALSE
# # E(net)$col[bb$backbone] <- TRUE
# # 
# # 
# # net %>% 
# #   as_tbl_graph() %>% 
# #   mutate(community = as.factor(group_infomap())) %>% 
# #   # ggraph(layout = 'kk') + 
# #   geom_edge_link(aes(alpha = after_stat(index)), show.legend = FALSE) + 
# #   geom_node_point(aes(colour = community), size = 1) + 
# #   theme_graph()
