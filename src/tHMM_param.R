options(warn = 1)
suppressPackageStartupMessages(library(tidyverse))
library(abind)
library(RcppArmadillo)
source("src/tHMM_utils.R")
source("src/tHMM_updates.R")
source("src/tHMM_gibbs.R")
Rcpp::sourceCpp('src/cpp_utils.cpp')
Rcpp::sourceCpp('src/ArgientoPaci/code/gibbs_utility.cpp')
Rcpp::sourceCpp('src/ArgientoPaci/code/hyperg2.cpp')

PPE = readRDS("output/tHMM/PPE_Gnedin20019.RDS")
res = readRDS("output/tHMM/res_Gnedin20019.RDS")
Y = res$input$Y
m = res$input$m
par_lik = res$input$par$par_likelihood


seedrun = 20019
urntype = "Gnedin"
nburn = 5000
nchain = 15000
printstep = 1000

out = tHMM_gibbs_parameters(
  PPE = PPE, Y = Y, m = m,
  par_likelihood = par_lik,
  ctr_mcmc = list(seed = seedrun, nburnin = nburn, nchain = nchain, print_step = printstep, 
                  verbose = "1"),
  ctr_save = list(save = FALSE, filepath = "output/tHMM/",
                  filename = paste("param_", urntype, seedrun, ".RDS", sep = ""))
)


# mu.mode = lapply(out$mu, function(x) apply(x, c(1, 2), moda))
# str(mu.mode)
# 
# 
# ghe = readRDS("data/rds/dataGHE.RDS") %>%
#   unite("ID", c(CountryN, Sex), remove = FALSE)
# ghe.lab = unique(ghe %>% select(c("Age", "CauseS"))) %>% 
#   group_by(Age) %>% 
#   mutate(CauseI = as.integer(factor(CauseS)))
# 
# mu.mode_df = mu.mode %>% reshape2::melt() %>% 
#   rename(Cluster = Var1, Age = Var2, CauseI = value, Year = L1) %>%
#   left_join(., ghe.lab %>% mutate(Age = as.integer(Age)), by = c("Age", "CauseI"))
# 
# 
# plt = mu.mode_df %>% 
#   filter(Year == 1) %>%
#   ggplot(aes(x = Age, y = Cluster, fill = CauseS)) +
#   geom_raster() + 
#   facet_wrap(~Year) + 
#   geom_text(aes(label = CauseS), size = 2) +
#   scale_fill_manual(values = palette) +
#   theme(legend.position = "none",
#         panel.grid = element_blank(),
#         plot.background = element_rect(fill = "white"),
#         panel.background = element_rect(fill = "white"),
#         axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1))
# 
# sigma.mean = lapply(out$sigma, function(x) apply(x, c(1, 2), mean))
# 
