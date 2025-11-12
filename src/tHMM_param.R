options(warn = 1)
suppressPackageStartupMessages(library(tidyverse))
library(abind)
library(RcppArmadillo)
source("src/utils.R")
source("src/tHMM_utils.R")
source("src/tHMM_updates.R")
source("src/tHMM_gibbs.R")
Rcpp::sourceCpp('src/cpp_utils.cpp')
Rcpp::sourceCpp('src/ArgientoPaci/code/gibbs_utility.cpp')
Rcpp::sourceCpp('src/ArgientoPaci/code/hyperg2.cpp')


SEX = "female"
SEED = 20019
NAME = paste0("Gnedin", SEED)
OUTFOLDER = paste0("output/tHMM/5yrs/PT/", SEX, "/")
res = readRDS(paste0(OUTFOLDER, "res_", NAME, ".RDS"))
PPE = readRDS(paste0(OUTFOLDER, "PPE_all.RDS"))[[as.character(SEED)]]

Y = res$Y
m = res$m
par_lik = res$par$par_likelihood


seedrun = SEED
urntype = "Gnedin"
nburn = 5000
nchain = 15000
printstep = 1000

out = tHMM_gibbs_parameters(
  PPE = PPE, Y = Y, m = m,
  par_likelihood = par_lik,
  ctr_mcmc = list(seed = seedrun, nburnin = nburn, nchain = nchain, print_step = printstep, 
                  verbose = "1"),
  ctr_save = list(save = TRUE, filepath = OUTFOLDER,
                  filename = paste("param_", urntype, seedrun, ".RDS", sep = ""))
)


mu.mode = lapply(out$mu, function(x) apply(x, c(1, 2), moda))
str(mu.mode)


ghe = readRDS("data/rds/dataGHE.RDS") %>%
  unite("ID", c(CountryN, Sex), remove = FALSE)
ghe.lab = unique(ghe %>% select(c("Age", "CauseS"))) %>%
  group_by(Age) %>%
  mutate(CauseI = as.integer(factor(CauseS)))

mu.mode_df = mu.mode %>% reshape2::melt() %>%
  rename(Cluster = Var1, Age = Var2, CauseI = value, Year = L1) %>%
  left_join(., ghe.lab %>% mutate(Age = as.integer(Age)), by = c("Age", "CauseI"))


saveRDS(mu.mode_df, paste0(OUTFOLDER, "muestdf_", SEED, ".RDS"))
