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

n = dim(ghe.array)[1]
p = dim(ghe.array)[2]
.T = dim(ghe.array)[3]
u = rep(4, p)
v = rep(0.25, p)

eta = 0.405
{
  probs_gnedin = HGnedin(n, 1:n, gamma = eta)
  cat("prior expected ncl with eta ", eta, " is ", round(sum(1:n*probs_gnedin), 2), "\n", sep = "")
}

seedrun = 20148
nburn = 15000
nchain = 15000
printstep = 500
urntype = "Gnedin"

out = tHMM_gibbs(
  Y = ghe.array,
  m = m,
  # Parameters for the likelihood-specific local parameters
  par_likelihood = list(u = u, # First hyperparam HIG
                        v = v # Second hyperparam HIG
  ),
  fix_alpha.flag = TRUE, fix_alpha.value = 0.01,
  # Parameters for the tRPM process prior
  par_tRPM = list(a_alpha = 1, b_alpha = 1, # Beta prior on alpha
                  urn_type = urntype, # Type of urn process
                  eta = eta # Gnedin parameter
  ),
  # Control parameters for MCMC settings
  ctr_mcmc = list(seed = seedrun, nburnin = nburn, nchain = nchain, print_step = printstep, verbose = "1"),
  # Control parameters for result storing
  ctr_save = list(save = TRUE, filepath = "output/tHMM/",
                  filename = paste("res_GnedinFixAlpha", seedrun, ".RDS", sep = ""))
)
