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

cmline = c("female", 12345) #commandArgs(trailingOnly = TRUE)
cat("Argument from command line: ", cmline, sep = "")
sex = cmline[1]
seed = cmline[2] %>% as.integer

if (sex == "female"){
  ghe = readRDS("data/rds/GHEdf_female_5yrs.RDS")
} else if (sex == "male"){
  ghe = readRDS("data/rds/GHEdf_male_5yrs.RDS")
} else {
  stop("Wrong sex value")
}

ghe.lab = unique(ghe %>% select(c("Age", "CauseS"))) %>% 
  group_by(Age) %>% 
  mutate(CauseI = as.integer(factor(CauseS)))
m = ghe.lab %>% group_by(Age) %>% summarise(m = n()) %>% pull(m)

ghe.w = ghe %>% 
  left_join(ghe.lab, by = c("Age", "CauseS")) %>% 
  select(c("Country", "Age", "Year", "CauseI"))  %>% 
  arrange(Country, Age, Year)   %>% 
  pivot_wider(names_from = Age, values_from = CauseI)
ghe.array = abind(map(split(ghe.w, ghe.w$Year), 
                      \(x) x %>% select(-Year) %>% column_to_rownames("Country")), along = 3)

n = dim(ghe.array)[1]
p = dim(ghe.array)[2]
.T = dim(ghe.array)[3]
u = rep(4, p)
v = rep(0.25, p)

eta_values = seq(0.45, 0.55, by = 0.01)
exp_ncl = sapply(eta_values, function(e) {
  probs_gnedin = HGnedin(n, 1:n, gamma = e)
  round(sum(1:n*probs_gnedin))
})
cbind(eta = eta_values, exp_ncl = exp_ncl)
eta = 0.55
{
  probs_gnedin = HGnedin(n, 1:n, gamma = eta)
  cat("prior expected ncl with eta ", eta, " is ", round(sum(1:n*probs_gnedin), 2), "\n", sep = "")
}

seedrun = seed
nburn = 2
nchain = 5
printstep = 1
nclinit = 50
urntype = "Gnedin"
SAVE = TRUE

# PT parameters
.nrep = 5
.ntry = floor(.nrep / 2)
.start_swap = 1
.adaptive = TRUE
.deterministic = TRUE
.prob_target = 0.65

out = tHMM_gibbs_adaptive_PT(
  Y = ghe.array,
  m = m,
  # Parameters for the likelihood-specific local parameters
  par_likelihood = list(u = u, # First hyperparam HIG
                        v = v # Second hyperparam HIG
  ),
  # Parameters for the tRPM process prior
  par_tRPM = list(a_alpha = 1, b_alpha = 1, # Beta prior on alpha
                  urn_type = urntype, # Type of urn process
                  eta = eta # Gnedin parameter
  ),
  # Control parameters for MCMC settings
  ctr_mcmc = list(seed = seedrun, nburnin = nburn, nchain = nchain, print_step = printstep, 
                  ncl.init = nclinit, verbose = "1",
                  parallel = FALSE, nthreads = .nrep),
  # Control parameters for Parallel Tempering
  ctr_swap = list(nrep = .nrep, ntry = .ntry, start = .start_swap, deterministic = .deterministic, adaptive = .adaptive,
                  prob_target = .prob_target),
  # Control parameters for result storing
  ctr_save = list(save = SAVE, filepath = paste0("output/tHMM/5yrs/", sex, "/"),
                  filename = paste("res_", urntype, seedrun, ".RDS", sep = "")),
)
