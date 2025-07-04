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

seed = commandArgs(trailingOnly = TRUE) %>% as.numeric()
res = readRDS(paste0("output/tHMM/simstudy_simplified/summary_Gnedin", seed, ".RDS"))

PPE = res$ppe
Y = res$data$Y.converted
m = apply(Y, 2, function(x) length(unique(c(x))))
n = res$data$info$n
p = res$data$info$p
nt = res$data$info$nT
u = c(5, rep(4, p-2), 4.5)
v = rep(0.2, p)
par_lik = list(u = u, v = v)


oglab = res$data$info$unique$original_lab
convlab = res$data$info$unique$converted_lab

seedrun = seed
urntype = "Gnedin"
nburn = 5000
nchain = 15000
printstep = 5000

out = tHMM_gibbs_parameters(
  PPE = PPE, Y = Y, m = m,
  par_likelihood = par_lik,
  ctr_mcmc = list(seed = seedrun, nburnin = nburn, nchain = nchain, print_step = printstep, 
                  verbose = "1"),
  ctr_save = list(save = FALSE)
)

mu.mode = lapply(out$mu, function(x) apply(x, c(1, 2), moda))
mu.mode_oglab = list()

for (t in 1:nT){
  tmp = array(NA, dim = dim(mu.mode[[t]]))
  for (j in 1:p){
    tmp[,j] = oglab[[j]][mu.mode[[t]][, j]]
  }
  mu.mode_oglab[[t]] = tmp
}

prec_est = rep(NA, nT)
for (t in 1:8){
  M = mu.mode_oglab[[t]]
  prec_est[t] = mean(res$data$Y[,,t] ==  M[PPE[,t], ])
}

saveRDS(list(NMI = res$NMI, prec_est = prec_est),
        file = paste0("output/tHMM/simstudy_simplified/metrics_", urntype, seedrun, ".RDS"))