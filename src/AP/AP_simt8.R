options(warn = 1)
suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
  library(abind)
  library(RcppArmadillo)
})
source("src/utils.R")
source("src/tHMM_utils.R")
source("src/tHMM_updates.R")
source("src/tHMM_gibbs.R")
source('src/ArgientoPaci/code/gibbs_sampler2.R', echo=TRUE)
source('src/ArgientoPaci/code/complement_functions.R')
Rcpp::sourceCpp('src/cpp_utils.cpp')
Rcpp::sourceCpp('src/ArgientoPaci/code/gibbs_utility.cpp')
Rcpp::sourceCpp('src/ArgientoPaci/code/hyperg2.cpp')

# Seed ----
seed = commandArgs(trailingOnly = TRUE) %>% as.numeric()

# Simulate data ----
sim = simdata(seed = seed)
Y = sim$Y
p = sim$info$p
m = apply(sim$Y, 2, function(x) length(unique(c(x))))
n = sim$info$n
nT = sim$info$nT

# Fit model ----
u = c(5, rep(4, p-2), 4.5)
v = rep(0.2, p)

# number of clusters a priori
Kstar = 9
Lambda = Kstar
gam = AntMAN::AM_find_gamma_Pois(n=n, Lambda=Lambda,Kstar=Kstar)
prior = AntMAN::AM_prior_K_Pois(n=n, gam, Lambda = Lambda)


## initial values
k.init = ceiling(0.5*n)
M.na.init = n-k.init
seed_init = seed
set.seed(seed_init)  
C.init = sample(1:k.init,n,replace=T)
cent.init = matrix(sample(1:m[1],k.init+M.na.init,replace=T),nrow = (k.init+M.na.init), ncol=p)
sigma.init = matrix(0.5,nrow = (k.init+M.na.init), ncol=p)


### RUN GIBBS SAMPLER ###
burn = 5
G = 30*10^3

seed_sim = seed
set.seed(seed_sim)
sim_ghe = gibbs_mix_con(G=G,
                        burnin = burn,
                        data=Y[ , , 8],
                        u=u,v=v,
                        Lambda = Lambda,
                        gam = gam,
                        C.init = C.init,
                        k.init = k.init,
                        cent.init = cent.init,
                        sigma.init = sigma.init,
                        M.na.init =M.na.init,
                        M.max=n)

saveRDS(list(ppe = ppe, NMI = NMI, data = sim),
        paste0("output/tHMM/simstudy_APt8/summary_AP", seedrun, ".RDS"))


ppe = salso::salso(sim_ghe$C[-(1:10002), ]) %>% as.vector()
NMI = aricode::NMI(ppe, sim$C[,8])

saveRDS(list(ppe = ppe, NMI = NMI, data = sim),
        paste0("output/tHMM/simstudy_APt8/summary_AP", seedrun, ".RDS"))

