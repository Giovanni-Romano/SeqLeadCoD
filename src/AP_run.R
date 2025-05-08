source("renv/activate.R")
library(mcclust)
library(mcclust.ext)
library(AntMAN)
library(klaR)
library(TraMineR)
library(readr)
suppressPackageStartupMessages(library(tidyverse))
source('src/ArgientoPaci/code/gibbs_sampler2.R', echo=TRUE)
source('src/ArgientoPaci/code/complement_functions.R')


#GHE data ----
ghe = readRDS("data/rds/dataGHE.RDS") %>% 
  unite("ID", c(CountryN, Sex), remove = F) %>% 
  select(Year, ID, Country, CountryN, Age, CauseS, CauseT, CauseC, Region, Sex)


## Select single year ----
year = 1999+as.numeric(commandArgs(trailingOnly = TRUE))

ghe.lab = unique(ghe %>% filter(Year == year) %>% select(c("Age", "CauseS"))) %>% 
  group_by(Age) %>% 
  mutate(CauseI = as.integer(factor(CauseS)))

nIDs = length(unique(ghe$ID))
ages = levels(ghe$Age)
nages = length(ages)

ghe.chr = ghe %>% filter(Year == year) %>% select(c("ID", "Age", "CauseS")) %>% 
  arrange(Age) %>%
  pivot_wider(names_from = Age, values_from = CauseS) %>% 
  arrange(ID) %>% 
  column_to_rownames("ID") %>% as.matrix()

ghe.data = matrix(0, nrow = nIDs, ncol = nages)

for (j in 1:nages) {
  lab.tmp = ghe.lab %>% filter(Age == ages[j])
  
  for (i in 1:nIDs) {
    ghe.data[i,j] = lab.tmp$CauseI[lab.tmp$CauseS == ghe.chr[i,j]]
  }
  
}


p = nages
n = nIDs
m = ghe.lab %>% 
  group_by(Age) %>% 
  summarise(m = length(unique(CauseS))) %>% 
  pull(m)


###########################################################################

# number of clusters a priori
Kstar = 15

# prior hyperparameters
Lambda = Kstar
gam = AntMAN::AM_find_gamma_Pois(n=n, Lambda=Lambda,Kstar=15)
# gam = 1
prior = AM_prior_K_Pois(n=n, gam, Lambda = Lambda)
# u = rep(4, 12)#c(rep(4.5, 4), rep(4, 4), rep(3.5, 4))
# v = rep(0.25, length(m_unique))
# val = cbind(m_unique, u, v)
# nsig = 10^4
# S  = apply(val, 1, function(x) rhyper_sig2(n = nsig, c = x[2], d = x[3], m = x[1]), simplify = T)
# 
# G_max <- 1 - 1 / m_unique
# 
# gini_values = {
#   expS = (exp(2 / S) + (matrix(rep(m_unique - 1, each = nrow(S)), ncol = length(m_unique)))) /
#     (exp(1 / S) + (matrix(rep(m_unique - 1, each = nrow(S)), ncol = length(m_unique))))^2
#   expS[S < 0.01] = 1
#   gini_values = t( (1 - t(expS)) / G_max )
# }
# 
# par(mfrow = c(3, 4))
# for (j in 1:ncol(gini_values)){
#   hist(gini_values[ , j], freq = F, breaks = 21, main = paste0("N-Gini, v=", round(u[j], 2), ", w=", round(v[j], 2)),
#        xlab = "Gini values", ylab = "Density", col = "lightblue", border = "black")
# }
u = rep(4, p)
v = rep(0.25, p)

## initial values
k.init = ceiling(0.9*n)
M.na.init = n-k.init
seed_init = 20146
set.seed(seed_init)  
C.init = sample(1:k.init,n,replace=T)
cent.init = matrix(sample(1:m[1],k.init+M.na.init,replace=T),nrow = (k.init+M.na.init), ncol=p)
sigma.init = matrix(0.5,nrow = (k.init+M.na.init), ncol=p)


### RUN GIBBS SAMPLER ###
burn = 5
G = 35*10^3

seed_sim = 20146
set.seed(seed_sim)
sim_ghe = gibbs_mix_con(G=G,
                        burnin = burn,
                        data=ghe.data,
                        u=u,v=v,
                        Lambda = Lambda,
                        gam = gam,
                        C.init = C.init,
                        k.init = k.init,
                        cent.init = cent.init,
                        sigma.init = sigma.init,
                        M.na.init =M.na.init,
                        M.max=n)


g.idx = 2+(0:G)
# posterior K
post_k = table(sim_ghe$k[g.idx])/length(g.idx)

# posterior similarity matrix
psm = comp.psm(sim_ghe$C[g.idx, ])

# estimated clusters
pred_VI = minVI(psm)$cl
table(pred_VI)


foldersave = "output/AP_prova/"
if (!dir.exists(foldersave)) {
  dir.create(foldersave, recursive = TRUE)
}
saveRDS(list(
  sim_ghe = sim_ghe,
  param = list("v" = u, "w" = v, "Lambda" = Lambda, "gam" = gam, "Kstar" = Kstar),
  init_values = list("K" = k.init, "M.na" = M.na.init, "C" = C.init,
                     "cent" = cent.init, "sigma" = sigma.init),
  seed = list("init" = seed_init, "sim" = seed_sim),
  inference = list("post_k" = post_k,
                   "psm" = psm,
                   "pred_VI" = pred_VI),
  data = list("year" = year, "data" = ghe.data, "ghe.lab" = ghe.lab, "ghe.chr" = ghe.chr)
), file = paste0(foldersave, "/APout_", year, ".rds"))

###### parameters estimation #########
source("src/ArgientoPaci/code/gibbs_param_estim2.R")
Kest = length(unique(pred_VI))
cent.last = sigma.last = matrix(NA, Kest, p)
for(i in 1:Kest){
  cent.last[i,] = sim_ghe$Cent[[i]][nrow(sim_ghe$Cent[[i]]),]
  sigma.last[i,] = sim_ghe$Sigma[[i]][nrow(sim_ghe$Cent[[i]]),]
}

seed_param = 5
set.seed(seed_param)
gibbs_param = gibbs_param_estim(G=5*10^3,
                                gam=gam,
                                data=ghe.data,
                                u=u,
                                v=v,
                                estim_Truth = pred_VI,
                                cent.init = cent.last,
                                sigma.init = sigma.last,
                                Sm.init = sim_ghe$Sm[[length(sim_ghe$Sm)]][1:Kest]
)

saveRDS(list(
  sim_ghe = sim_ghe,
  param = list("v" = u, "w" = v, "Lambda" = Lambda, "gam" = gam, "Kstar" = Kstar),
  init_values = list("K" = k.init, "M.na" = M.na.init, "C" = C.init,
                     "cent" = cent.init, "sigma" = sigma.init),
  seed = list("init" = seed_init, "sim" = seed_sim),
  inference = list("post_k" = post_k,
                   "psm" = psm,
                   "pred_VI" = pred_VI,
                   "Kest" = Kest),
  gibbs_param = gibbs_param,
  data = list("year" = year, "data" = ghe.data, "ghe.lab" = ghe.lab, "ghe.chr" = ghe.chr)
), file = paste0(foldersave, "/APout_", year, ".rds"))
