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


# {
#   val = cbind(m, u, v)
#   nsig = 10^4
#   S  = apply(val, 1, function(x) rhyper_sig2(n = nsig, c = x[2], d = x[3], m = x[1]), simplify = T)
# 
#   Gmax <- 1 - 1 / m
# 
#   gini_values = {
#     expS = (exp(2 / S) + (matrix(rep(m - 1, each = nrow(S)), ncol = length(m)))) /
#       (exp(1 / S) + (matrix(rep(m - 1, each = nrow(S)), ncol = length(m))))^2
#     expS[S < 0.01] = 1
#     gini_values = t( (1 - t(expS)) / Gmax )
#   }
# 
#   par(mfrow = c(2, 3))
#   for (j in 1:ncol(gini_values)){
#     hist(gini_values[ , j], freq = F, breaks = 21, main = paste0("N-Gini, v=", u[j], ", w=", v[j], ", m=", m[j]),
#          xlab = "Gini values", ylab = "Density", col = "lightblue", border = "black")
#   }
# }

eta_values = seq(0.45, 0.55, by = 0.025)
exp_ncl = sapply(eta_values, function(e) {
  probs_gnedin = HGnedin(n, 1:n, gamma = e)
  round(sum(1:n*probs_gnedin), 2)
})
cbind(eta = eta_values, exp_ncl = exp_ncl)
eta = 0.5
{
  probs_gnedin = HGnedin(n, 1:n, gamma = eta)
  cat("prior expected ncl with eta ", eta, " is ", round(sum(1:n*probs_gnedin), 2), "\n", sep = "")
}

seedrun = seed
nburn = 15000
nchain = 15000
printstep = 1000
nclinit = 50
urntype = "Gnedin"

out = tHMM_gibbs(
  Y = Y,
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
                  ncl.init = nclinit, verbose = "1"),
  # Control parameters for result storing
  ctr_save = list(save = TRUE, filepath = "output/tHMM/simstudy/",
                  filename = paste("res_", urntype, seedrun, ".RDS", sep = "")),
  ctr_alpha = list(fix_alpha.flag = FALSE)
)

ncltrace = apply(out$output$C, c(2, 3), function(x) length(unique(x)))
ncltrue = data.frame(Year = 1:nT, ncl = apply(sim$C, 2, function(x) length(unique(x))))
ncltrace %>% 
  melt(varnames = c("Year", "Iteration"), value.name = "ncl") %>% 
  filter(Iteration > 5) %>% 
  ggplot(aes(x = Iteration, y = ncl, group = Year)) +
  geom_line(aes(color = as.factor(Year)), alpha = 0.5) + 
  facet_wrap(~ Year, ncol = 2) +
  geom_hline(data = ncltrue, aes(yintercept = ncl), linetype = "dashed", color = "red")

ppe = apply(out$output$C[ , , -c(1:nburn)], 2, function(x) salso::salso(t(x)))
NMI = sapply(1:nT, function(t) aricode::NMI(sim$C[, t], ppe[, t]))

saveRDS(list(ppe = ppe, NMI = NMI, data = sim),
        paste0("output/tHMM/simstudy/summary_", urntype, seedrun, ".RDS"))

# nT = sim$info$nT
# lapply(1:nT, function(t) table(true = sim$C[,t], ppe = ppe[, t]))
# 
# ppesort = ppe
# ppesort[, 1] = refsort(ppe[, 1], sim$C[, 1])
# ppesort = seqsort(ppesort, dimsort = FALSE)
# lapply(1:nT, function(t) table(true = sim$C[,t], ppe = ppesort[, t]))
# 
# ppeflow = ppesort %>% melt(varnames = c("Subject", "Year"), value.name = "Cluster") %>%
#   mutate(freq = 1) %>% 
#   mutate(Cluster = factor(Cluster, levels = 1:max(ppe, na.rm = T)))
# 
# Cflow = sim$C %>% melt(varnames = c("Subject", "Year"), value.name = "Cluster") %>%
#   mutate(freq = 1) %>% 
#   mutate(Cluster = factor(Cluster, levels = 1:max(sim$C, na.rm = T)))
# 
# library(ggalluvial)
# library(patchwork)
# (ppeflow %>% ggplot(aes(x = Year, stratum = Cluster, alluvium = Subject,
#                                 y = freq,
#                                 fill = Cluster)) +
#   geom_flow() +
#   geom_stratum(color = "white") + 
#   geom_text(stat = "stratum", aes(label = after_stat(stratum)),
#             size = 3) + 
#   coord_cartesian(expand = TRUE, xlim = c(0.95, 8.05)) +
#   scale_x_continuous(breaks = 1:nT) +
#   theme(axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_blank())) / 
#   (Cflow %>% ggplot(aes(x = Year, stratum = Cluster, alluvium = Subject,
#                                 y = freq,
#                                 fill = Cluster)) +
#   geom_flow() +
#   geom_stratum(color = "white") +
#   geom_text(stat = "stratum", aes(label = after_stat(stratum)),
#             size = 3) +
#   coord_cartesian(expand = TRUE, xlim = c(0.95, 8.05)) +
#   scale_x_continuous(breaks = 1:nT) +
#   theme(axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_blank()))
#     