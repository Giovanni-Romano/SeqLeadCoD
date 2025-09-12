suppressPackageStartupMessages(library(tidyverse))
rdirichlet = function(n, alpha){
  Gam <- matrix(0, n, length(alpha))
  for (i in 1:length(alpha)) Gam[, i] <- rgamma(n, shape = alpha[i])
  Gam/rowSums(Gam)
}

simdata = function(n, p, nT,
                   ncl = NULL, Pi = NULL, 
                   m = NULL, m.important = NULL, values = NULL,
                   seed = NULL){
  
  # Check dimensions
  if (!(length(ncl) == nT)) {
    stop("ncl must have length equal to nT")
  }
  
  if (!(length(Pi) == nT)) {
    stop("Pi must have length equal to nT")
  }
  
  if (!(length(m) == p)) {
    stop("m must have length equal to p")
  }
  
  if (!(length(m.important) == p)) {
    stop("m.important must have length equal to p")
  }
  
  if (!all(sapply(values, length) == m)){
    stop("values must be a list of vectors each with lengths equal to corresponding element in m")
  }
  
  # Check if ncl is a vector of integers
  if (!is.numeric(ncl) || any(ncl <= 0) || any(ncl > n)) {
    stop("ncl must be a vector of positive integers smaller than n")
  }
  
  # Check if first element of Pi is a ncl[1]-vector of probabilities
  if (!is.numeric(Pi[[1]]) || length(Pi[[1]]) != ncl[1] || any(Pi[[1]] < 0) || sum(Pi[[1]]) != 1) {
    stop("Pi[[1]] must be a vector of probabilities of length ncl[1]")
  }
  
  # Check if each Pi[[t]] is a ncl[t-1] x ncl[t] matrix of transition probabilities
  for (t in 2:nT) {
    if (!is.numeric(Pi[[t]]) || 
        !is.matrix(Pi[[t]]) || 
        nrow(Pi[[t]]) != ncl[t-1] || 
        ncol(Pi[[t]]) != ncl[t] || 
        any(Pi[[t]] < 0) || 
        any(round(rowSums(Pi[[t]]), 6) != 1)) {
      stop(paste("Pi[[", t, "]] must be a matrix of transition probabilities with dimensions ncl[t] x ncl[t-1]", sep = ""))
    }
  }
  
  # Set seed
  if (is.null(seed)){
    seed = rpois(1, 1000000)
  }
  set.seed(seed)
  
  C = matrix(NA_integer_, nrow = n, ncol = nT)
  dimnames(C) = list(Subject = 1:n, Time = 1:nT)
  
  # Sample cluster allocation for time 1
  C[,1] = sample(x = 1:ncl[1], size = n, replace = TRUE, prob = Pi[[1]])
  
  # Sample cluster allocation for subsequent time points
  for (t in 2:nT) {
    for (i in 1:n){
      C[i,t] = sample(x = 1:ncl[t], size = 1, replace = TRUE, prob = Pi[[t]][C[i,t-1], ])
    }
  }
  
  # Generate multinom params
  Prob = list()
  for (j in 1:p){
    Prob[[j]] = matrix(NA, nrow = max(ncl), ncol = m[j])
    counter = 1
    for (cl in 1:max(ncl)){
      top = counter
      tmp = rep(0.1/(m[j]-1), m[j])
      tmp[top] = 0.9
      Prob[[j]][cl, ] = tmp
      counter = ifelse(counter + 1 > m.important[j], 1, counter + 1)
    }
  }
  
  # Generate data
  Y = array(NA, dim = c(n, p, nT), 
            dimnames = list(Subject = 1:n, 
                            Variable = 1:p, 
                            Time = 1:nT))
  for (t in 1:nT) {
    cl.obs = sort(unique(C[, t]))
    
    for (cl in cl.obs) {
      idx = which(C[, t] == cl)
      for (j in 1:p) {
        Y[idx, j, t] = sample(x = values[[j]], size = length(idx), replace = TRUE, prob = Prob[[j]][cl, ])
      }
    }
    
  }
  
  out = list(Y = Y, Prob = Prob, C = C, seed = seed)
}


n = 366
p = 19
nT = 20
m = c(10, 25, 50, 50, 45, 45, 
      40, 40, 40, 30, 30, 25, 
      25, 25, 25, 20, 15, 15, 15)
m.important = c(4, 6, 10, 10, rep(8, 12), rep(6, 3))
values = list(1:10, 2:26, 2:51, 2:51, 16:60, 21:65,
              26:65, 31:70, 31:70, 31:60, 31:60, 36:60,
              36:60, 41:65, 41:65, 51:70, 51:65, 51:65, 51:65)

ncl = c(15, 15, 13, 16, 14,
        12, 10, 11, 10, 12, 
        10, 8, 8, 9, 10, 
        12, 14, 13, 15, 15)
phi[1] = 0.8
phi[-1] = c(rep(0.1, 3), rep(0.2, 2), rep(0.1, 3),
            rep(0.1, 3), rep(0.2, 3), rep(0.15, 3),
            0.5, 0.4)+0.05

Pi = list()
pi1.tmp = phi[1]^(1:ncl[1])
Pi[[1]] = pi1.tmp/sum(pi1.tmp)

for (t in 2:nT) {
  nc1 = ncl[t-1]
  nc2 = ncl[t]
  ncmin = min(nc1, nc2)
  
  Pi[[t]] = matrix(NA, nrow = nc1, ncol = nc2)
  
  for (k in 1:ncmin){
    pik.tmp = rep(0, nc2)
    pik.tmp = phi[t] * nc2:1
    pik.tmp[k] = 0
    pik.tmp = pik.tmp/sum(pik.tmp)*phi[t]
    pik.tmp[k] = 1 - phi[t]
    Pi[[t]][k, ] = pik.tmp
  }
  
  if (nc1>nc2){
    pik.tmp = rep(0, nc2)
    pik.tmp = (nc2:1)^2
    pik.tmp = pik.tmp/sum(pik.tmp)
    Pi[[t]][(nc2+1):nc1, ] = matrix(pik.tmp, nrow = nc1-nc2, ncol = nc2, byrow = TRUE)
  }
}

pdf("Pi_sim.pdf", width = 12, height = 9)
for (t in 2:nT){
  plt = pheatmap::pheatmap(Pi[[t]], cluster_rows = F, cluster_cols = F, silent = T,
                           breaks = seq(0, 1, length = 101)) %>% ggplotify::as.ggplot()
  print(plt)
}
dev.off()

out = simdata(n = n, p = p, nT = nT, 
              ncl = ncl, Pi = Pi, 
              m = m, m.important = m.important, values = values, seed = 1)
C = out$C
Y = out$Y

library(ggalluvial)
Cflow = C %>% reshape2::melt(value.name = "cluster") %>% 
  mutate(freq = 1) %>% 
  mutate(cluster = factor(cluster, levels = 1:max(C)))

Cflow %>% ggplot(aes(x = Time, stratum = cluster, alluvium = Subject,
                     y = freq,
                     fill = cluster)) +
  geom_flow() +
  geom_stratum(color = "white") + 
  scale_x_discrete(expand = c(.1, .1)) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3)

ggsave("simulated_data_flow.pdf", width = 10, height = 6)

palette = c(pals::alphabet(26)[-(4:5)], pals::alphabet2(26), pals::glasbey(32)[-c(4, 18)])[1:max(sapply(values, max))] %>% unname
names(palette) = 1:max(sapply(values, max))

pdf("heatmap_sim.pdf", width = 12, height = 9)

for (t in 1:nT){
  plt = Y %>% 
    reshape2::melt(value.name = "value") %>% 
    mutate(value = factor(value)) %>% 
    left_join(reshape2::melt(C) %>% rename(Cluster = value), by = c("Subject", "Time")) %>%
    mutate(Subject = factor(Subject, levels = 1:n)) %>%
    filter(Time == t) %>% 
    ggplot(aes(x = Variable, y = Subject, fill = value)) +
    geom_raster() +
    scale_fill_manual(values = palette) +
    theme_minimal() +
    facet_wrap(~ Cluster, ncol = 4, scales = "free_y")
  print(plt)
}


dev.off()
