# Packages ----
options(warning = 1)
suppressPackageStartupMessages(
  {
    library(tidyverse)
    library(mcclust) # comp.psm()
    library(salso)
  }
)

source("src/utils.R")
source("src/tHMM_updates.R")

# Command line arguments ----
cmndargs = commandArgs(trailingOnly = TRUE)

# Load objects ----
SEX = "female" #cmndargs[1]
SEEDS = c(20010, 20019, 20146, 20148, 76137)
NAMES = paste0("Gnedin", SEEDS)
OUTFOLDER = paste0("output/tHMM/5yrs/PT_5/", SEX, "/")
IMGFOLDER = paste0("img/tHMM/5yrs/PT_5/", SEX, "/compare_runs/")
if (!dir.exists(IMGFOLDER)){
  dir.create(IMGFOLDER, recursive = T)
}

## Model output ----
RDS_list = lapply(paste0(OUTFOLDER, "res_", NAMES, ".RDS"), readRDS)
C_list = lapply(RDS_list, function(RDS) RDS$traces[[1]]$C)

## Data ----
ghe = readRDS(paste0("data/rds/GHEdf_", SEX, "_5yrs.RDS"))

## Define dimensions ----
ages = levels(ghe$Age)
nages = length(ages)
years = sort(unique(ghe$Year))
nyears = length(years)
transitions = paste(years[-nyears], years[-1], sep = "-")
countries = ghe %>% pull(CountryN) %>% unique %>% as.character
ncountries = length(countries)
ntransitions = length(transitions)
nruns = length(SEEDS)
nburnin = RDS_list[[1]]$ctr$ctr_mcmc$nburnin
niter = (RDS_list[[1]]$ctr$ctr_mcmc$nchain + nburnin)
thin = 10
idx = seq(nburnin+thin, niter, by = thin)

# Compute and save PSM and PPE ----
# PSM: Posterior Similarity Matrix
PSM = lapply(C_list, function(C) {
  tmp  = apply(C[ , , idx], 2, function(x) comp.psm(t(x)), simplify = FALSE) %>% 
    abind::abind(., along = 3)
  dimnames(tmp) = list("Country1" = countries, "Country2" = countries, "Year"  = years)
  tmp
})
names(PSM) = NAMES
# saveRDS(PSM, file = paste0(OUTFOLDER, "PSM_all.RDS"))


# PPE: Partition Point Estimate
PPE = salso_out = list()
for (run in seq_along(SEEDS)){
  cat("Run ", run, ": ")
  tmp_PPE = matrix(NA, ncountries, nyears)
  tmp_out = list()
  C = C_list[[run]]
  for (y in 1:nyears){
    cat(y, "\t")
    tmp_out[[y]] = salso(C[ , y, idx] %>% t, 
                         loss = VI(),
                         maxNClusters = 90, maxZealousAttempts = 30,
                         nRuns = 100,
                         nCores = 0, # use all the cores
    )
    tmp_PPE[ , y] = tmp_out[[y]] %>% as.integer()
    gc()
  }
  dimnames(tmp_PPE) = list("Country" = countries, "year" = years)
  names(tmp_out) = years
  PPE[[run]] = tmp_PPE
  salso_out[[run]] = tmp_out
  cat("\n")
}
names(PPE) = names(salso_out) = SEEDS
# saveRDS(PPE, file = paste0(OUTFOLDER, "PPE_all.RDS"))

# SALSO summary ----
salso_expLoss = sapply(salso_out, function(out) 
  sapply(out, function(x) attr(x, "info")["expectedLoss"] %>% as.numeric))

# Some summary plots ----
# Number of clusters
ncl = lapply(C_list, function(C) apply(C[ , , -1], c(2:3), max)) %>% abind::abind(along = 3)
dimnames(ncl) = list("Year" = years, "Iterations" = 2:niter, "Run" = SEEDS)

pdf(paste0(IMGFOLDER, "ncl.pdf"), width = 12, height = 12)
par(mfrow = c(5, 5))
for (y in 1:nyears){
  min_ncl = min(ncl[y, ,])
  max_ncl = max(ncl[y, , ])
  for (run in seq_along(SEEDS)){
    est_ncl = PPE[[run]][ , y] %>% max
    
    plot(ncl[y, -(1:15), run], type = "l", col = 2,
         ylim = c(min_ncl, max_ncl),
         ylab = "", xlab = "Iter", main = paste0(years[y], " - Nr. Cl. (Run ", SEEDS[run], ")"),
         yaxt = "n")
    axis(side = 2, at = min_ncl:max_ncl)
    abline(v = nburnin, lty = 2, col = "blue")
    abline(h = est_ncl, lty = 2, col = "gold")
    legend("topright", lty = c(1, 2, 2), col = c(2, "gold", "blue"), 
           legend = c("Nr. Cl.", "Est.",  "burnin"))
  }
  
}
dev.off()

# NMI consecutive partitions
NMI_cons = array(NA, dim = c(nruns, nyears, niter/thin))
for (r in 1:nruns){
  cat("Run ", r, ": ")
  for (y in 1:nyears){
    cat(years[y], "\t")
    for (i in 2:(niter/thin)){
      NMI_cons[r, y, i-1] = aricode::NMI(C_list[[r]][, y, i*thin], C_list[[r]][, y, (i-1)*thin])
    }
  }
  cat("\n")
}
dimnames(NMI_cons) = list("Run" = SEEDS, "Year" = years, "Iter" = 1:(niter/thin))

pdf(paste0(IMGFOLDER, "NMI_chain.pdf"), width=12, height=12)
par(mfrow=c(5, 5))
for (y in 1:nyears){
  for (r in 1:nruns){
    plot(seq(thin, niter, by = thin), NMI_cons[r, y, ], type = "l", ylim = c(0.5, 1), xlab = "Iter", ylab = "", 
         main = paste0(years[y], " - NMI at different lags"))
    lines(seq(50, niter, by = 50), NMI_cons[r, y, seq(1, dim(NMI_cons)[3], by = 5)], col = 2)
    legend("bottomright", lty = 1, col = 1:2, legend = c(paste("thin = ", 10), "thin = 50"))
  }
}
dev.off()

# Alpha e gamma agreement
pdf(paste0(IMGFOLDER, "alphagamma.pdf"), width=15, height=12)
par(mfrow=c(4, 5))
for (y in 2:nyears){
  for (r in 1:nruns){
    plot(seq(thin, niter, by = thin), RDS_list[[r]]$traces[[1]]$alpha[y, seq(thin, niter, by = thin)], 
         type = "l", ylim = c(0, 1), 
         xlab = "Iter", ylab = "", 
         main = paste0(years[y], "- seed ", SEEDS[r], " (thinned by ", thin, ")"))
    lines(seq(thin, niter, by = thin), RDS_list[[r]]$traces[[1]]$gamma[y, seq(thin, niter, by = thin)]/ncountries, col = 2)
    abline(v = nburnin, lty = 2, col = "blue")
    legend("bottomright", lty = c(1, 1, 2), col = c(1:2, "blue"), 
           legend = c("alpha", "% gamma = 1", "burnin"))
  }
}
dev.off()


# Loglikelihood
log_lik = sapply(RDS_list, function(x) x$traces[[1]]$log_lik)
dimnames(log_lik) = list("Iteration" = 1:niter, "Run" = SEEDS)

plt_logLik = log_lik %>% reshape2::melt() %>% 
  mutate(Run = factor(Run)) %>% 
  filter(Iteration %in% seq(thin, niter, by = thin)) %>% 
  ggplot() + 
  geom_vline(aes(xintercept = nburnin), lty = 2, lwd = 1) +
  geom_line(aes(x = Iteration, y = value, col = Run)) +
  labs(y = "logLik", title = "logLik during MCMC scheme")
ggsave(path = IMGFOLDER,
       plot = plt_logLik,
       filename = "loglik.pdf",
       width = 12, height = 6)


best_run = which(SEEDS == 20019)

# pheatmap PSM
PSM_plotlist = list()
for (y in 1:nyears){
  ph = pheatmap::pheatmap(PSM[[best_run]][ , , y], treeheight_row=0, treeheight_col=0,
                          clustering_distance_rows = as.dist(1-PSM[[best_run]][ , , y]),
                          clustering_distance_cols = as.dist(1-PSM[[best_run]][ , , y]),
                          clustering_method = "ward.D",
                          silent = T)
  order_y = ph$tree_row$order
  
  for (r in 1:nruns){
    run = c(1, 2, 3, 4, 5)[r]
    df_ann = as.data.frame(PPE[[run]][ , y])
    colnames(df_ann) = "Cluster"
    
    PSM_plotlist[[(y-1)*(nruns) + r]] = pheatmap::pheatmap(PSM[[run]][order_y , order_y, y], 
                                                           cluster_rows = F, cluster_cols = F, 
                                                           show_rownames = F, show_colnames = F,
                                                           border_color = "grey60",
                                                           clustering_method = "ward.D",
                                                           annotation_row = df_ann, annotation_col = df_ann,
                                                           annotation_colors = list(Cluster = pals::glasbey(max(df_ann))),
                                                           annotation_legend = FALSE,
                                                           annotation_names_col = FALSE,
                                                           annotation_names_row = FALSE,
                                                           
                                                           main = paste0("Year ", years[y], " - Run ", SEEDS[run]),
                                                           
                                                           legend = FALSE,
                                                           silent = T) %>% ggplotify::as.ggplot()
  }
  
  # PSM_plotlist[[(y-1)*(nruns + 1) + nruns + 1]] = ggplot() + theme_void()
}


ggsave(
  path = IMGFOLDER,
  filename = "PSM.pdf", 
  plot = gridExtra::marrangeGrob(PSM_plotlist, nrow = 1, ncol = 5, top = NULL),
  width = 20, height = 4)



# NMI PPE
NMI_PPE = array(NA, c(nyears, nruns, nruns))
dimnames(NMI_PPE) = list(Year = years, Run1 = SEEDS, Run2 = SEEDS)
for (r1 in 1:length(RDS_list)){
  for (r2 in 1:length(RDS_list)){
    for (y in 1:nyears)
      NMI_PPE[y, r1, r2] = aricode::NMI(PPE[[r1]][ , y], PPE[[r2]][ , y])
  }
}

NMI_plot = NMI_PPE %>% 
  reshape2::melt() %>% 
  mutate(Run1 = factor(Run1), Run2 = factor(Run2)) %>% 
  ggplot(aes(x = Run1, y = Run2))+
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_viridis_c("NMI", limits = c(0, 1)) +
  facet_wrap(~Year) +
  labs(title = "NMI between partition point estimates across runs (with parallel tempering)")

ggsave(
  path = IMGFOLDER,
  filename = "NMI_PPE.pdf", 
  plot = NMI_plot,
  width = 9, height = 6)




# --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Evaluation of PT ----
nchains = RDS_list[[1]]$traces %>% length
diagnostics = lapply(RDS_list, function(x) x$diagnostics)

att_swap = sapply(diagnostics, function(d) d$attempted_swap)
succ_swap = sapply(diagnostics, function(d) d$succesfull_swap)
dimnames(att_swap) = dimnames(succ_swap) = list(Chain = paste0("Chain ", 1:nchains), Run = SEEDS)


prob_move = lapply(diagnostics, function(d) {
  tmp = t(d$prob_move[ , 2:niter])
  dimnames(tmp) = list(Iteration = 2:niter, Swap = paste(2:nchains, 1:(nchains-1), sep = "-"))
  tmp %>% as.data.frame() %>% rownames_to_column("Iteration") %>% pivot_longer(-Iteration, names_to = "Chain", values_to = "Prob_move")
})
names(prob_move) = SEEDS

prob_rolling = lapply(diagnostics, function(d) {
  tmp = apply(diagnostics[[1]]$prob_move[ , 2:niter], 1, function(x) cumsum(x)/(1:(niter-1)))
  dimnames(tmp) = list(Iteration = 2:niter, Swap = paste(2:nchains, 1:(nchains-1), sep = "-"))
  tmp %>% as.data.frame() %>% rownames_to_column("Iteration") %>% pivot_longer(-Iteration, names_to = "Chain", values_to = "Prob_rolling")
})
names(prob_rolling) = SEEDS

rho = lapply(diagnostics, function(d) {
  tmp = t(d$rho[ , 2:niter])
  dimnames(tmp) = list(Iteration = 2:niter, Swap = paste(2:nchains, 1:(nchains-1), sep = "-"))
  tmp %>% as.data.frame() %>% rownames_to_column("Iteration") %>% pivot_longer(-Iteration, names_to = "Chain", values_to = "Rho")
})
names(rho) = SEEDS

temperature =  lapply(diagnostics, function(d) {
  tmp = t(apply(d$rho[ , 2:niter], 2, update_temp_schedule))
  dimnames(tmp) = list(Iteration = 2:niter, Chain = 1:nchains)
  tmp %>% as.data.frame() %>% rownames_to_column("Iteration") %>% pivot_longer(-Iteration, names_to = "Chain", values_to = "Temperature")
})
names(temperature) = SEEDS

tempering_df = bind_rows(prob_move, .id = "Run") %>%
  full_join(bind_rows(prob_rolling, .id = "Run"), by = c("Iteration", "Chain", "Run")) %>% 
  full_join(bind_rows(rho, .id = "Run"), by = c("Iteration", "Chain", "Run")) %>%
  mutate(Iteration = as.integer(Iteration))

prob_target = RDS_list[[1]]$ctr$ctr_swap$prob_target

# Plot of temperatures
temperatures_plot = bind_rows(temperature, .id = "Run") %>%
  mutate(Iteration = as.integer(Iteration)) %>%
  filter(Iteration %in% seq(1, niter, by = thin)) %>% 
  ggplot(aes(x = Iteration, y = Temperature)) + 
  geom_line(aes(color = Chain)) +
  facet_grid(~Run) +
  labs(title = "Evolution of temperature across chains and runs [thinning 10]") +
  scale_y_continuous(limits = c(0.8, 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(path = IMGFOLDER, 
       filename = "temperatures.pdf", 
       plot = temperatures_plot, width = 12, height = 4)

# Plot of prob. of swap and its rolling mean
tempering_plot = tempering_df %>%
  select(-Rho) %>% 
  pivot_longer(c("Prob_move", "Prob_rolling"), names_to = "Quantity", values_to = "value") %>%
  filter(Iteration %in% seq(1, niter, by = thin)) %>% 
  ggplot() +
  geom_vline(aes(xintercept = nburnin), lty = 2) +
  geom_line(aes(x = Iteration, y = value, col = Quantity, alpha = Quantity)) +
  geom_hline(aes(yintercept = prob_target), col = "gold", lty = 2, lwd = 1) +
  scale_alpha_discrete(limits = "Prob_move", range = c(0.3, 0.3), guide = "none") +
  labs(title = "Iter-specific and rolling average probability of swap [target in yellow, burnin in black]",
       y = "Prob swap") +
  facet_grid(Chain~Run)

ggsave(path = IMGFOLDER, 
       filename = "tempering_prob.pdf", 
       plot = tempering_plot, width = 16, height = 9)


# Plot of evolution chains positions
chains_track = lapply(diagnostics, function(d) {
  tmp = t(d$chains_tracking)
  dimnames(tmp) = list(Iteration = 1:niter, Chain = 1:nchains)
  tmp %>% as.data.frame() %>% rownames_to_column("Iteration") %>% pivot_longer(-Iteration, names_to = "Chain", values_to = "Position")
})
names(chains_track) = SEEDS

position_evolution_plot = bind_rows(chains_track, .id = "Run") %>% 
  mutate(Iteration = as.integer(Iteration)) %>% 
  filter(Iteration %in% seq(1, niter, by = thin)) %>%
  ggplot(aes(x = Iteration, y = Position)) +
  geom_line()+
  facet_grid(Run ~ Chain) +
  labs(title = "Position of each parallel chain (columns) during the PT-MCMC scheme [thinning by 10]")

ggsave(path = IMGFOLDER, 
       filename = "trajectories_PT.pdf", 
       plot = position_evolution_plot, width = 16, height = 6)

# Plot of proportion of iterations spent in each position
prop_time_plot = bind_rows(chains_track, .id = "Run") %>% 
  mutate(Iteration = as.integer(Iteration)) %>% 
  filter(Iteration %in% seq(nburnin + thin, niter, by = thin)) %>% 
  group_by(Run, Chain, Position) %>% 
  summarise(Prop_time = n()/((niter-nburnin)/thin)) %>% 
  ggplot(aes(x = Chain, y = Position)) +
  geom_tile(aes(fill = Prop_time)) +
  geom_text(aes(label = round(Prop_time *100)), col = "gray95") +
  scale_fill_viridis_c("Iter %", limits = c(0, 1)) +
  facet_wrap(~Run) +
  labs(title = "% of each parallel chain (x-axis) spent at each position (y-axis) [thinning by 10, after 5k burnin]")


ggsave(path = IMGFOLDER, 
       filename = "iterperc_PT.pdf", 
       plot = prop_time_plot, width = 12, height = 6)



# Compare PPE with non-PT ones
PPE_nonPT = readRDS(paste0("output/tHMM/5yrs/no_PT/", SEX, "/PPE_all.RDS"))
PPE_PT8 = readRDS(paste0("output/tHMM/5yrs/PT/", SEX, "/PPE_all.RDS"))

PPE_all = append(append(PPE, PPE_PT8), PPE_nonPT)
names(PPE_all) = c(paste0(SEEDS, "_PT5"), paste0(SEEDS, "_PT8"), paste0(SEEDS, "_nonPT"))
NMI_PPEall = array(NA, dim = c(nyears, 3*nruns, 3*nruns), 
                   dimnames = list(Year = years, Run1 = names(PPE_all), Run2 = names(PPE_all)))
for (r1 in 1:length(PPE_all)){
  for (r2 in 1:length(PPE_all)){
    for (y in 1:nyears){
      NMI_PPEall[y, r1, r2] = aricode::NMI(PPE_all[[r1]][ , y], PPE_all[[r2]][ , y])
    }
  }
}


NMIall_plotlist = lapply(years, function(yyyy) 
  NMI_PPEall %>% 
    reshape2::melt() %>% 
    filter(Year == yyyy) %>% 
    mutate(Run1 = factor(Run1), Run2 = factor(Run2),
           Group1 = str_split_i(Run1, "_", 2),
           Group2 = str_split_i(Run2, "_", 2),) %>% 
    ggplot(aes(x = Run1, y = Run2))+
    geom_tile(aes(fill = value)) +
    geom_text(aes(label = round(value, 2))) +
    scale_fill_viridis_c("NMI", limits = c(0.5, 1)) +
    ggh4x::facet_nested(Group2 ~ Year+Group1, scales = "free") +
    # facet_grid(Group2~Group1, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = "NMI between partition point est. across runs (both with and without PT)")
)

ggsave(
  path = IMGFOLDER,
  filename = "NMI_PPEall.pdf", 
  plot = gridExtra::grid.arrange(grobs = NMIall_plotlist, nrow = 2),  
  width = 32*0.75, height = 18*0.75)
