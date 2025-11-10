# Packages ----
options(warning = 1)
# source("renv/activate.R")
suppressPackageStartupMessages(
  {
    library(tidyverse)
    library(mcclust) # comp.psm()
    library(salso)
  }
)

source("src/utils.R")

# Command line arguments ----
cmndargs = commandArgs(trailingOnly = TRUE)

# Load objects ----
SEX = "male" #cmndargs[1]
SEEDS = c(20010, 20019, 20146, 20148, 76137)
NAMES = paste0("Gnedin", SEEDS)
OUTFOLDER = paste0("output/tHMM/5yrs/", SEX, "/")
IMGFOLDER = paste0("img/tHMM/5yrs/", SEX, "/compare_runs/")
if (!dir.exists(IMGFOLDER)){
  dir.create(IMGFOLDER, recursive = T)
}

## Model output ----
RDS_list = lapply(paste0(OUTFOLDER, "res_", NAMES, ".RDS"), readRDS)
C_list = lapply(RDS_list, function(RDS) RDS$output$C)

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
nburnin = RDS_list[[1]]$input$ctr$ctr_mcmc$nburnin
niter = (RDS_list[[1]]$input$ctr$ctr_mcmc$nchain + nburnin)
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
cat("Done! \n")

# SALSO summary ----
salso_expLoss = sapply(salso_out, function(out) 
  sapply(out, function(x) attr(x, "info")["expectedLoss"] %>% as.numeric))

# Some summary plots ----
# Number of clusters
ncl = lapply(RDS_list, function(RDS) apply(RDS$output$C, c(2:3), max)) %>% abind::abind(along = 3)
dimnames(ncl) = list("Year" = years, "Iterations" = 1:niter, "Run" = SEEDS)

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
      NMI_cons[r, y, i-1] = aricode::NMI(RDS_list[[r]]$output$C[, y, i*thin], RDS_list[[r]]$output$C[, y, (i-1)*thin])
    }
  }
  cat("\n")
}
dimnames(NMI_cons) = list("Run" = SEEDS, "Year" = years, "Iter" = 1:(niter/thin))

pdf(paste0(IMGFOLDER, "NMI_chain.pdf"), width=12, height=12)
par(mfrow=c(5, 5))
for (y in 1:nyears){
  for (r in 1:nruns){
    plot(seq(thin, niter, by = thin), NMI_cons[r, y, ], type = "l", ylim = c(0.7, 1), xlab = "Iter", ylab = "", 
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
    plot(seq(thin, niter, by = thin), RDS_list[[r]]$output$alpha[y, seq(thin, niter, by = thin)], 
         type = "l", ylim = c(0, 1), 
         xlab = "Iter", ylab = "", 
         main = paste0(years[y], "- seed ", SEEDS[r], " (thinned by ", thin, ")"))
    lines(seq(thin, niter, by = thin), RDS_list[[r]]$output$gamma[y, seq(thin, niter, by = thin)], col = 2)
    abline(v = nburnin, lty = 2, col = "blue")
    legend("bottomright", lty = c(1, 1, 2), col = c(1:2, "blue"), 
           legend = c("alpha", "% gamma = 1", "burnin"))
  }
}
dev.off()


# Loglikelihood
log_lik = lapply(RDS_list, function(x) x$output$loglik) %>% abind::abind(along = 3)
dimnames(log_lik) = list("Year" = years, "Iteration" = 1:niter, "Run" = SEEDS)

plt_logLik = bind_rows(log_lik %>% reshape2::melt() %>%
                         mutate(Year = factor(Year, levels = c(years, "Total"))),
                       log_lik %>% reshape2::melt() %>% 
                         group_by(Run, Iteration) %>%
                         summarise(value = sum(value), .groups = "drop") %>%
                         mutate(Year = "Total") %>% 
                         mutate(Year = factor(Year, levels = c(years, "Total")))) %>%
  mutate(Run = factor(Run)) %>% 
  filter(Iteration %in% seq(thin, niter, by = thin)) %>% 
  ggplot() + 
  geom_vline(aes(xintercept = nburnin), lty = 2, lwd = 1) +
  geom_line(aes(x = Iteration, y = value, col = Run)) +
  facet_wrap(~Year, scales = "free_y") +
  labs(y = "logLik", title = "logLik during MCMC scheme")
ggsave(path = IMGFOLDER,
       plot = plt_logLik,
       filename = "loglik.pdf",
       width = 12, height = 6)


best_run = which(SEEDS == 20148)

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
    run = c(4, 3, 5, 2, 1)[r]
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

NMI_PPE %>% 
  reshape2::melt() %>% 
  ggplot() +
  geom_tile(aes(x = Run1, y =Run2))
