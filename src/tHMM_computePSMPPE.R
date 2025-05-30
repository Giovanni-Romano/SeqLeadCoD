# Packages ----
options(warning = 1)
suppressPackageStartupMessages(
  {
    library(tidyverse)
    library(mcclust) # comp.psm()
    library(salso)
  }
)

# Command line arguments ----
cmndargs = commandArgs(trailingOnly = TRUE)

# Load objects ----
TYPE = cmndargs[1] #"Gnedin"
SEED = cmndargs[2]
NAME = paste0(TYPE, SEED)
OUTFOLDER = "output/tHMM/"

## Model output ----
RDS = readRDS(paste0(OUTFOLDER, "res_", NAME, ".RDS"))
C = RDS$output$C

## Data ----
ghe = readRDS("data/rds/dataGHE.RDS") %>% unite("ID", c(CountryN, Sex), remove = FALSE) %>% 
  unite("IDshort", c(Country, Sex), remove = FALSE)

## Define dimensions ----
ages = levels(ghe$Age)
nages = length(ages)
years = 2000:2021
nyears = length(years)
transitions = paste(years[-nyears], years[-1], sep = "-")
IDs = ghe %>% pull(IDshort) %>% unique %>% sort
nIDs = length(IDs)
ncauses = length(levels(ghe$CauseS))
ntransitions = length(transitions)
niter = (RDS$input$ctr$ctr_mcmc$nchain + RDS$input$ctr$ctr_mcmc$nburnin)
nburnin = RDS$input$ctr$ctr_mcmc$nburnin

# Compute and save PSM and PPE ----
# PSM: Posterior Similarity Matrix
PSM = apply(C[ , , -(1:nburnin)], 2, function(x) comp.psm(t(x)), simplify = FALSE)
PSM = lapply(PSM, function(x){
  dimnames(x) = list("ID1" = IDs, "ID2" = IDs)
  x
})
names(PSM) = years
saveRDS(PSM, file = paste0("output/tHMM/PSM_", NAME, ".RDS"))
# PPE: Partition Point Estimate
PPE = list()
for (y in 1:nyears){
  cat(y, "\t")
  PPE[[y]] = salso(C[ , y, -(1:nburnin)] %>% t)
  gc()
}
dimnames(PPE) = list("ID" = IDs, "year" = years)
saveRDS(PPE, file = paste0("output/tHMM/salso_", RDS$input$ctr$ctr_mcmc$seed, ".RDS"))
