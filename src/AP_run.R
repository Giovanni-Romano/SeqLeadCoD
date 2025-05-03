library(mcclust)
library(mcclust.ext)
# library(AntMAN)
library(klaR)
library(TraMineR)
library(readr)
suppressPackageStartupMessages(library(tidyverse))
source('src/ArgientoPaci/code/gibbs_sampler2.R', echo=TRUE)
source('src/ArgientoPaci/code/complement_functions.R')




#GHE data ----
causes = read_csv2("data/raw/causes.csv", col_types = cols(.default = "character"))
countries = read_csv2("data/raw/countries.csv", col_types = cols(.default = "character")) %>% 
  mutate(ParentCode = gsub("R", "", ParentCode) )
ghe = read_csv2("data/raw/WHO_GHE_Top1.csv", 
                col_types = cols(DIM_GHECAUSE_CODE = col_character())) %>% 
  select(-VAL_DTHS_RATE100K_NUMERIC) %>% 
  rename_with(~ gsub("DIM_", "", .), starts_with("DIM_")) %>% 
  rename(Country = COUNTRY_CODE, Year = YEAR_CODE,
         CauseC = GHECAUSE_CODE, CauseT = GHECAUSE_TITLE,
         Age = AGEGROUP_CODE, Sex = SEX_CODE) %>% 
  filter(!(Age %in% c("D0T27", "M1T11", "TOTAL"))) %>% 
  left_join(causes, by = c("CauseC", "CauseT")) %>%
  left_join(countries %>% 
              rename(CountryN = Title,
                     Region = ParentCode) %>% 
              select(Code, CountryN, Region), by = c("Country" = "Code")) %>% 
  mutate(Sex = if_else(Sex == "FEMALE", "F", "M"),
         Age = case_when(Age == "YGE_85" ~ "85+",
                         .default = gsub("Y(\\d{1,2})T(\\d{1,2})", "\\1-\\2", Age)),
         mutate(across(where(is.character), as.factor))) %>% 
  select(Year, Country, CountryN, Age, Sex, CauseC, CauseS, CauseT, Region)
ghe$Age = factor(ghe$Age, levels = c("0-1", "1-4", "5-9", 
                                     "10-14", "15-19", "20-24", 
                                     "25-29", "30-34", "35-39", 
                                     "40-44", "45-49", "50-54", 
                                     "55-59", "60-64", "65-69", 
                                     "70-74", "75-79", "80-84", 
                                     "85+"))

ghe.listchr = lapply(split(ghe %>% select(Year, CountryN, Sex, Age, CauseS), ghe$Year),
                     function(df) df %>% select(-Year) %>% 
                       pivot_wider(names_from = Age, values_from = CauseS) %>%
                       unite("ID", c(CountryN, Sex)) %>% column_to_rownames("ID") %>% 
                       select(levels(ghe$Age)) %>% 
                       as.matrix())
ghe.lab = levels(ghe$CauseS)

ghe.list = replicate(length(ghe.listchr), matrix(0, nrow = nrow(ghe.listchr[[1]]), ncol = ncol(ghe.listchr[[1]])), 
                     simplify = F)
for (j in seq_along(ghe.list)) {
  for (i in seq_along(ghe.lab)) {
    ghe.list[[j]][ghe.listchr[[j]] == ghe.lab[i]] = i
  }
  
  rownames(ghe.list[[j]]) = rownames(ghe.listchr[[j]])
  colnames(ghe.list[[j]]) = colnames(ghe.listchr[[j]])
}

## Select single year ----
year = commandArgs(trailingOnly = TRUE) %>% as.numeric()
ghe.data = ghe.list[[year]]
ages=colnames(ghe.data)
p = ncol(ghe.data)
n = nrow(ghe.data)
m = rep(length(ghe.lab), p)

###########################################################################

# number of clusters a priori
Kstar = 15

# prior hyperparameters
Lambda = Kstar
gam = 1
prior = AM_prior_K_Pois(n=n, gam, Lambda = Lambda)
u = rep(55, p)
v = rep(1e-3, p)

## initial values
k.init = 20
M.na.init = 30
set.seed(20146)
C.init = sample(1:k.init,n,replace=T)
cent.init = matrix(sample(1:m[1],k.init+M.na.init,replace=T),nrow = (k.init+M.na.init), ncol=p)
sigma.init = matrix(0.5,nrow = (k.init+M.na.init), ncol=p)


### RUN GIBBS SAMPLER ###
set.seed(20146)
sim_ghe = gibbs_mix_con(G=1000,
                        burnin = 1000,
                        data=ghe.data,
                        u=u,v=v,
                        Lambda = Lambda,
                        gam = gam,
                        C.init = C.init,
                        k.init = k.init,
                        cent.init = cent.init,
                        sigma.init = sigma.init,
                        M.na.init =M.na.init,
                        M.max=100)
# posterior K
post_k = table(sim_ghe$k[2:1002])/length(2:1002)

# posterior similarity matrix
psm = comp.psm(sim_ghe$C[2:1002,])

# estimated clusters
pred_VI = minVI(psm)$cl
table(pred_VI)


###### parameters estimation #########
source("src/ArgientoPaci/code/gibbs_param_estim2.R")
Kest = length(unique(pred_VI))
cent.last = sigma.last = matrix(NA, Kest, p)
for(i in 1:Kest){
  cent.last[i,] = sim_ghe$Cent[[i]][nrow(sim_ghe$Cent[[i]]),]
  sigma.last[i,] = sim_ghe$Sigma[[i]][nrow(sim_ghe$Cent[[i]]),]
}

set.seed(5)
gibbs_param = gibbs_param_estim(G=1000,
                                gam=gam,
                                data=ghe.data,
                                u=u,
                                v=v,
                                estim_Truth = pred_VI,
                                cent.init = cent.last,
                                sigma.init = sigma.last,
                                Sm.init = sim_ghe$Sm[[length(sim_ghe$Sm)]][1:Kest]
)

