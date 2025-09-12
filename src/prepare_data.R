suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# My regions
{
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  postsoviet = c("ARM", "AZE", "BLR", "EST", "GEO", "KAZ", "KGZ", "LVA", 
                 "LTU", "MDA", "RUS", "TJK", "TKM", "UKR", "UZB") 
  myregions = world %>% 
    as.data.frame() %>%
    filter(iso_a3_eh != "-99",
           type != "Dependency") %>% 
    select(iso_a3_eh, geounit, region_wb, subregion) %>% 
    mutate(postsoviet = if_else(iso_a3_eh %in% postsoviet, "Yes", "No")) %>%
    mutate(myregions = case_when(postsoviet == "Yes" ~ "Post-Sovietic",
                                 iso_a3_eh == "TUR" ~ "Middle East & North Africa",
                                 region_wb == "Europe & Central Asia" ~ "Europe (Non P-S)",
                                 iso_a3_eh %in% c("USA", "CAN", "AUS", "NZL") ~ "Western Offshoots",
                                 TRUE ~ region_wb)) %>% 
    mutate(myregions = factor(myregions,
                              levels = c("Western Offshoots", "Latin America & Caribbean",
                                         "Sub-Saharan Africa", "Middle East & North Africa", 
                                         "Europe (Non P-S)", "Post-Sovietic",
                                         "South Asia", "East Asia & Pacific")),
           Region = case_when(myregions == "Western Offshoots" ~ "WO",
                              myregions == "Latin America & Caribbean" ~ "LAC",
                              myregions == "Sub-Saharan Africa" ~ "SSA", 
                              myregions == "Middle East & North Africa" ~ "MENA", 
                              myregions == "Europe (Non P-S)" ~ "EU", 
                              myregions == "Post-Sovietic" ~ "PS",
                              myregions == "South Asia" ~ "SA", 
                              myregions == "East Asia & Pacific" ~ "EAP")) %>%
    select(iso_a3_eh, geounit, myregions, Region) %>%
    rename(CountryN = geounit, RegionN = myregions) %>% 
    unique() 
}

# Causes with short name
causes = readr::read_csv2("data/raw/causes.csv", 
                          col_types = cols(.default = "character")) %>% 
  select(CauseC:CauseT)

# Import CSV and join with Causes and Regions info
raw_df = read_csv2("data/raw/WHO_GHE_Top1.csv", 
                   col_types = cols(DIM_GHECAUSE_CODE = col_character())) %>% 
  select(-VAL_DTHS_RATE100K_NUMERIC) %>% 
  rename_with(~ gsub("DIM_", "", .), starts_with("DIM_")) %>% 
  rename(Country = COUNTRY_CODE, Year = YEAR_CODE,
         CauseC = GHECAUSE_CODE, CauseT = GHECAUSE_TITLE,
         Age = AGEGROUP_CODE, Sex = SEX_CODE) %>% 
  filter(!(Age %in% c("D0T27", "M1T11", "TOTAL"))) %>% 
  left_join(myregions, by = c("Country" = "iso_a3_eh")) %>% 
  left_join(causes, by = c("CauseC", "CauseT")) %>% 
  mutate(Sex = if_else(Sex == "FEMALE", "F", "M"),
         Age = case_when(Age == "YGE_85" ~ "85+",
                         .default = gsub("Y(\\d{1,2})T(\\d{1,2})", "\\1-\\2", Age)),
         mutate(across(where(is.character), as.factor))) %>% 
  select(Year, Country, CountryN, Age, Sex, CauseC, CauseS, CauseT, Region, RegionN) %>% 
  mutate(Age = factor(Age, levels = c("0-1", "1-4", "5-9", 
                                      "10-14", "15-19", "20-24", 
                                      "25-29", "30-34", "35-39", 
                                      "40-44", "45-49", "50-54", 
                                      "55-59", "60-64", "65-69", 
                                      "70-74", "75-79", "80-84", 
                                      "85+")))

# Split female and male
GHEdf_male = filter(raw_df, Sex == "M")
GHEdf_female = filter(raw_df, Sex == "F")
nrow(GHEdf_female) + nrow(GHEdf_male); nrow(raw_df)

saveRDS(GHEdf_male, "data/rds/GHEdf_male.RDS")
saveRDS(GHEdf_female, "data/rds/GHEdf_female.RDS")
