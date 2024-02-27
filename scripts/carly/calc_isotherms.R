library(here)
library(tidyverse)
library(data.table)
library(purrr)
library(broom)

# Read in necessary data
poldat.stats <- readRDS(here("processed-data","poldat.stats.rds"))
eqdat.stats <- readRDS(here("processed-data","eqdat.stats.rds")) 
soda.stats <- readRDS(here("processed-data","soda_stats.rds")) 
# oisst.stats <- readRDS(here("processed-data","oisst_stats.rds"))
# hadisst.stats <- readRDS(here("processed-data","hadisst_stats.rds"))

# oisst.neus <- readRDS(here("processed-data","oisst_neus.rds"))
# hadisst.neus <- readRDS(here("processed-data","hadisst_neus.rds"))
soda.neus <- readRDS(here("processed-data","soda_neus.rds"))

# Linear model functions for SODA 
lm.soda <- soda.neus %>% 
  nest(-year_match) %>% 
  mutate(
    model = purrr::map(data, ~lm(surftemp ~ y, data = .x)), 
    tidymodel = purrr::map(model, tidy)
  ) %>% 
  unnest(tidymodel, .drop=TRUE)
get.soda.lat <- function(temp, est_year) {
  tmp <- lm.soda[lm.soda$year_match== est_year,]
  out <- (temp-tmp[[1,3]])/tmp[[2,3]]
  return(out)
}
get.soda.temp <- function(x, est_year) {
  tmp <- lm.soda[lm.soda$year_match== est_year,]
  out <- x*tmp[[2,3]]+tmp[[1,3]]
  return(out)
}

# poleward edge isotherm 
poldat.iso.soda.prep <- poldat.stats %>% 
  filter(est_year >= min(soda.neus$year_match)) %>% 
  group_by(comname) %>% 
  arrange(est_year) %>% 
  slice(1:3) %>% # get first 3 years that species is observed
  ungroup()

# calculate temperature values for first three years when species was observed, using the linear model for each dataset and each year 
poldat.iso.soda.prep2 <- NULL

for(i in unique(poldat.iso.soda.prep$comname)){
  # year.range <- poldat.iso.soda.prep[poldat.iso.soda.prep$comname==i,]$est_year
  # year.range <- as.numeric(unique(year.range))
  for(j in year.range){
    spp.lat95 <- poldat.iso.soda.prep[poldat.iso.soda.prep$est_year==j & poldat.iso.soda.prep$comname==i,]$spp.lat95
    est.soda <- get.soda.temp(spp.lat95, j)
    out <- cbind(i, j, spp.lat95, est.soda)
    poldat.iso.soda.prep2 <- rbind(out, poldat.iso.soda.prep2)} 
  rm(i, j, year.range, spp.lat95, est.soda, out)}

# tidy data frame to get the baseline temperature values for each species
poldat.iso.soda.prep3 <- poldat.iso.soda.prep2 %>% 
  as.data.frame() %>% 
  rename(comname=i, est_year=j) %>% 
  group_by(comname) %>% 
  mutate(est.soda = as.numeric(as.character(est.soda)), 
         est.edge.temp.soda = mean(est.soda)) %>% 
  dplyr::select(comname, est.edge.temp.soda) %>% 
  ungroup() %>% 
  distinct() %>% 
  full_join(poldat.stats, by="comname") %>% 
  filter(est_year >= min(soda.neus$year_match))

# calculate future latitudes for the species-specific isotherms

poldat.iso.soda.prep4 <- NULL

for(i in unique(poldat.iso.soda.prep3$latinname)) {
  year.range <- poldat.iso.soda.prep3[poldat.iso.soda.prep3$latinname==i,]$year 
  est.edge.temp.soda <- poldat.iso.soda.prep3 %>% filter(latinname==i) %>% group_by(latinname) %>% dplyr::select(est.edge.temp.soda) %>% distinct() %>% pull()
  for(j in year.range) {
    est.edge.lat.soda <- get.soda.lat(est.edge.temp.soda, j)
    out <- cbind(i, j, est.edge.lat.soda, est.edge.temp.soda)
    poldat.iso.soda.prep4 <- rbind(out, poldat.iso.soda.prep4)
  }
  rm(i, j, year.range, est.edge.temp.soda, est.edge.lat.soda, out)
}

poldat.iso.soda.prep5 <- 
  poldat.iso.soda.prep4 %>% 
  as.data.frame() %>% 
  rename(latinname=i, year=j) %>% 
  mutate(latinname = as.character(latinname),
         year = as.integer(as.character(year)),
         est.edge.lat.soda = as.numeric(as.character(est.edge.lat.soda))) 


poldat.stats.iso <- poldat.stats %>% 
  # left_join(poldat.iso.had.prep5, by=c('year','latinname')) %>% 
  left_join(poldat.iso.soda.prep5, by=c('year','latinname')) %>% 
  # left_join(poldat.assembl.iso.had, by=c('year','assemblage.lat95')) %>% 
  left_join(poldat.assembl.iso.soda, by=c('year','assemblage.lat95'))

write_rds(poldat.stats.iso, here("processed-data","poldat.stats.iso.rds"))