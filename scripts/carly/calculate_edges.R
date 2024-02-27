library(here)
library(tidyverse)
library(sf)
library(Hmisc)

coastdistdat <- readRDS(here("processed-data","coastdistdat.rds"))
poldat <- readRDS(here("processed-data","poldat.rds"))
eqdat <- readRDS(here("processed-data","eqdat.rds"))

# NOTE that this matches the points of the survey to points on the coast using simple distance minimization. 
# This means that a point that is close to a coastline to the north and slightly further offshore to the west will be assigned a coastal distance further north than if I had just measured the latitude of the observation and matched it to the coast. 
# However, I think this is more accurate because it assigns observations to the closest part of the shelf. 

get_length <- function(lon, lat, distdf) {
  tmp <- distdf %>% 
    mutate(abs.diff.x2 = abs(x-lon)^2,
           abs.diff.y2 = abs(y-lat)^2,
           abs.diff.xy = sqrt(abs.diff.x2 + abs.diff.y2
           )) %>% 
    filter(abs.diff.xy == min(abs.diff.xy)) %>% 
    dplyr::select(lengthfromhere) %>% 
    pull()
  return(tmp)
}

poldat.stats <- poldat %>% 
  rename("lat" = "decdeg_beglat",
         "lon" = "decdeg_beglon") %>%
  rowwise() %>% 
  mutate(coastdist_km = (get_length(lon=lon, lat=lat, distdf = coastdistdat))/1000) %>% 
  ungroup() %>% 
  group_by(est_year) %>% 
  mutate(
    assemblage.dist95 = quantile(coastdist_km, 0.95),
    assemblage.lat95 = quantile(lat, 0.95),
  ) %>% 
  ungroup() %>% 
  group_by(est_year, comname) %>% 
  mutate(
    spp.lat95 = quantile(lat, 0.95), 
    spp.lat95round = round((spp.lat95 + 0.25) * 2) / 2 - 0.25,
    spp.dist90 = quantile(coastdist_km, 0.90), 
    spp.dist95 = quantile(coastdist_km, 0.95), 
    spp.dist95.wt = Hmisc::wtd.quantile(coastdist_km, weights=biomass_kg, probs=0.95, normwt=FALSE),
    spp.dist99 = quantile(coastdist_km, 0.99), 
    spp.distmax = max(coastdist_km),
    depth.mean = mean(avgdepth),
    depth.mean.wt = weighted.mean(avgdepth, w=biomass_kg)) %>%
  ungroup() %>% 
  dplyr::select(est_year, comname, numyears, numobs, numobsyear, meanobsyear, assemblage.dist95, spp.dist90, spp.dist95, spp.dist99, spp.distmax, biomass_kg, depth.mean, depth.mean.wt, assemblage.lat95, spp.lat95, spp.lat95round, spp.dist95.wt) %>% 
  distinct() 


# equatorial 
eqdat.stats <- eqdat %>% 
  rename("lat" = "decdeg_beglat",
         "lon" = "decdeg_beglon") %>%
  rowwise() %>% 
  mutate(coastdist_km = (get_length(lon=lon, lat=lat, distdf = coastdistdat))/1000) %>% 
  ungroup() %>% 
  # group_by(est_year) %>% 
  # mutate(
  #   assemblage.dist05 = quantile(coastdist_km, 0.05),
  #   assemblage.lat05 = quantile(lat, 0.05),
  # ) %>% 
  # ungroup() %>% 
  group_by(est_year, comname) %>% 
  mutate(
    spp.lat05 = quantile(lat, 0.05), 
    spp.lat05round = round((spp.lat05 + 0.25) * 2) / 2 - 0.25,
    spp.dist05 = quantile(coastdist_km, 0.05), 
    # spp.dist95 = quantile(coastdist_km, 0.95), 
    spp.dist05.wt = Hmisc::wtd.quantile(coastdist_km, weights=biomass_kg, probs=0.05, normwt=FALSE),
    spp.dist01 = quantile(coastdist_km, 0.01), 
    spp.distmax = max(coastdist_km),
    depth.mean = mean(avgdepth),
    depth.mean.wt = weighted.mean(avgdepth, w=biomass_kg)) %>%
    ungroup() %>% 
    #dplyr::select(est_year, comname, numyears, numobs, numobsyear, meanobsyear, assemblage.dist95, spp.dist90, spp.dist95, spp.dist99, spp.distmax, biomass_kg, depth.mean, depth.mean.wt, assemblage.lat95, spp.lat95, spp.lat95round, spp.dist95.wt) %>% 
    distinct() 
    
write_rds(eqdat.stats, here("processed-data","eqdat.stats.rds"))
write_rds(poldat.stats, here("processed-data","poldat.stats.rds"))

rm(list=ls())
