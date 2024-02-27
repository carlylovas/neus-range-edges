## load libraries 
library(here)
library(tidyverse)
library(gmRi)
install.packages("matrixStats")
library(matrixStats)

# Load NEFSC Bottom Trawl Survey data ####
clean_survey <- gmri_survdat_prep(
  survdat_source ="most recent",
  box_location ="cloudstorage")

clean_survey  <- clean_survey %>% 
  mutate(haulid = paste(est_year, cruise6, stratum, station, sep="-")) %>%
  distinct(est_year, survey_area, stratum, tow, haulid, est_towdate, season, comname, catchsex, .keep_all = T) %>%
  group_by(est_year, survey_area, stratum, tow, haulid, est_towdate, season, 
           avgdepth, surftemp, bottemp, decdeg_beglat, decdeg_beglon, comname, abundance) %>% 
  summarise(biomass_kg = sum(biomass_kg, na.rm = T), .groups = "drop")

# Weight by biomass
grouped_center_bio <- function(clean_survey, ...){
  clean_survey %>% 
    group_by(comname, ...) %>% 
    summarise(
      # Un-weighted averages
      total_biomass   = sum(biomass_kg),
      avg_biomass     = mean(biomass_kg),
      biomass_sd      = sd(biomass_kg),
      # All below are weighted by biomass
      avg_depth       = weightedMean(avgdepth, w = biomass_kg, na.rm = T),
      avg_bot_temp    = weightedMean(bottemp, w = biomass_kg, na.rm = T),
      avg_sur_temp    = weightedMean(surftemp, w = biomass_kg, na.rm = T),
      avg_lat         = weightedMean(decdeg_beglat, w = biomass_kg, na.rm = T),
      avg_lon         = weightedMean(decdeg_beglon, w = biomass_kg, na.rm = T),
      depth_sd        = weightedSd(avgdepth, w = biomass_kg, na.rm = T),
      temp_sd         = weightedSd(bottemp, w = biomass_kg, na.rm = T),
      lat_sd          = weightedSd(decdeg_beglat, w = biomass_kg, na.rm = T),
      lon_sd          = weightedSd(decdeg_beglon, w = biomass_kg, na.rm = T),
      .groups = "drop") 
}

weighted_survey_data <- grouped_center_bio(clean_survey, est_year, season) %>%
  mutate(decade = 10*est_year %/% 10) # not sure we need this yet, this technically is center of biomass?

## Range edges ####
# data preferences
num_obs_year_cutoff <- 10 # how many times does a species need to be observed in a year for that species*year combo to be included in analysis? 
numyears_cutoff <- 10 
northerncutoff <- 42 # max lat at which a range edge can start and still be classified a northern edge (beyond this, the species could extend into Cape Cod / Canada )
southerncutoff <- 32 # min lat

# poleward edge dataframe
poldat.firstyear <- clean_survey %>% # calculate first year when species was seen
  group_by(comname) %>%
  mutate(firstyear = min(est_year)) %>% 
  filter(est_year == firstyear) %>% 
  ungroup() %>% 
  group_by(comname, est_year) %>% 
  mutate(firstlat_max = max(decdeg_beglat)) %>% 
  ungroup() %>% 
  dplyr::select(comname, firstyear, firstlat_max) %>%
  distinct() 

poldat <- clean_survey %>% 
  left_join(poldat.firstyear, by="comname") %>% 
  filter(firstlat_max < 42) %>% 
  group_by(comname) %>% 
  mutate(numobsyear = length(unique(haulid))) %>% 
  ungroup() %>% 
  filter(numobsyear >= num_obs_year_cutoff) %>% 
  group_by(comname) %>% 
  mutate(numobs = length(unique(haulid)),
         numyears = length(unique(est_year)),
         meanobsyear = numobs / numyears) %>% 
  ungroup() %>% 
  filter(numyears >= numyears_cutoff)

poldat.summ <- poldat %>% 
  group_by(comname, numobs, numyears, meanobsyear) %>% 
  summarise()

# equatorial edge (following poleward code since I don't have access to aquamaps)
eqdat.firstyear <- clean_survey %>% # calculate first year when species was seen
  group_by(comname) %>%
  mutate(firstyear = min(est_year)) %>% 
  filter(est_year == firstyear) %>% 
  ungroup() %>% 
  group_by(comname, est_year) %>% 
  mutate(firstlat_min = min(decdeg_beglat)) %>% 
  ungroup() %>% 
  dplyr::select(comname, firstyear, firstlat_min) %>%
  distinct() 

eqdat <- clean_survey %>% 
  left_join(eqdat.firstyear, by="comname") %>% 
  filter(firstlat_min > 32) %>% # not sure about this 
  group_by(comname) %>% 
  mutate(numobsyear = length(unique(haulid))) %>% 
  ungroup() %>% 
  filter(numobsyear >= num_obs_year_cutoff) %>% 
  group_by(comname) %>% 
  mutate(numobs = length(unique(haulid)),
         numyears = length(unique(est_year)),
         meanobsyear = numobs / numyears) %>% 
  ungroup() %>% 
  filter(numyears >= numyears_cutoff)

eqdat.summ <- eqdat %>% 
  group_by(comname, numobs, numyears, meanobsyear) %>% 
  summarise()


# final check
intersect(poldat.summ$comname, eqdat.summ$comname)

# save out 
saveRDS(eqdat, "processed-data/eqdat.rds")
saveRDS(poldat, "processed-data/poldat.rds")

write_csv(eqdat.summ, "processed-data/eqdat_summary.csv")
write_csv(poldat.summ, "processed-data/poldat_summary.csv")

rm(list=ls())
