## Combining A.Fredston's range edges and L.Carlson's latitudinal accumulation code
## load libraries 
library(here)
library(tidyverse)
library(gmRi)
install.packages("matrixStats")
library(matrixStats)
install.packages("ggridges")
library(ggridges)

# species data
species <- read_csv(here("data", "speciesList_inNECC.csv")) %>%
  rename("comname" = "Species_comnam") %>%
  mutate(comname   = tolower(comname)) %>%
  select(comname)

# Load NEFSC Bottom Trawl Survey data ####
trawl_data <- gmri_survdat_prep(
  survdat_source ="most recent",
  box_location ="cloudstorage")

clean_survey  <- trawl_data %>% 
  distinct(svspp, est_year, survey_area, stratum, tow, id, est_towdate, season, comname, catchsex, .keep_all = T) %>%
  group_by(svspp, est_year, survey_area, stratum, tow, id, est_towdate, season, 
           avgdepth, surftemp, bottemp, decdeg_beglat, decdeg_beglon, comname, abundance) %>% 
  filter(comname %in% species$comname) %>% 
  summarise(biomass_kg = sum(biomass_kg, na.rm = T), .groups = "drop")

# Sum total biomass by species in each year (L. Carlson)
sum_biomass_year <- trawl_data %>%
  group_by(est_year, svspp) %>% 
  summarise(annual_species_biomass = sum(biomass_kg))

# Total catch by unique haul (L. Carlson)
sum_biomass_haul <- trawl_data %>% 
  distinct(id, season, comname, catchsex, biomass_kg, abundance, .keep_all = T) %>%
  group_by(svspp, id) %>%
  summarise(species_biomass = sum(biomass_kg), species_abundance = sum(abundance))

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

weighted_survey_data <- grouped_center_bio(clean_survey, est_year) # leaving out season 


weighted_survey_data %>%
  filter(comname == "acadian redfish") %>% 
  mutate(decade = 10*est_year %/% 10) %>%
  ggplot() +
  geom_density_ridges(aes(x = avg_lat, y = decade, fill = as.factor(decade)), alpha = .9) +
  guides(fill = guide_legend(title = "Decade")) +
  # scale_y_continuous(breaks = c("1970", "2010")) +
  s# cale_y_reverse() +
  # ylim(c(2010, 1960)) +
  scale_fill_gmri() +
  theme_gmri(legend.position = "left") 

# Latitudinal distribution  
lat_dist <- weighted_survey_data %>% 
  mutate(decade = 10*est_year %/% 10) %>%
  group_by(comname) %>% 
  nest() %>% 
  mutate(plot = map2(data, comname, function(x,y){
    out <- ggplot(x) +
      geom_density_ridges(aes(x = avg_lat, y = decade, fill = as.factor(decade)), alpha = .9) +
      guides(fill = guide_legend(title = "Decade")) +
      ylab("Decade") + xlab("Biomass-weighted average latitude") + ggtitle(str_to_sentence(comname)) +
      ylim(c(1970, NA)) +
      scale_fill_gmri() +
      coord_flip() +
      theme_gmri()
  }))

lat_dist$plot[9] # I like this

# Center of Biomass 
centerBiomass <- weighted_survey_data %>% 
  select(comname, est_year, avg_lat, avg_lon, avg_biomass) %>% 
  group_by(comname) %>% 
  nest() %>% 
  mutate(lat_lm   = map(data, ~lm(avg_lat ~ est_year, data = .x)),
         lon_lm   = map(data, ~lm(avg_lon ~ est_year, data = .x)),
         tidy_lat = map(lat_lm, broom::tidy),
         tidy_lon = map(lon_lm, broom::tidy))

lat_centroid <- centerBiomass %>% 
  select(comname, lat_lm, tidy_lat) %>% 
  # mutate(slope_lat = tidy_lat %>% map_dbl(function(x) x$estimate[2]))
  unnest(tidy_lat) %>%
  filter(!term == "(Intercept)") %>% 
  select(comname, lat_lm, estimate, p.value) %>%
  rename("slope" = "estimate") %>%
  drop_na() %>%
  mutate(movement = ifelse(slope > 0, "T","F")) %>% 
  mutate(significant = ifelse(p.value < 0.05, "T","F")) %>% 
  arrange(desc(movement), desc(significant), comname)


lat_centroid$trend = NA
lat_centroid$trend[lat_centroid$movement == "T" & lat_centroid$significant == "T"] = "Northward" 
lat_centroid$trend[lat_centroid$movement == "T" & lat_centroid$significant == "F"] = "Stable"
lat_centroid$trend[lat_centroid$movement == "F" & lat_centroid$significant == "F"] = "Stable"
lat_centroid$trend[lat_centroid$movement == "F" & lat_centroid$significant == "T"] = "Southward"

# Calculate and plot 5%, 10%, 25%, 75%, 90%, and 95% biomass-weighed percentiles
grouped_quantiles <- function(clean_survey, ...){
  clean_survey %>% 
    group_by(comname, ...) %>% 
    summarise(
      # Un-weighted averages
      total_biomass   = sum(biomass_kg),
      avg_biomass     = mean(biomass_kg),
      avg_lat         = mean(decdeg_beglat),
      # Weight quantiles
      `5%`  = Hmisc::wtd.quantile(decdeg_beglat, weights = biomass_kg, probs = 0.05, na.rm = T),
      `10%` = Hmisc::wtd.quantile(decdeg_beglat, weights = biomass_kg, probs = 0.10, na.rm = T), 
      `25%` = Hmisc::wtd.quantile(decdeg_beglat, weights = biomass_kg, probs = 0.25, na.rm = T),
      `50%` = Hmisc::wtd.quantile(decdeg_beglat, weights = biomass_kg, probs = 0.50, na.rm = T),
      `75%` = Hmisc::wtd.quantile(decdeg_beglat, weights = biomass_kg, probs = 0.75, na.rm = T), 
      `90%` = Hmisc::wtd.quantile(decdeg_beglat, weights = biomass_kg, probs = 0.90, na.rm = T),
      `95%` = Hmisc::wtd.quantile(decdeg_beglat, weights = biomass_kg, probs = 0.95, na.rm = T),
      .groups = "drop") %>%
    mutate(across(where(is.numeric), round, 4))
}

quantiles <- grouped_quantiles(clean_survey, est_year)

quantiles %>% 
  filter(comname == "black sea bass") %>% 
  mutate(decade = 10*est_year %/% 10) %>%
  ggplot() + 
  geom_col(aes(x = avg_biomass, y = as.factor(`5%`)), fill = "#00608A") +
  geom_col(aes(x = avg_biomass , y = as.factor(`25%`)), fill = "#EACA00") + 
  geom_col(aes(x = avg_biomass , y = as.factor(`50%`)), fill = "#535353") +
  geom_col(aes(x = avg_biomass, y = as.factor(`75%`)), fill = "#EACA00") +
  geom_col(aes(x = avg_biomass, y = as.factor(`95%`)), fill = "#00608A") +
   scale_y_discrete(breaks= c(36, 37, 38, 39, 40)) +
  facet_wrap(~decade, scales = "free_y") +
  theme_gmri() # I don't like this, refer back to Lindsey's code
  

# Latitudinal shift percentiles 
lat_shift <- quantiles %>% 
  pivot_longer(cols = 6:12, names_to = "quantile", values_to = "lat") %>%
  select(comname, est_year, quantile, lat) %>%
  group_by(comname, quantile) %>%
  mutate(rollmean = zoo::rollapplyr(lat, width = 5, FUN = mean, align = "center", partial = T)) %>%
  group_by(comname) %>% 
  nest() %>% 
  mutate(plot = map2(data, comname, function(x,y){
    out <- ggplot(x) +
      geom_line(aes(x = est_year, y = rollmean, color = quantile), linetype = 2) +
      geom_smooth(aes(x = est_year, y = rollmean, color = quantile), method = "lm", se = F) +
      ylab("Latitude") + xlab("Year") + ggtitle(str_to_sentence(comname)) + 
      scale_color_gmri() +
      theme_gmri()
    }))

lat_shift$plot[[1]]

# 95th Percentile
quantiles %>% 
  filter(comname == "black sea bass") %>%
  mutate(decade = 10*est_year %/% 10) %>% 
  ggplot() +
  geom_density_ridges(aes(x = `95%`, y = decade, fill = as.factor(decade)), alpha = .9) +
  guides(fill = guide_legend(title = "Decade")) +
  ylab("Decade") + xlab("95th Percentile") + ggtitle("Black sea bass") +
  ylim(c(1970, NA)) +
  scale_fill_gmri() +
  coord_flip() +
  theme_gmri()
