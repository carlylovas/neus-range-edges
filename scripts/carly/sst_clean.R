## Cleaning temp data
library(here) 
library(tidyverse)
library(raster)
library(sf)
library(oceanmap)
library(data.table)
source(here("functions","sfc_as_cols.R"))

# create shelf mask
lonrange <- c(-77, -66)
latrange <- c(35, 45)

bathy <- get.bathy(lon = lonrange, lat = latrange, visualize = F, res = 15) 

bathy.crs <- bathy %>% 
  as("SpatialPolygonsDataFrame") %>% 
  st_as_sf() %>% 
  st_crs()

eezs <- st_read(here("data/World_EEZ_v12_20231025","eez_v12.shp")) 

useez <- eezs %>% 
  dplyr::filter(SOVEREIGN1 == "United States") %>% 
  st_transform(crs=bathy.crs) # reproject to match bathymetry 
#slow

#also slow
bathy.mask <- bathy %>% 
  as("SpatialPolygonsDataFrame") %>% 
  st_as_sf() %>% # retains CRS of bathy.raw
  dplyr::filter(layer <= 300) %>% # get rid of values over 300m deep
  st_intersection(st_buffer(st_union(useez), 0)) # keep only points within the EEZ; crop out lakes, Canada; st_buffer() prevents error from sf about points not matching exactly

# will need to pull SODA from Box
box_path <- "/Users/clovas/Library/CloudStorage/Box-Box/RES_Data/SODA/"
soda <- paste0(box_path, "SODA_Temp_Red.nc")

soda.neus <- raster::stack(soda)

soda.neus.crop <- soda.neus %>%
   raster::mask(bathy.mask) %>% 
   raster::crop(extent(bathy.mask)) 

unstack(soda.neus.crop) 
unstack(soda.neus)

soda.neus <- soda.neus.crop %>%
  raster::as.data.frame(xy=TRUE, long=TRUE)  %>%
  mutate(
  time = lubridate::ymd(str_sub(str_replace(layer, "X", ""), 1, 10)),
  year = lubridate::year(time),
  month = lubridate::month(time)) %>%
  dplyr::select(-layer) %>%
  rename(surftemp = value) %>%
  filter(!is.na(surftemp)) %>%
  mutate(year_measured = ifelse(month %in% c(1,2), year-1, year),
         year_match = year_measured + 1) %>%
  filter(year_match > min(year_match)) # ditch year_match=1980; it is misleading because it is actually only the first 3 months of the previous year

# load raw temperature datasets, crop to extent of shelf, change back into dataframes 
hadisst.neus <- read_rds(here("processed-data","hadisst_raw.rds")) %>% 
  filter(!is.na(sst)) %>% # do early to make the object smaller
  st_as_sf(coords=c("longitude","latitude"), crs = bathy.crs) %>% 
  st_join(bathy.mask, left=FALSE) %>% 
  rename(bathymetry = layer) %>% 
  sfc_as_cols() %>% 
  as.data.table() %>% 
  mutate(
    time = lubridate::ymd(str_sub(time, 1, 10)),
    year = lubridate::year(time),
    month = lubridate::month(time)
  ) %>% 
  mutate(year_measured = ifelse(month %in% c(1,2), year-1, year),
         year_match = year_measured + 1) %>% 
  filter(year_match > min(year_match))


# the trawl survey data is from spring (starts in april), and the coldest temperatures are often found in early spring
# so it's not appropriate to compare a survey from spring in one year to either all temperatures in that calendar year or last year
# rather, we're "redefining" the year to run from march to february, and then comparing edge position to the previous 12 months of temperaturue data 

oisst.neus <- read_rds(here("processed-data","oisst_raw.rds")) %>% 
  filter(!is.na(sst)) %>% # do early to make the object smaller
  st_as_sf(coords=c("longitude","latitude"), crs = bathy.crs) %>% 
  st_join(bathy.mask, left=FALSE) %>% 
  rename(bathymetry = layer) %>% 
  sfc_as_cols() %>% 
  as.data.table() %>% 
  mutate(
    time = lubridate::ymd(str_sub(time, 1, 10)),
    year = lubridate::year(time),
    month = lubridate::month(time)
  ) %>% 
  mutate(year_measured = ifelse(month %in% c(1,2), year-1, year),
         year_match = year_measured + 1) %>% 
  filter(year_match > min(year_match))

# save files--they are slow to generate!
write_rds(oisst.neus, here("processed-data","oisst_neus.rds"))
write_rds(hadisst.neus, here("processed-data","hadisst_neus.rds"))
write_rds(soda.neus, here("processed-data","soda_neus.rds"))


oisst.stats <- oisst.neus %>% 
  dplyr::select(-zlev, -bathymetry) %>% 
  group_by(year, month, x, y) %>% 
  mutate(
    cell.month.mean = mean(sst)
  ) %>% 
  ungroup() %>% 
  dplyr::select(year, month, year_measured, year_match, x, y, cell.month.mean) %>% # trim down to resolution of other datasets 
  distinct() %>% 
  group_by(year_measured) %>% 
  mutate(
    year.month.mean = mean(cell.month.mean),
    year.month.sd = sd(cell.month.mean),
    year.month.max = max(cell.month.mean),
    year.month.min = min(cell.month.mean)
  ) %>% 
  ungroup() %>% 
  group_by(year_measured, y) %>% 
  mutate(
    lat.year.month.mean = mean(cell.month.mean),
    lat.year.month.sd = sd(cell.month.mean),
    lat.year.month.max = max(cell.month.mean),
    lat.year.month.min = min(cell.month.mean) 
  ) %>% 
  ungroup() %>% 
  dplyr::select(y, year_measured, year_match, year.month.mean, year.month.max, year.month.sd, year.month.min, lat.year.month.mean, lat.year.month.sd, lat.year.month.max, lat.year.month.min) %>% 
  distinct()

hadisst.stats <- hadisst.neus %>% 
  group_by(year_measured) %>% 
  mutate(
    year.month.mean = mean(sst),
    year.month.sd = sd(sst),
    year.month.max = max(sst),
    year.month.min = min(sst)
  ) %>% 
  ungroup() %>% 
  group_by(year_measured, y) %>% 
  mutate(
    lat.year.month.mean = mean(sst),
    lat.year.month.sd = sd(sst),
    lat.year.month.max = max(sst),
    lat.year.month.min = min(sst) 
  ) %>% 
  ungroup() %>% 
  dplyr::select(y, year_measured, year_match, year.month.mean, year.month.max, year.month.sd, year.month.min, lat.year.month.mean, lat.year.month.sd, lat.year.month.max, lat.year.month.min) %>% 
  distinct()

soda.stats <- soda.neus %>% 
  group_by(year_measured) %>% 
  mutate(
    year.month.mean = mean(surftemp),
    year.month.sd = sd(surftemp),
    year.month.max = max(surftemp),
    year.month.min = min(surftemp)
  ) %>% 
  ungroup() %>% 
  group_by(year_measured, y) %>% 
  mutate(
    lat.year.month.mean = mean(surftemp),
    lat.year.month.sd = sd(surftemp),
    lat.year.month.max = max(surftemp),
    lat.year.month.min = min(surftemp) 
  ) %>% 
  ungroup() %>% 
  dplyr::select(y, year_measured, year_match, year.month.mean, year.month.max, year.month.sd, year.month.min, lat.year.month.mean, lat.year.month.sd, lat.year.month.max, lat.year.month.min) %>% 
  distinct()

write_rds(oisst.stats, here("processed-data","oisst_stats.rds"))
write_rds(hadisst.stats, here("processed-data","hadisst_stats.rds"))
write_rds(soda.stats, here("processed-data","soda_stats.rds"))
rm(list=ls())
