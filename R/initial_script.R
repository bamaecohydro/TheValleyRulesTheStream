#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Fuck around and find out
#Coder: Nate et al
#Date Created: 3/1/2023
#Purpose: Develop initial function to characterize valley characteristics
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Shared notes doc: https://docs.google.com/document/d/1OrxS69xVh48elGpXuBDVH0tdwRh-sJSfpwKd-3UEsqg/edit#
#Valley Delineation Method: https://doi.org/10.5194/esurf-5-369-2017

#Steps
#   Download NED Data
#   Define River Corridor
#   Define Flow Network
#   Identify study reach
#   Develop Geomorphic Metrics 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.0 Setup workspace ----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#clear workspace (sorry JP!)
remove(list=ls())

#load packages
library(tidyverse)
library(raster)
library(sf)
library(nhdplusTools)
library(elevatr)
library(whitebox)
library(mapview)

#Create temp dir
temp_dir <- tempdir()

#read in gagesII shapefile
gage <- get_gagesII(id = "06879650")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.0 Download NED data --------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Download NED
r <- get_elev_raster(gage, z=14) 

#plot
mapview(r) + mapview(gage)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.0 Define flow net ----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4.0 Define valley bottom -----------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.1 Create hand raster --------------------------------------------------------
#Write raster to temp_dir
writeRaster(r,paste0(temp_dir, '\\r.tif'), overwrite=T)

#Smooth the dem
wbt_fast_almost_gaussian_filter(
  input = "r.tif",
  output = "r_smooth.tif", 
  sigma = 1.8, 
  wd = temp_dir
)

#breach depressions (for flownet analysis only)
wbt_breach_depressions(
  dem = "r_smooth.tif",
  output = "r_breach.tif", 
  wd= temp_dir
)

#flow accumulation
wbt_d8_flow_accumulation(
  input = 'r_breach.tif',
  output = 'r_fac.tif',
  pntr = F,
  wd = temp_dir
)

#Create stream raster
wbt_extract_streams(
  flow_accum = "r_fac.tif", 
  output = "streams.tif", 
  threshold = 1000,
  wd=temp_dir
)

#Estiamte Height above stream
wbt_elevation_above_stream(
  dem = "r_breach.tif", 
  stream = "streams.tif", 
  output = "hand.tif",
  wd = temp_dir
)

#Pull relative releif (hand?) and channel rasters in to R environment
hand <- raster(paste0(temp_dir, "\\hand.tif"))
channel <- raster(paste0(temp_dir, "\\streams.tif"))

#4.2 Identify valley bottom ----------------------------------------------------
#Isolate hand values
hand_values <- values(hand) %>% as_vector() %>% na.omit()

#limit to 1st and 3rd quartiles
hand_q25 <- quantile(hand_values, probs = 0.25)
hand_q75 <- quantile(hand_values, probs = 0.75)
hand_subset <- hand_values[hand_values>hand_q25]
hand_subset <- hand_values[hand_values<hand_q75]

#Estiamte distributional parameters
mean <- mean(hand_subset)
sd   <- sd(hand_subset)


#Create of vector of hand values (if conformed to normal dist)
hand_expected <- 
  rnorm(
    n=length(hand_values), 
    mean = mean(hand_values), 
    sd = sd(hand_values))

#Create tibble to look at deviation from qqline
threshold<-tibble(
  theoretical_quantiles = qnorm(seq(0, 1, by = 0.001)),
  emperical_quantiles   = quantile(hand_values, probs = seq(0,1, by = 0.001)),
  normal_quant          = qnorm(mean, sd, p = seq(0, 1, by = 0.001)), 
  diff = (emperical_quantiles - normal_quant)/normal_quant*100) %>% 
  #Filter to where values are within 1% of dist
  dplyr::filter(abs(diff)<1) %>% 
  #Identify threshold
  slice(1) %>% dplyr::select(emperical_quantiles ) %>% pull()

#Create binary raster
valley <- hand < threshold
valley[valley==0] <- NA

#Convert to polygon using wbt
writeRaster(valley, paste0(temp_dir, "\\valley.tif"), overwrite=T)
wbt_raster_to_vector_polygons(
  input = "valley.tif", 
  output = "output.shp", 
  wd=temp_dir)
valley_shp <- st_read(paste0(temp_dir,"\\output.shp")) 
st_crs(valley_shp) = 4269

#Plot for funzies
mapview(gage)+mapview(r)+mapview(valley_shp)

#Isolate NHDplus Reach ---------------------------------------------------------


#Estimate valley characteristics ------------------------------------------------
