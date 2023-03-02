#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Fuck around and find out
#Coder: Nate et al
#Date Created: 3/1/2023
#Purpose: Develop initial function to characterize valley characteristics
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Shared notes doc: https://docs.google.com/document/d/1OrxS69xVh48elGpXuBDVH0tdwRh-sJSfpwKd-3UEsqg/edit#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setup workspace --------------------------------------------------------------
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

#readin gagesII shapefile
gages <- read_sf("data/gagesII/gagesII_9322_sept30_2011.shp")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Download NED data ------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Konza gage for testing
gage <- gages %>% filter(STAID == "06879650")

#Download NED
r <- get_elev_raster(gage, z=14) 

#plot
mapview(r) + mapview(gage)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define valley bottom and reach reach -----------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Use delineation method here: https://doi.org/10.5194/esurf-5-369-2017

#Create hand raster ------------------------------------------------------------
#Create temp dir
temp_dir <- tempdir()

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


#Identify valley bottom --------------------------------------------------------
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
  slice(1) %>% select(emperical_quantiles ) %>% pull()

#Create bunary raster
valley <- hand < threshold
valley[valley==0] <- NA

#Plot for funzies
mapview(gage)+mapview(r)+mapview(valley)

#Isolate Reach -----------------------------------------------------------------

#Estimate valley characteristics ------------------------------------------------
