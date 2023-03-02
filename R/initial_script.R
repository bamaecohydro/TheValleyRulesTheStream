#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Fuck around and find out
#Coder: Nate et al
#Date Created: 3/1/2023
#Purpose: Develop initial function to characterize valley characteristics
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Shared notes doc: https://docs.google.com/document/d/1OrxS69xVh48elGpXuBDVH0tdwRh-sJSfpwKd-3UEsqg/edit#


#Steps to work on
# 1) dowlnoad gagesII
# 2) downlaod ned raster
# 3) snag valley charactersitics


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


hand<-raster(paste0(temp_dir, "\\hand.tif"))
channel <- raster(paste0(temp_dir, "\\streams.tif"))
qqnorm(values(hand))







