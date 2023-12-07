#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Topographic controls on stream expansion and contraction
#Coder: Nate Jones (cnjones7@ua.edu)
#Date Created: 12/6/2023
#Purpose: Recreate Prancevich and Kirchner (https://doi.org/10.1029/2018GL081799)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setup workspace --------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#clear workspace (sorry JP!)
remove(list=ls())

#load packages
library(tidyverse)
library(raster)
library(sf)
library(elevatr)
library(whitebox)
library(stars)
library(mapview)
library(tmap)

#Create temp dir
temp_dir <- "C:\\WorkspaceR\\TheValleyRulesTheStream\\scratch\\"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Gather data ------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Define watershed outlet
outlet<-tibble(
  lat = 33.76218, 
  lon = -85.5955) %>% 
  st_as_sf(
    coords = c("lon", "lat"), 
    crs = '+proj=longlat +datum=WGS84 +no_defs') %>% 
  st_transform(crs = 3160)

sf_use_s2(FALSE)

#Download DEM
dem <- get_elev_raster(outlet, z=14)

#Export data to temp directory
writeRaster(dem,paste0(temp_dir, 'dem.tif'), overwrite=T)
st_write(outlet, paste0(temp_dir, "outlet.shp"), append=T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create fdr and face ----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Smooth the dem
wbt_fast_almost_gaussian_filter(
  input = "dem.tif",
  output = "dem_smooth.tif", 
  sigma = 1.8, 
  wd = temp_dir
)

#breach depressions 
wbt_breach_depressions(
  dem = "dem_smooth.tif",
  output = "dem_breach.tif", 
  wd= temp_dir
)

#flow direction
wbt_d8_pointer(
  dem = "dem_breach.tif",
  output = "fdr.tif",
  wd = temp_dir)

#flow accumulation
wbt_d8_flow_accumulation(
  input = 'dem_breach.tif',
  output = 'fac.tif',
  pntr = F,
  wd = temp_dir
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delineate Watershed-----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Snap pour point
wbt_snap_pour_points(
  pour_pts = "outlet.shp",
  flow_accum = "fac.tif",
  snap_dist = 100,
  output = "snap.shp",
  wd = temp_dir
)

#Check snap pour point with mapview (you may need to change the search dist)
fac <- raster(paste0(temp_dir, "//fac.tif"))
fdr <- raster(paste0(temp_dir, "//fdr.tif"))
snap<- st_read(paste0(temp_dir,"//snap.shp"))

#Create watershed
wbt_watershed(
  d8_pntr  = "fdr.tif",
  pour_pts = "snap.shp",
  output   = "watershed.tif",
  wd       = temp_dir
)

#read into R
watershed <- raster(paste0(temp_dir, "//watershed.tif"))

#Convert raster to vector
watershed <- watershed %>% st_as_stars() %>% st_as_sf(., merge = TRUE)

#Crop fac and fdr
fac <- mask(fac, watershed)
fac <- crop(fac, watershed)
fdr <- mask(fdr, watershed)
fdr <- crop(fdr, watershed)

#write to workspace
writeRaster(fac, paste0(temp_dir,"fac_shed.tif"), overwrite=TRUE)
writeRaster(fdr, paste0(temp_dir,"fdr_shed.tif"), overwrite=TRUE)

#plot for funzies
mapview(log10(fac))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Estimate scaling factors  ----------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Estimate alpha based on accumulation area threshold for channel heads ---------
#Create function
fun <- function(n_cells){
  #Create flow net
  wbt_extract_streams(
    flow_accum = "fac_shed.tif",
    output = "flow_net.tif",
    threshold = n_cells, 
    wd = temp_dir)
  
  #convert to vector
  wbt_raster_streams_to_vector(
    streams = "flow_net.tif", 
    d8_pntr = "fdr_shed.tif", 
    output = "streams.shp", 
    wd = temp_dir)
  
  #Read into R environement
  flow_net<-st_read(
    paste0(temp_dir, "\\streams.shp"), 
    crs = st_crs(dem))
  
  #Quantify flownet length
  stream_length <- st_length(flow_net) %>% sum(., na.rm=T)
  
  #Quantify flow density
  stream_density <- stream_length/st_area(watershed)
  
  #Create export 
  tibble(
    A_channel_head_m2 = n_cells*res(fac)[1]*res(fac)[2], 
    stream_length_m = stream_length %>% paste0() %>% as.numeric(),
    stream_density_m_m2 = stream_density %>% paste0() %>% as.numeric()) 
}

#Apply function
output <- lapply(
  X = seq(10, 7500, 10), 
  FUN = fun
) %>% bind_rows()

#ugly plot (for now)
plot(output$A_channel_head_m2, output$stream_density_m_m2, log="xy")
