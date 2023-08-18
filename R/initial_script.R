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

#Need to figure out soils data download + confining layer depth
#  perhaps https://cran.r-project.org/web/packages/soilDB/soilDB.pdf

#Watershed storage data:
# (1) Watershed Shapefiles Obtained from John Hammond (during the DryRiversRCN)
# (2) Depth to bedrock: https://doi.org/10.1371/journal.pone.0169748
#web address: https://data.isric.org/geonetwork/srv/eng/catalog.search#/metadata/f36117ea-9be5-4afd-bb7d-7a3e77bf392a
# (3) porosity:  https://doi.org/10.1002/2014GL059856
#web address: https://dataverse.scholarsportal.info/dataset.xhtml?persistentId=doi:10.5683/SP2/DLGXYO

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
plot(r) + plot(gage %>% st_geometry(), add=T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.0 Define flow net ----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create function to create stream layer ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stream_fun<-function(r, threshold_n_cells, temp_dir){
  
  #Load libraries of interest
  library(tidyverse) #join the cult!
  library(raster)
  library(sf)
  library(whitebox)

  #Export DEM to scratch workspace
  writeRaster(r, paste0(temp_dir,"\\r.tif"), overwrite=T)
  
  #Smooth DEM
  wbt_fast_almost_gaussian_filter(
    input = "r.tif",
    output = "r_smooth.tif", 
    sigma = 1.8, 
    wd = temp_dir
  )
  
  #breach depressions
  wbt_breach_depressions(
    dem =    "r_smooth.tif",
    output = "r_breached.tif",
    fill_pits = F,
    wd = temp_dir)
  
  #Flow direction raster
  wbt_d8_pointer(
    dem= "r_breached.tif",
    output ="fdr.tif",
    wd = temp_dir
  )
  
  #Flow accumulation raster
  wbt_d8_flow_accumulation(
    input = "r_breached.tif",
    output = "fac.tif",
    wd = temp_dir
  )
  
  #Create Stream Layer
  stream_grd<-raster(paste0(temp_dir,"\\fac.tif"))
  stream_grd[stream_grd<threshold_n_cells]<-NA
  writeRaster(stream_grd, paste0(temp_dir,"\\stream.tif"), overwrite=T)
  
  #Identify links
  wbt_stream_link_identifier(
    d8_pntr = 'fdr.tif',
    streams = 'stream.tif',
    output = 'stream_link.tif',
    wd = temp_dir
  )
  
  #Convert stream to vector
  wbt_raster_streams_to_vector(
    streams = "stream_link.tif",
    d8_pntr = "fdr.tif",
    output = "streams.shp",
    wd = temp_dir)
  
  #Read streams layer in 
  stream_shp<-st_read(paste0(temp_dir,"\\streams.shp"), crs=st_crs(r))
  
  #Export streams
  list(stream_grd, stream_shp)
}

#Apply streams function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Run function
streams<-stream_fun(r, threshold_n_cells=10000, temp_dir)

#Define components of output list
stream_grd <- streams[[1]]
stream_shp <- streams[[2]]

#Plot for funzies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mapview(stream_shp) + mapview(gage)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4.0 Define valley bottom -----------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4.1 Create function to define valley ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
valley_fun <- function(r, threshold_n_cells, temp_dir){

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
    threshold = threshold_n_cells,
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
  valley_grd <- hand < threshold
  valley_grd[valley_grd==0] <- NA
  
  #Convert to polygon using wbt
  writeRaster(valley_grd, paste0(temp_dir, "\\valley.tif"), overwrite=T)
  wbt_raster_to_vector_polygons(
    input = "valley.tif", 
    output = "output.shp", 
    wd=temp_dir)
  valley_shp <- st_read(paste0(temp_dir,"\\output.shp"), crs=st_crs(r)) 
  
  #Export Valley Bottom
  list(valley_grd, valley_shp)
}

#Execute function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Run function
valleys <- valley_fun(r, 10000, temp_dir)

#Define components of output list
valley_grd <- valleys[[1]]
valley_shp <- valleys[[2]]

#Plot for funzies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mapview(stream_shp) + mapview(valley_shp) + mapview(gage)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5.0 Reproject spatial data into planar coordinates ---------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Define utm zone
utm_zone <- floor((st_coordinates(gage)[1] + 180)/6)+1

#Define CRS
crs <- paste0("+proj=utm +zone=",utm_zone," +datum=WGS84 +units=m +no_defs")

#reproject raster datasets of interest
r          <- projectRaster(r,          crs = crs)
stream_grd <- projectRaster(stream_grd, crs = crs)
valley_grd <- projectRaster(valley_grd, crs = crs)

#reproject vector datasets of interest
gage       <- st_transform(gage,       crs = st_crs(crs))
stream_shp <- st_transform(stream_shp, crs = st_crs(crs))
valley_shp <- st_transform(valley_shp, crs = st_crs(crs))

#plot for funzies 
mapview(r) + mapview(stream_shp) + mapview(valley_shp) + mapview(gage)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6.0 Define study reach -------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Identify stream reach ---------------------------------------------------------
#Write files to temp data
st_write(gage, paste0(temp_dir,'\\gage.shp'), append = F)
writeRaster(r, paste0(temp_dir, "\\r.tif"), overwrite = T)

#Create Flow accumulation raster
#Smooth DEM
wbt_fast_almost_gaussian_filter(
  input = "r.tif",
  output = "r_smooth.tif", 
  sigma = 1.8, 
  wd = temp_dir
)

#breach depressions
wbt_breach_depressions(
  dem =    "r_smooth.tif",
  output = "r_breached.tif",
  fill_pits = F,
  wd = temp_dir)

#Flow accumulation raster
wbt_d8_flow_accumulation(
  input = "r_breached.tif",
  output = "fac.tif",
  wd = temp_dir
)

#Snap gage to channel
wbt_snap_pour_points(
  pour_pts = "gage.shp", 
  flow_accum = "fac.tif", 
  output = "snap.shp",
  snap_dist = 100, 
  wd = temp_dir)

#pull into R env
gage_snap <- st_read(paste0(temp_dir,"//snap.shp"), crs = st_crs(crs))

#identify study reach
reach <- stream_shp[st_buffer(gage_snap, 5),]

mapview(reach) + mapview(gage_snap)

#valley reach ------------------------------------------------------------------
#Create XS at ends of study reach
xs_fun<-function(reach, width=1000){
  
  #Define points along reach
  reach_pnts <- st_cast(reach, "POINT")
  
  #Define reach end points
  reach_end_pnts <- bind_rows(reach_pnts[1,],reach_pnts[nrow(reach_pnts),])

  #Estimate planar slope along reach
  reach_coord<-st_coordinates(reach_pnts)
  reach_slope<-(reach_coord[1,"Y"]-reach_coord[nrow(reach_coord),"Y"])/(reach_coord[1,"X"]-reach_coord[nrow(reach_coord),"X"])
  
  #Estimate inverse slope
  xs_slope <- -1/reach_slope
  
  #Estimate endpoints of XS_1
  xs_coord_1 <- st_coordinates(reach_end_pnts[1,])
  xs_coord_1 <-rbind(
    xs_coord_1, 
    matrix(0, nrow=2, ncol=2)
  )
  xs_coord_1[2,"X"] <- xs_coord_1[1,"X"] + width/2*cos(atan(xs_slope))
  xs_coord_1[2,"Y"] <- xs_coord_1[1,"Y"] + width/2*sin(atan(xs_slope))
  xs_coord_1[3,"X"] <- xs_coord_1[1,"X"] - width/2*cos(atan(xs_slope))
  xs_coord_1[3,"Y"] <- xs_coord_1[1,"Y"] - width/2*sin(atan(xs_slope))
  xs_coord_1<-xs_coord_1[-1,]
  
  #Create XS
  xs_1<-xs_coord_1 %>%  
    as_tibble() %>% 
    st_as_sf(coords = c("X","Y")) %>% 
    st_coordinates() %>% 
    st_linestring() %>% 
    st_sfc(.) %>% 
    st_set_crs(st_crs(stream_shp)) %>% 
    st_as_sf() 
  
  #Estimate endpoints of XS_1
  xs_coord_2 <- st_coordinates(reach_end_pnts[2,])
  xs_coord_2 <-rbind(
    xs_coord_2, 
    matrix(0, nrow=2, ncol=2)
  )
  xs_coord_2[2,"X"] <- xs_coord_2[1,"X"] + width/2*cos(atan(xs_slope))
  xs_coord_2[2,"Y"] <- xs_coord_2[1,"Y"] + width/2*sin(atan(xs_slope))
  xs_coord_2[3,"X"] <- xs_coord_2[1,"X"] - width/2*cos(atan(xs_slope))
  xs_coord_2[3,"Y"] <- xs_coord_2[1,"Y"] - width/2*sin(atan(xs_slope))
  xs_coord_2<-xs_coord_2[-1,]
  
  #Create XS
  xs_2<-xs_coord_2 %>%  
    as_tibble() %>% 
    st_as_sf(coords = c("X","Y")) %>% 
    st_coordinates() %>% 
    st_linestring() %>% 
    st_sfc(.) %>% 
    st_set_crs(st_crs(stream_shp)) %>% 
    st_as_sf() 
  
  #Export XS Shape
  xs <- bind_rows(xs_1, xs_2)
}
xs <- xs_fun(reach)
xs <- st_combine(xs)
xs <- st_cast(xs, 'MULTILINESTRING')

#Split polygon by xs
valley_chopped_shp <- st_split(valley_shp, xs) %>%  st_collection_extract(c("POLYGON"))
valley_reach <- valley_chopped_shp[gage_snap,]

#Plot for funzies
mapview(valley_reach) + mapview(gage)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7.0 Extract metrics ----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Extract metrics





#Semi code ---------------------------------------------------------------------
#Snap gage to stream link
#Identify stream link and create XS at begining and end
#clip valley shape with XS
#Define Valley slope
#Define width -- perhaps average of ~100 XS or something similar
#Define valley length (area divided by width)
#Define levation change based on stream elevation at start and stop
#Define 
