#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Topographic controls on stream expansion and contraction
#Coder: Nate Jones (cnjones7@ua.edu)
#Date Created: 12/6/2023
#Purpose: Recreate Prancevich and Kirchner (https://doi.org/10.1029/2018GL081799)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.0 Setup workspace ----------------------------------------------------------
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
library(plotly)

#Create temp dir
temp_dir <- "C:\\WorkspaceR\\TheValleyRulesTheStream\\scratch\\"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.0 Gather data --------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Define watershed outlet
outlet<-tibble(
  lat = 33.7621789,   
  lon = -85.5955205) %>% 
  st_as_sf(
    coords = c("lon", "lat"), 
    crs = '+proj=longlat +datum=WGS84 +no_defs') %>% 
  st_transform(crs = 3160)

#Turn off 3D distance calcs
sf_use_s2(FALSE)

#Download DEM
dem <- get_elev_raster(st_buffer(outlet,1000), z=14)

#Export data to temp directory
writeRaster(dem,paste0(temp_dir, 'dem.tif'), overwrite=T)
st_write(outlet, paste0(temp_dir, "outlet.shp"), append=T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.0 Create fdr and fac ------------------------------------------------------
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
  out_type = "catchment area",
  wd = temp_dir
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4.0 Delineate Watershed-----------------------------------------------------------
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
mapview(fac)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5.0 Estimate flownetwork  ----------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Define threshold for channel head in m2
threshold <- 1.5 * 10000

#Extract streams based on threshold
wbt_extract_streams(
  flow_accum = "fac_shed.tif",
  output = "flow_net.tif",
  threshold = threshold, 
  wd = temp_dir)

#Create points along stream
streams <- raster(paste0(temp_dir,"flow_net.tif"))
streams <- rasterToPoints(streams) %>% 
  as_tibble() %>% 
  filter(flow_net == 1) %>% 
  st_as_sf(coords = c("x", "y"), crs = st_crs(fac))

#Plot for vizual inspection
mapview(streams)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6.0 Estimate scaling factors  ------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#6.1 Estimate alpha based on accumulation area threshold for channel heads -----
#Create function
fun <- function(n){
  #Create flow net
  wbt_extract_streams(
    flow_accum = "fac_shed.tif",
    output = "flow_net.tif",
    threshold = n, 
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
  stream_density <- stream_length/(sum(st_area(watershed)))
  
  #Create export 
  tibble(
    A_channel_head_m2 = n, 
    stream_length_m = stream_length %>% paste0() %>% as.numeric(),
    stream_density_m_m2 = stream_density %>% paste0() %>% as.numeric()) 
}

#Apply function
output <- lapply(
  X = seq(0.5*10000, 2.5*10000, 1000), 
  FUN = fun
) %>% bind_rows()

#export results for safe keeping
write.csv(output, "scratch//output_PaintRock.csv")
output<-read_csv("scratch//output_PaintRock.csv")

#ugly plot (for now)
alpha_plot <- output %>% 
  mutate(A_channel_head_ha = A_channel_head_m2/10000) %>% 
  # filter(A_channel_head_ha > 0.5) %>% 
  # filter(A_channel_head_ha < 1.5) %>% 
  ggplot(aes(x=A_channel_head_ha, y=stream_density_m_m2)) + 
  geom_point() + 
  theme_bw()+
  theme(
    axis.title = element_text(size = 14), 
    axis.text  = element_text(size = 10)
  ) + 
  xlab("Threshold Area [ha]") +
  ylab("Stream Density [m^2/m]") +
  scale_x_log10() + 
  scale_y_log10()

#plot for funzies
ggplotly(alpha_plot)

#Filter results to area of interest
output <- output %>% 
  mutate(A_channel_head_ha = A_channel_head_m2/10000) %>% 
  filter(A_channel_head_ha > 0.5) %>% 
  filter(A_channel_head_ha < 1.7)

alpha <- lm(log10(output$stream_density_m_m2)~log10(output$A_channel_head_m2))
summary(alpha)
alpha <- summary(alpha)$coefficients[2]*-1

#6.2 Estimate theta (scaling factor between slope and area) ------------------------
#Extract streams based on threshold
wbt_extract_streams(
  flow_accum = "fac_shed.tif",
  output = "flow_net.tif",
  threshold = threshold, 
  wd = temp_dir)

#Estimate curvature 
wbt_slope(
  dem = "dem_smooth.tif",
  output = "slope.tif", 
  wd = temp_dir
)
slope <- raster(paste0(temp_dir, "slope.tif"))

#Create points along stream
flow_grid <- raster(paste0(temp_dir,"flow_net.tif"))
flow_pnts <- rasterToPoints(flow_grid) %>% 
  as_tibble() %>% 
  filter(flow_net == 1) %>% 
  st_as_sf(coords = c("x", "y"), crs = st_crs(fac))

#Estimate fac and curvature for each point along stream
output <- flow_pnts %>% 
  mutate(
    contributing_area_m2 = raster::extract(fac, flow_pnts), 
    slope                = raster::extract(slope, flow_pnts)) 

#ugly plot (for now)
theta_plot <- output %>% 
  mutate(contributing_area_ha = contributing_area_m2/10000) %>% 
  ggplot(aes(x=contributing_area_ha, y=slope)) + 
  geom_point() + 
  theme_bw()+
  theme(
    axis.title = element_text(size = 14), 
    axis.text  = element_text(size = 10)
  ) + 
  xlab("Area [ha]") +
  ylab("Slope") +
  scale_x_log10() +
  scale_y_log10() 

ggplotly(theta_plot)

#esimate scaling coefficient
theta <- lm(log10(output$slope)~log10(output$contributing_area_m2))
summary(theta)
theta <- summary(theta)$coefficients[2]*-1
#6.3 Estimate delta (scaling factor between curviture and area) ---------------------
#Extract streams based on threshold
wbt_extract_streams(
  flow_accum = "fac_shed.tif",
  output = "flow_net.tif",
  threshold = threshold, 
  wd = temp_dir)

#Estimate curvature 
wbt_mean_curvature(
  dem = "dem_smooth.tif",
  output = "curv.tif", 
  wd = temp_dir
)
curv <- raster(paste0(temp_dir, "curv.tif"))
curv <- curv*-1

#Create points along stream
flow_grid <- raster(paste0(temp_dir,"flow_net.tif"))
flow_pnts <- rasterToPoints(flow_grid) %>% 
  as_tibble() %>% 
  filter(flow_net == 1) %>% 
  st_as_sf(coords = c("x", "y"), crs = st_crs(fac))

#Estimate fac and curvature for each point along stream
output <- flow_pnts %>% 
  mutate(
    contributing_area_m2 = raster::extract(fac, flow_pnts), 
    curvature            = raster::extract(curv, flow_pnts)) 

#ugly plot (for now)
output %>% 
  mutate(contributing_area_ha = contributing_area_m2/10000) %>% 
  mutate(curvature = curvature*1000) %>% 
  ggplot(aes(x=contributing_area_ha, y=curvature)) + 
  geom_point() + 
  theme_bw()+
  theme(
    axis.title = element_text(size = 14), 
    axis.text  = element_text(size = 10)
  ) + 
  xlab("Area [ha]") +
  ylab("Curviture [1/km]") +
  scale_x_log10() +
  scale_y_log10()

#esimate scaling coefficient
delta <- lm((output$curvature*1000)~log10(output$contributing_area_m2/(1000^2)))
summary(delta)
delta <- summary(delta)$coefficients[2]

#6.4 Estimate gamma based on Prancevich and Kirchner derivation -----------------
#Estimate beta
k1 <- 1.5
k2 <- -0.51
beta_predicted <- alpha/((k1*delta) + k2+ theta+1)
beta_predicted <- if_else(beta_predicted < 0.04, 0.05, beta_predicted)

#Estimate gamma
gamma <- 1 + theta + (alpha/beta_predicted)

#6.5 Export results ------------------------------------------------------------
tibble(
  site = "TAL",
  alpha,
  theta,
  delta,
  gamma, 
  beta_predicted
)

#site  alpha theta delta gamma beta_estimated
#TAL   0.418 0.405 -2.84  11.8         -0.124

