
### try bio-oracle data? its 0.05 degree. but not sure how near shore.
library(raster)
library(terra)
library(dplyr)

clim <- rast("analysis/environmental_data/bio-oracle/temp_mean_thetao_baseline_2000_2019_depthsurf_4e3e_1426_a71d_U1738005641788.nc")
#clim <- project(clim, "EPSG:4326")

plot(clim)

location <- read.csv("Tursiops_RADseq_Metadata_new.csv")
coords<-data.frame(lon=location$Long, lat=location$Lat)

library(maps)
head(coords)
plot(clim, xlim=c(-100, -63), ylim=c(20,45))
points(coords, pch=21, bg="grey65", col="black")
# all working as expected

#

#make df to store output:
dfout <- as.data.frame(matrix(ncol=5, nrow=nrow(coords)))
colnames(dfout) <- c("id", "lon", "lat", "sample_date", "annual_mean_temp")
dfout$id <- location$Lab.ID
dfout$lon <- location$Lon
dfout$lat <- location$Lat
dfout$sample_date <-location$Collection.Date.1

# pull out temperature
for(i in 1:nrow(dfout)){
  val<-terra::extract(x=clim, y=coords[i,])
  dfout$annual_mean_temp[i] <- val[1,2]
}


dfout$annual_mean_temp

dfout_oracle <- dfout

plot(clim, xlim=c(-100, -63), ylim=c(20,45))
points(coords, pch=21, bg="grey65", col="black")
points(coords[which(is.na(dfout_oracle$annual_mean_temp)),], pch=21, bg="red", col="black")

sum(is.na(dfout_oracle$annual_mean_temp))
# 64
# they all near shore I think. 
# zoom in
plot(clim, xlim=c(-79, -77), ylim=c(33.5,34.5))
points(coords, pch=21, bg="grey65", col="black")
points(coords[which(is.na(dfout_oracle$annual_mean_temp)),], pch=21, bg="grey65", col="orange")

# some do not fall in a grid. we can move them. Not ideal, but most are very close.
# this should work fine for temp, but maybe a problem for salinity? 
  # bc moving from river to coast... not sure we can do better

dfout_corrected <- dfout

#skip_index <- c(28) 

missing_index <- which(is.na(dfout$annual_mean_temp))
missingdat <- data.frame(lon = coords$lon[missing_index],
                         lat = coords$lat[missing_index])

clim_wgs84 <- clim
for(i in 1:length(missing_index)){
  cat("Starting index", i, "\n")
  # first, convert to SpatVector, then do the extraction
  sample_vect <- vect(missingdat[i,1:2], geom = c("lon", "lat"),
                      crs="WGS84")
  sample_ext <- distance(x=clim_wgs84, y=sample_vect)
  df_ext <- values(sample_ext)
  
  # returns distance in meters
  # the problem is that it identifies those even with NA. so get next closest
  minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[1],df_ext[,1])
  # get these coordinates:
  close_coords <- as.data.frame((crds(sample_ext)[minrast_index,]))
  close_coords_vect <- vect(t(close_coords), 
                            crs="WGS84") 
  # get the temp
  val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
  # some move to another empty cell. if this happens, go to 2nd match
  start_val <- 2
  while(is.na(val[1,2]) | val[1,2] > 50){
    cat("Depth still 0, starting start_val:",(start_val), "\n")
    minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[start_val],df_ext[,1])
    # get these coordinates:
    close_coords <- as.data.frame((crds(sample_ext)[minrast_index,]))
    close_coords_vect <- vect(t(close_coords), 
                              crs="WGS84") 
    val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
    start_val <- start_val + 1
  }
  
  dfout_corrected$annual_mean_temp[missing_index[i]] <- val[1,2]
  dfout_corrected$lon[missing_index[i]] <- close_coords[1,1]
  dfout_corrected$lat[missing_index[i]] <- close_coords[2,1]
  cat("done with index", i, "\n")
  cat("\n")
}


#check manually, because sometimes wacky things happen. 
# 38 especially. 40, 46
buffer <- 0.2
i=46
# get resolution
r <- res(clim_wgs84)

# round to nearest grid cell
round_to_grid <- function(x, res) {
  round(x/res) * res
}

plot(clim_wgs84, 
     xlim=c(round_to_grid(dfout$lon[missing_index[i]] - buffer*1.3, r[1]),
            round_to_grid(dfout$lon[missing_index[i]] + buffer*1.3, r[1])), 
     ylim=c(round_to_grid(dfout$lat[missing_index[i]] - buffer, r[2]),
            round_to_grid(dfout$lat[missing_index[i]] + buffer, r[2])))

e <- ext(clim_wgs84)
r <- res(clim_wgs84)

# Add grid lines
abline(v=seq(e[1], e[2], by=r[1]), col="black", lwd=0.5)
abline(h=seq(e[3], e[4], by=r[2]), col="black", lwd=0.5)

points(dfout[missing_index[i],c(2,3)], pch=21, lwd=2, col="red")
points(dfout_corrected[missing_index[i],c(2,3)], pch=21, lwd=2, col="orange")

dfout[missing_index[i],]
dfout_corrected[missing_index[i],]

# add corrected latlon to df
# but drop date
colnames(dfout_corrected) <- c("id", "lon_corrected", "lat_corrected","sample_date", "annual_mean_temp_corrected")
dfout_corrected <- dfout_corrected %>% select(-sample_date)


alldat <- merge(dfout, dfout_corrected, by="id")
nrow(alldat) == nrow(dfout)

write.csv(alldat, file="analysis/environmental_data/annual_mean_temp.csv", row.names = F, quote = F)



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# run other variables.

# bc I don't know what grids will have NA, need to again find closest values for the missing ones. 


library(terra)
library(dplyr)

#### write some functions to make this easier

# pull out intitial env variable
extract_env_variable <- function(raster_path, location_data, variable_name) {
  # Read raster data
  clim <- rast(raster_path)
  
  # make coordinates dataframe
  coords<-data.frame(lon=location_data$Long, lat=location_data$Lat)
  
  
  # Initialize output dataframe
  dfout <- data.frame(
    id = location_data$Lab.ID,
    lon = coords$lon,
    lat = coords$lat
  )
  dfout[[variable_name]] <- NA
  
  # Extract values
  for(i in 1:nrow(dfout)) {
    val <- terra::extract(x = clim, y = coords[i,])
    dfout[[variable_name]][i] <- val[1,2]
  }
  
  return(dfout)
}

# function to correct missing values
correct_missing_values <- function(raster_path, dfout, variable_name, extreme_val) {
  # read in env variables
  clim <- rast(raster_path)
  
  # get missing values
  missing_index <- which(is.na(dfout[[variable_name]]))
  
  # make new df for corrections
  missingdat <- data.frame(
    lon = dfout$lon[missing_index],
    lat = dfout$lat[missing_index]
  )
  
  dfout_corrected <- dfout
  df_startvals <- data.frame(index=missing_index, cells_searched = NA)
  # Process each missing value
  for(i in 1:length(missing_index)) {
    sample_vect <- vect(missingdat[i,1:2], geom = c("lon", "lat"), crs = "WGS84")
    sample_ext <- distance(x = clim, y = sample_vect)
    df_ext <- values(sample_ext)
    # returns distance in meters
    # the problem is that it identifies those even with NA. so get next closest
    minrast_index <- match(sort(df_ext[,1], decreasing = FALSE)[1], df_ext[,1])
    # get these coordinates:
    close_coords <- as.data.frame((crds(sample_ext)[minrast_index,]))
    close_coords_vect <- vect(t(close_coords), 
                              crs="WGS84") 
    # pull out the value
      val <- terra::extract(x = clim, y = close_coords_vect)
      # some move to another empty cell. if this happens, go to 2nd match
      start_val <- 2      
      while(is.na(val[1,2]) | val[1,2] > extreme_val){
        cat("Value still NA or Extreme: ",val[1,2], " Starting start_val:",(start_val), "\n")
        minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[start_val],df_ext[,1])
        # get these coordinates:
        close_coords <- as.data.frame((crds(sample_ext)[minrast_index,]))
        close_coords_vect <- vect(t(close_coords), 
                                  crs="WGS84") 
        val<-terra::extract(x=clim, y=close_coords_vect)
        start_val <- start_val + 1
      }
      # add to df
      dfout_corrected[[variable_name]][missing_index[i]] <- val[1,2]
      dfout_corrected$lon[missing_index[i]] <- close_coords[1,1]
      dfout_corrected$lat[missing_index[i]] <- close_coords[2,1]
      cat("done with index", i, "\n")
      cat("\n")
      df_startvals$cells_searched[i] <- start_val
    }
  
  
  # Rename columns with name + corrected
  colnames(dfout_corrected)[colnames(dfout_corrected) == "lon"] <- "lon_corrected"
  colnames(dfout_corrected)[colnames(dfout_corrected) == "lat"] <- "lat_corrected"
  colnames(dfout_corrected)[colnames(dfout_corrected) == variable_name] <- paste0(variable_name, "_corrected")
  
  return(dfout_corrected)
  return(df_startvals)
}



# location <- read.csv("Tursiops_RADseq_Metadata_new.csv")
#
#Initial extraction
temp_data <- extract_env_variable(
   raster_path = "analysis/environmental_data/bio-oracle/temp_mean_thetao_baseline_2000_2019_depthsurf_4e3e_1426_a71d_U1738005641788.nc",
   location_data = location,
   variable_name = "annual_mean_temp"
 )

# PLOT!!!!!
dfout <- temp_data
dfout$annual_mean_temp

plot(clim, xlim=c(-100, -63), ylim=c(20,45))
points(coords, pch=21, bg="grey65", col="black")
points(coords[which(is.na(dfout_oracle$annual_mean_temp)),], pch=21, bg="red", col="black")

sum(is.na(dfout_oracle$annual_mean_temp))
# 64
# they all near shore I think. 
# zoom in
plot(clim, xlim=c(-79, -77), ylim=c(33.5,34.5))
points(coords, pch=21, bg="grey65", col="black")
points(coords[which(is.na(dfout_oracle$annual_mean_temp)),], pch=21, bg="grey65", col="orange")

# some do not fall in a grid. we can move them. Not ideal, but most are very close.
# this should work fine for temp, but maybe a problem for salinity? 
# bc moving from river to coast... not sure we can do better


#skip_index <- c(28) 

# correct missing
temp_data_corrected <- correct_missing_values(
   raster_path = "analysis/environmental_data/bio-oracle/temp_mean_thetao_baseline_2000_2019_depthsurf_4e3e_1426_a71d_U1738005641788.nc",
   dfout = temp_data,
   variable_name = "annual_mean_temp",
   extreme_val = 50
   )



head(temp_data_corrected)
head(dfout_corrected)


#

























#----------------------------------------------------------------------------------



setwd("C:/Users/Reid.Brennan/Downloads")
dat <- terra::rast("oisst-avhrr-v02r01.20240901.nc")


dat <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")
as.date(dat$collection.dat)
years <- format(as.Date(dat$collection.dat, format="%m/%d/%Y"),"%Y")

uniq_years <- unique(years)




# read in data, get overlaps:

##### https://rfunctions.blogspot.com/2017/08/extracting-data-from-rasters-using.html
library(raster)
library(terra)

clim <- rast("analysis/environmental_variables/temperature_oisst/sst.week.mean.nc")
plot(clim)

time(clim)

location <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")

coords<-data.frame(lon=location$long, lat=location$lat)
# need to adjust these to match those used in the raster file
coords$lon <- 360 +(coords$lon)


library(maps)
head(coords)
plot(clim[[1]])
points(coords, pch=16)

# need to loop over the coords, get the dates, find the cooresponding temp

# adjust the date formats
location$date_correct <- as.Date(location$collection.date, "%m/%d/%Y")

# find closest week to each collection date
# if there's tie, should return the earlier week, which seems correct
date_index <- outer(location$date_correct, time(clim), `-`) |> abs() |> apply(1, which.min)

#make df to store output:
dfout <- as.data.frame(matrix(ncol=5, nrow=nrow(coords)))
colnames(dfout) <- c("id", "lon", "lat", "sample_date", "weekly_mean_temp")
dfout$id <- location$Sample
dfout$lon <- location$long
dfout$lat <- location$lat
dfout$sample_date <-location$date_correct
  
for(i in 1:length(date_index)){
  tmp_clim <- clim[[date_index[i]]]
  val<-terra::extract(x=tmp_clim, y=coords[i,])
  dfout$weekly_mean_temp[i] <- val[1,2]
}


# need to fix the above for those with NA
## first plot and zoom in the make sure there is no data for these grids:
missing_index <- which(is.na(dfout$weekly_mean_temp))
length(missing_index)
missingdat <- data.frame(lon = coords$lon[missing_index],
                         lat = coords$lat[missing_index])

climplt <- clim[[1]]
plot(climplt, xlim=c(277,278), ylim=c(26.5,27.5))
points(missingdat, pch=21, lwd=1)

# yes, definitely falling outside of a raster. 
# find the 2nd closest raster- they're very close, just outside.
missingdat$lon_new <- NA
missingdat$lat_new <- NA


for(i in 1:length(missing_index)){

  # first, convert to SpatVector, then do the extraction
  sample_vect <- vect(missingdat[i,1:2], geom = c("lon", "lat"),
                      crs="+proj=longlat +datum=WGS84")
  sample_ext <- distance(x=clim[[date_index[1]]], y=sample_vect)
  df_ext <- values(sample_ext)

  # returns distance in meters
  head(df_ext)
  # the problem is that it identifies those even with NA. so get next closest 
  minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[2],df_ext[,1])
  # get these coordinates:
  close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
  # take this index, and pull it from the original 
  # use the same closest date as above
  tmp_clim <- clim[[date_index[missing_index[i]]]]
  # use the new coords
  val<-terra::extract(x=clim_wgs84, y=close_coords,
                      method="bilinear")
  val<-terra::extract(x=clim_wgs84, y=close_coords)
  dfout$annual_mean_temp[missing_index[i]] <- val[1,2]
  
  # save the new lon lat
  missingdat$lon_new[i] <- close_coords$x
  missingdat$lat_new[i] <- close_coords$y
}

# plot and look at them, to make sure there are no mistakes
climplt <- clim[[1]]
plot(climplt, xlim=c(277,278), ylim=c(26.5,27.5))
points(missingdat, pch=19, lwd=1)
points(missingdat[,3:4], pch=19, lwd=1, col="red")

# all seems good
# write the tmp output
write.csv(dfout, file="analysis/environmental_variables/weekly_mean_temp.csv", 
              quote=F, row.names=F)




#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# weekly temperature anomaly
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# need to handle this differently, bc daily data. 
# find the date, select that day and 7 days prior. take mean anomaly

location <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")
# adjust the date formats
location$date_correct <- as.Date(location$collection.date, "%m/%d/%Y")
location$year <- format(as.Date(location$date_correct, format="%d/%m/%Y"),"%Y")


coords<-data.frame(lon=location$long, lat=location$lat)
# need to adjust these to match those used in the raster file
coords$lon <- 360 +(coords$lon)


dfout <- as.data.frame(matrix(ncol=7, nrow=nrow(coords)))
colnames(dfout) <- c("id", "lon", "lat", "sample_date", "anom_temp_mean", "anom_temp_median", "n_days")
dfout$id <- location$Sample
dfout$lon <- location$long
dfout$lat <- location$lat
dfout$sample_date <-location$date_correct

# note that 6 samples still fall just outside of the grid. move them to the same grid as above. 
#
for(i in 1:length(missing_index)){
  coords$lon[missing_index[i]] <-missingdat$lon_new[i]
  coords$lat[missing_index[i]] <- missingdat$lat_new[i]
}

# loop over each sample

for(i in 1:nrow(location)){
  # get the year index
  year_tmp <- location$year[i]
  # read in the data
  anom_nm <- paste0("analysis/environmental_variables/temperature_anom/sst.day.anom.", year_tmp, ".nc")
  anom_in <- rast(anom_nm)
  # get date
  date_tmp <- location$date_correct[i]
  date_index <- which(time(anom_in) == date_tmp)
  date_week <- c(date_index, date_index-seq(1, 6, 1))
  # get the anom vals for each day:
  val_tmp <- rep(NA, 7)  
  for(day in 1:length(date_week)){
    tmp_clim <- anom_in[[date_week[day]]]
    anom_tmp <- terra::extract(x=tmp_clim, y=coords[i,])
    val_tmp[day] <- anom_tmp[1,2]
  }
  # get median for the week and save to output
  dfout$anom_temp_mean[i] <- median(val_tmp, na.rm=T)
  dfout$anom_temp_median[i] <-mean(val_tmp, na.rm=T)
  dfout$n_days[i] <- sum(!is.na(val_tmp))
    
}

# check for missing data:
dfout[which(dfout$n_days<7),]

plot(dfout$anom_temp_mean, dfout$anom_temp_median)

# write output:
write.csv(dfout[,1:5], file="analysis/environmental_variables/weekly_anomaly_temp.csv", 
          quote=F, row.names=F)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# gather all data so far:
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
dat_depth <- read.csv("analysis/environmental_variables/depth_distance.csv")
dat_anom <- read.csv("analysis/environmental_variables/weekly_anomaly_temp.csv")
dat_temp <- read.csv("analysis/environmental_variables/weekly_mean_temp.csv")


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# long-term data, annual avg
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# ------------------------------
# temperature
# ------------------------------



# read in data, get overlaps:

##### https://rfunctions.blogspot.com/2017/08/extracting-data-from-rasters-using.html
library(raster)
library(terra)
library(maps)

clim <- rast("analysis/environmental_variables/annual_mean_sst_1995_2004/gom_95A4_t00_10.nc")
plot(clim)

location <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")

coords<-data.frame(lon=location$long, lat=location$lat)


plot(clim[[1]], ylim=c(20,32), xlim=c(-99, -79))
points(coords, pch=21, col="black", bg="red2", cex=1.3)


library(maps)
library(sf)

# Get states data
usa <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

# Convert raster to dataframe for ggplot
clim_df <- as.data.frame(clim, xy = TRUE)
names(clim_df)[3] <- "value"  # rename the value column

# Create the plot
p <- ggplot() +
  # Plot the raster data
  geom_tile(data = clim_df, aes(x = x, y = y, fill = value)) +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  # Add the points
  
  geom_point(data = coords, aes(x = lon, y = lat),
             fill="red2",
             shape= 21, color="black", size = 2.5,
             alpha=0.5) +
  # Set the coordinate limits
  coord_sf(xlim = c(-100, -78), ylim = c(23, 32), expand = FALSE)+
  # Customize the theme
  theme_bw() +
  theme(
    #panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    #axis.text = element_blank(),
    #axis.ticks = element_blank()
    legend.position = "top",
    legend.title=element_blank()) +
  # Customize the color scale for the raster
  scale_fill_viridis_c() +
  # Labels
  labs(x = "Longitude", y = "Latitude", fill = "SST") +
  theme(panel.grid = element_blank(),
        legend.position = "top") +
  annotation_scale()

ggsave("figures/map_annualSST_pops.pdf", p, h=4, w=5)

p <- ggplot() +
  # Plot the raster data
  geom_tile(data = clim_df, aes(x = x, y = y, fill = value)) +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  # Set the coordinate limits
  coord_sf(xlim = c(-100, -78), ylim = c(23, 32), expand = FALSE)+
  # Customize the theme
  theme_bw() +
  theme(
    #panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    #axis.text = element_blank(),
    #axis.ticks = element_blank()
    legend.position = "top",
    legend.title=element_blank()) +
  # Customize the color scale for the raster
  scale_fill_viridis_c() +
  # Labels
  labs(x = "Longitude", y = "Latitude", fill = "SST") +
  theme(panel.grid = element_blank(),
        legend.position = "top") +
  annotation_scale()

ggsave("figures/map_annualSST.pdf", p, h=4, w=5)
ggsave("figures/map_annualSST.png", p, h=4, w=5)


# need to loop over the coords, find the cooresponding temp

# adjust the date formats
location$date_correct <- as.Date(location$collection.date, "%m/%d/%Y")

#make df to store output:
dfout <- as.data.frame(matrix(ncol=5, nrow=nrow(coords)))
colnames(dfout) <- c("id", "lon", "lat", "sample_date", "annual_mean_temp")
dfout$id <- location$Sample
dfout$lon <- location$long
dfout$lat <- location$lat
dfout$sample_date <-location$date_correct

for(i in 1:nrow(dfout)){
  val<-terra::extract(x=clim[[1]], y=coords[i,])
  dfout$annual_mean_temp[i] <- val[1,2]
}

# need to fix the above for those with NA
## first plot and zoom in the make sure there is no data for these grids:
missing_index <- which(is.na(dfout$annual_mean_temp))
length(missing_index)
missingdat <- data.frame(lon = coords$lon[missing_index],
                         lat = coords$lat[missing_index])

climplt <- clim[[1]]
plot(climplt, xlim=c(-83,-82), ylim=c(26.5,27.5))
points(missingdat, pch=21, lwd=2, col="red")

# yes, definitely falling outside of a raster. 
# find the 2nd closest raster- they're very close, just outside.
missingdat$lon_new <- NA
missingdat$lat_new <- NA

for(i in 1:length(missing_index)){
  
  # first, convert to SpatVector, then do the extraction
  sample_vect <- vect(missingdat[i,1:2], geom = c("lon", "lat"),
                      crs="EPSG:4326")
  clim_wgs84 <- project(clim[[1]], "EPSG:4326")
  sample_ext <- distance(x=clim_wgs84, y=sample_vect)
  df_ext <- values(sample_ext)
  
  # returns distance in meters
  head(df_ext)
  # the problem is that it identifies those even with NA. so get next closest 
  minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[2],df_ext[,1])
  # get these coordinates:
  close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
  close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                            crs="EPSG:4326") 
  # take this index, and pull it from the original 
  # use the same closest date as above
  #tmp_clim <- clim[[1]][missing_index[i]]
  # use the new coords
  val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
  # some move to another empty cell. if this happens, go to 2nd match
  if(is.na(val[1,2]) == TRUE){
    minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[3],df_ext[,1])
    # get these coordinates:
    close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
    close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                              crs="EPSG:4326") 
    val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
    
  }
  dfout$annual_mean_temp[missing_index[i]] <- val[1,2]
  
  # save the new lon lat
  missingdat$lon_new[i] <- close_coords$x
  missingdat$lat_new[i] <- close_coords$y
}

# plot each, to make sure there are no mistakes
plot(climplt, xlim=c(-83,-82), ylim=c(26.5,27.5))
points(missingdat, pch=19, lwd=1)
points(missingdat[,3:4], pch=19, lwd=1, col="red")

# all seems good
# write the tmp output
write.csv(dfout, file="analysis/environmental_variables/annual_mean_temp.csv", 
          quote=F, row.names=F)






# ------------------------------
# salinity
# ------------------------------

# read in data, get overlaps:

##### https://rfunctions.blogspot.com/2017/08/extracting-data-from-rasters-using.html
library(raster)
library(terra)
library(maps)

clim <- rast("analysis/environmental_variables/annual_mean_salinity_1995_2004/gom_95A4_s00_10.nc")
plot(clim)

location <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")

coords<-data.frame(lon=location$long, lat=location$lat)

head(coords)
plot(clim[[1]])
points(coords, pch=16)

# need to loop over the coords, find the cooresponding temp

# adjust the date formats
location$date_correct <- as.Date(location$collection.date, "%m/%d/%Y")

#make df to store output:
dfout <- as.data.frame(matrix(ncol=5, nrow=nrow(coords)))
colnames(dfout) <- c("id", "lon", "lat", "sample_date", "annual_mean_temp")
dfout$id <- location$Sample
dfout$lon <- location$long
dfout$lat <- location$lat
dfout$sample_date <-location$date_correct

for(i in 1:nrow(dfout)){
  val<-terra::extract(x=clim[[1]], y=coords[i,])
  dfout$annual_mean_temp[i] <- val[1,2]
}

# need to fix the above for those with NA
## first plot and zoom in the make sure there is no data for these grids:
missing_index <- which(is.na(dfout$annual_mean_temp))
length(missing_index)
missingdat <- data.frame(lon = coords$lon[missing_index],
                         lat = coords$lat[missing_index])

climplt <- clim[[1]]
plot(climplt, xlim=c(-83,-82), ylim=c(26.5,27.5))
points(missingdat, pch=21, lwd=2, col="red")

# yes, definitely falling outside of a raster. 
# find the 2nd closest raster- they're very close, just outside.
missingdat$lon_new <- NA
missingdat$lat_new <- NA

for(i in 1:length(missing_index)){
  
  # first, convert to SpatVector, then do the extraction
  sample_vect <- vect(missingdat[i,1:2], geom = c("lon", "lat"),
                      crs="EPSG:4326")
  clim_wgs84 <- project(clim[[1]], "EPSG:4326")
  sample_ext <- distance(x=clim_wgs84, y=sample_vect)
  df_ext <- values(sample_ext)
  
  # returns distance in meters
  head(df_ext)
  # the problem is that it identifies those even with NA. so get next closest 
  minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[2],df_ext[,1])
  # get these coordinates:
  close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
  close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                            crs="EPSG:4326") 
  # take this index, and pull it from the original 
  # use the same closest date as above
  #tmp_clim <- clim[[1]][missing_index[i]]
  # use the new coords
  val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
  # some move to another empty cell. if this happens, go to 2nd match
  if(is.na(val[1,2]) == TRUE){
    minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[3],df_ext[,1])
    # get these coordinates:
    close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
    close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                              crs="EPSG:4326") 
    val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
    
  }
  dfout$annual_mean_temp[missing_index[i]] <- val[1,2]
  
  # save the new lon lat
  missingdat$lon_new[i] <- close_coords$x
  missingdat$lat_new[i] <- close_coords$y
}

# plot each, to make sure there are no mistakes
plot(climplt, xlim=c(-83,-82), ylim=c(26.5,27.5))
points(missingdat, pch=19, lwd=1)
points(missingdat[,3:4], pch=19, lwd=1, col="red")

# all seems good
# write the tmp output
write.csv(dfout, file="analysis/environmental_variables/annual_mean_salinity.csv", 
          quote=F, row.names=F)









# ------------------------------
# oxygen
# ------------------------------

# read in data, get overlaps:

clim <- rast("analysis/environmental_variables/annual_mean_oxygen/gom_all_o00_01.nc")
plot(clim)

location <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")

coords<-data.frame(lon=location$long, lat=location$lat)

head(coords)
plot(clim[[1]])
points(coords, pch=16)

# need to loop over the coords, find the cooresponding temp

# adjust the date formats
location$date_correct <- as.Date(location$collection.date, "%m/%d/%Y")

#make df to store output:
dfout <- as.data.frame(matrix(ncol=5, nrow=nrow(coords)))
colnames(dfout) <- c("id", "lon", "lat", "sample_date", "annual_mean_temp")
dfout$id <- location$Sample
dfout$lon <- location$long
dfout$lat <- location$lat
dfout$sample_date <-location$date_correct

for(i in 1:nrow(dfout)){
  val<-terra::extract(x=clim[[1]], y=coords[i,])
  dfout$annual_mean_temp[i] <- val[1,2]
}

# need to fix the above for those with NA
## first plot and zoom in the make sure there is no data for these grids:
missing_index <- which(is.na(dfout$annual_mean_temp))
length(missing_index)
missingdat <- data.frame(lon = coords$lon[missing_index],
                         lat = coords$lat[missing_index])

climplt <- clim[[1]]
plot(clim[[1]], xlim=c(-83,-82), ylim=c(26.5,27.5))
plot(clim[[1]])
points(missingdat, pch=21, lwd=2, col="red")

# yes, definitely falling outside of a raster. 
# find the 2nd closest raster- they're very close, just outside.
missingdat$lon_new <- NA
missingdat$lat_new <- NA

for(i in 1:length(missing_index)){
  
  # first, convert to SpatVector, then do the extraction
  sample_vect <- vect(missingdat[i,1:2], geom = c("lon", "lat"),
                      crs="EPSG:4326")
  clim_wgs84 <- project(clim[[1]], "EPSG:4326")
  sample_ext <- distance(x=clim_wgs84, y=sample_vect)
  df_ext <- values(sample_ext)
  
  # returns distance in meters
  head(df_ext)

  minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[2], df_ext[,1])
  # get these coordinates:
  close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
  close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                            crs="EPSG:4326") 
  val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
  
  # if na still, fix
  j <- 3  # Start with the second closest point
  # some move to another empty cell. if this happens, go to 2nd match
  while(is.na(val[1,2])) {
    minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[j], df_ext[,1])
    # get these coordinates:
    close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
    close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                              crs="EPSG:4326") 
    val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
    #  counter if we need to try next point
    if(is.na(val[1,2])) {
      j <- j + 1
    }
    
  }
  dfout$annual_mean_temp[missing_index[i]] <- val[1,2]
  
  # save the new lon lat
  missingdat$lon_new[i] <- close_coords$x
  missingdat$lat_new[i] <- close_coords$y
}

# plot each, to make sure there are no mistakes
plot(clim[[1]])
points(missingdat, pch=19, lwd=1)
points(missingdat[,3:4], pch=19, lwd=1, col="red")

# all seems good
# write the tmp output
write.csv(dfout, file="analysis/environmental_variables/annual_mean_oxygen.csv", 
          quote=F, row.names=F)









# ------------------------------
# annual_mean_nitrate 
# ------------------------------

# read in data, get overlaps:

clim <- rast("analysis/environmental_variables/annual_mean_nitrate/gom_all_n00_01.nc")
plot(clim)

location <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")

coords<-data.frame(lon=location$long, lat=location$lat)

head(coords)
plot(clim[[1]])
points(coords, pch=16)

# need to loop over the coords, find the cooresponding temp

# adjust the date formats
location$date_correct <- as.Date(location$collection.date, "%m/%d/%Y")

#make df to store output:
dfout <- as.data.frame(matrix(ncol=5, nrow=nrow(coords)))
colnames(dfout) <- c("id", "lon", "lat", "sample_date", "annual_mean_temp")
dfout$id <- location$Sample
dfout$lon <- location$long
dfout$lat <- location$lat
dfout$sample_date <-location$date_correct

for(i in 1:nrow(dfout)){
  val<-terra::extract(x=clim[[1]], y=coords[i,])
  dfout$annual_mean_temp[i] <- val[1,2]
}

# need to fix the above for those with NA
## first plot and zoom in the make sure there is no data for these grids:
missing_index <- which(is.na(dfout$annual_mean_temp))
length(missing_index)
missingdat <- data.frame(lon = coords$lon[missing_index],
                         lat = coords$lat[missing_index])

climplt <- clim[[1]]
plot(clim[[1]], xlim=c(-83,-82), ylim=c(26.5,27.5))
plot(clim[[1]])
points(missingdat, pch=21, lwd=2, col="red")

# yes, definitely falling outside of a raster. 
# find the 2nd closest raster- they're very close, just outside.
missingdat$lon_new <- NA
missingdat$lat_new <- NA

for(i in 1:length(missing_index)){
  
  # first, convert to SpatVector, then do the extraction
  sample_vect <- vect(missingdat[i,1:2], geom = c("lon", "lat"),
                      crs="EPSG:4326")
  clim_wgs84 <- project(clim[[1]], "EPSG:4326")
  sample_ext <- distance(x=clim_wgs84, y=sample_vect)
  df_ext <- values(sample_ext)
  
  # returns distance in meters
  head(df_ext)
  
  minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[2], df_ext[,1])
  # get these coordinates:
  close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
  close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                            crs="EPSG:4326") 
  val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
  
  # if na still, fix
  j <- 3  # Start with the second closest point
  # some move to another empty cell. if this happens, go to 2nd match
  while(is.na(val[1,2])) {
    minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[j], df_ext[,1])
    # get these coordinates:
    close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
    close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                              crs="EPSG:4326") 
    val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
    #  counter if we need to try next point
    if(is.na(val[1,2])) {
      j <- j + 1
    }
    
  }
  dfout$annual_mean_temp[missing_index[i]] <- val[1,2]
  
  # save the new lon lat
  missingdat$lon_new[i] <- close_coords$x
  missingdat$lat_new[i] <- close_coords$y
}

# plot each, to make sure there are no mistakes
plot(clim[[1]])
points(missingdat, pch=19, lwd=1)
points(missingdat[,3:4], pch=19, lwd=1, col="red")

# all seems good
# write the tmp output
write.csv(dfout, file="analysis/environmental_variables/annual_mean_nitrate.csv", 
          quote=F, row.names=F)





# ------------------------------
# annual_mean_phosphate 
# ------------------------------

# read in data, get overlaps:

clim <- rast("analysis/environmental_variables/annual_mean_phosphate/gom_all_p00_01.nc")
plot(clim)

location <- read.csv("data/GoMx_Tursiops_Snicro_FarFromShore-NoStranding.csv")

coords<-data.frame(lon=location$long, lat=location$lat)

head(coords)
plot(clim[[1]])
points(coords, pch=16)

# need to loop over the coords, find the cooresponding temp

# adjust the date formats
location$date_correct <- as.Date(location$collection.date, "%m/%d/%Y")

#make df to store output:
dfout <- as.data.frame(matrix(ncol=5, nrow=nrow(coords)))
colnames(dfout) <- c("id", "lon", "lat", "sample_date", "annual_mean_temp")
dfout$id <- location$Sample
dfout$lon <- location$long
dfout$lat <- location$lat
dfout$sample_date <-location$date_correct

for(i in 1:nrow(dfout)){
  val<-terra::extract(x=clim[[1]], y=coords[i,])
  dfout$annual_mean_temp[i] <- val[1,2]
}

# need to fix the above for those with NA
## first plot and zoom in the make sure there is no data for these grids:
missing_index <- which(is.na(dfout$annual_mean_temp))
length(missing_index)
missingdat <- data.frame(lon = coords$lon[missing_index],
                         lat = coords$lat[missing_index])

climplt <- clim[[1]]
plot(clim[[1]], xlim=c(-83,-82), ylim=c(26.5,27.5))
plot(clim[[1]])
points(missingdat, pch=21, lwd=2, col="red")

# yes, definitely falling outside of a raster. 
# find the 2nd closest raster- they're very close, just outside.
missingdat$lon_new <- NA
missingdat$lat_new <- NA

for(i in 1:length(missing_index)){
  
  # first, convert to SpatVector, then do the extraction
  sample_vect <- vect(missingdat[i,1:2], geom = c("lon", "lat"),
                      crs="EPSG:4326")
  clim_wgs84 <- project(clim[[1]], "EPSG:4326")
  sample_ext <- distance(x=clim_wgs84, y=sample_vect)
  df_ext <- values(sample_ext)
  
  # returns distance in meters
  head(df_ext)
  
  minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[2], df_ext[,1])
  # get these coordinates:
  close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
  close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                            crs="EPSG:4326") 
  val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
  
  # if na still, fix
  j <- 3  # Start with the second closest point
  # some move to another empty cell. if this happens, go to 2nd match
  while(is.na(val[1,2])) {
    minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[j], df_ext[,1])
    # get these coordinates:
    close_coords <- as.data.frame(t((crds(sample_ext)[minrast_index,])))
    close_coords_vect <- vect(close_coords, geom=c("x", "y"), 
                              crs="EPSG:4326") 
    val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
    #  counter if we need to try next point
    if(is.na(val[1,2])) {
      j <- j + 1
    }
    
  }
  dfout$annual_mean_temp[missing_index[i]] <- val[1,2]
  
  # save the new lon lat
  missingdat$lon_new[i] <- close_coords$x
  missingdat$lat_new[i] <- close_coords$y
}

# plot each, to make sure there are no mistakes
plot(clim[[1]])
points(missingdat, pch=19, lwd=1)
points(missingdat[,3:4], pch=19, lwd=1, col="red")

# all seems good
# write the tmp output
write.csv(dfout, file="analysis/environmental_variables/annual_mean_phosphate.csv", 
          quote=F, row.names=F)






