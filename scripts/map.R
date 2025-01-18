#map

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(ggspatial)
library(marmap)
library(RColorBrewer)



# Get the world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# read in locations
dat <- read.csv("Tursiops_RADseq_Metadata_new.csv")

head(dat)

usa <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))


# 
p <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_point(data = dat, aes(x = Long, y = Lat),
             fill="sienna4",
             shape= 21, color="black", size = 2.5,
             alpha=0.8) +
  coord_sf() +
  theme_bw() +
  theme(
    #panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    #axis.text = element_blank(),
    #axis.ticks = element_blank()
    legend.position = "top",
    legend.title=element_blank()) +
  xlab("Longitude")+
  ylab("Latitude") +
  coord_sf(xlim = c(-100, -60), ylim = c(23, 47), expand = FALSE)+
  annotation_scale()
  
p

ggsave("figures/map.pdf", p, h=4, w=5)
ggsave("figures/map.png", p, h=4, w=5)


#----------------------------------------------
# add library colors

dat$RADSeq.Dataset.2 <- as.factor(dat$RADSeq.Dataset)
dat$RADSeq.Dataset.2 <- factor(dat$RADSeq.Dataset.2, levels=c("2019", "2018", "2018 & 2019"))
levels(dat$RADSeq.Dataset.2)

p <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_point(data = subset(dat, RADSeq.Dataset.2 == "2019"), 
             aes(x = Long, y = Lat, fill=RADSeq.Dataset.2, shape=RADSeq.Dataset.2),
             size = 2,
             shape=21,
             alpha=1) +
  # Then plot 2018 points on top
  geom_point(data = subset(dat, RADSeq.Dataset.2 == "2018"), 
             aes(x = Long, y = Lat, fill=RADSeq.Dataset.2, shape=RADSeq.Dataset.2),
             size = 2,
             shape=22,
             alpha=1) +
  # Finally plot 2018 & 2019 points if needed
  geom_point(data = subset(dat, RADSeq.Dataset.2 == "2018 & 2019"), 
             aes(x = Long, y = Lat, fill=RADSeq.Dataset.2, shape=RADSeq.Dataset.2),
             size = 2,
             shape=23,
             alpha=1) +
  coord_sf() +
  theme_bw() +
  theme(
    #panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    #axis.text = element_blank(),
    #axis.ticks = element_blank()
    legend.position = "top",
    legend.title=element_blank()) +
  xlab("Longitude")+
  ylab("Latitude") +
  scale_shape_manual(values=c(21,22,23)) +
  coord_sf(xlim = c(-100, -60), ylim = c(23, 47), expand = FALSE)+
  annotation_scale()

p

ggsave("figures/map_library.pdf", p, h=4, w=5)
ggsave("figures/map_library.png", p, h=4, w=5)


#--------------#
#
# Calculate least-cost distances to shore
#
#--------------#

# from https://github.com/Tom-Jenkins/seascape_rda_tutorial/blob/master/2.Prepare_spatial_data/2.prepare_spatial_data.R

# Import coordinates of sites

coords.gps = dplyr::select(dat, Long , Lat)

str(coords.gps)

# Get bathymetry data from NOAA using marmap 
bathydata = getNOAA.bathy(lon1 = min(coords.gps$Long) -1,
                          lon2 = max(coords.gps$Long)+1,
                          lat1 = min(coords.gps$Lat)-1,
                          lat2 = max(coords.gps$Lat) +1,
                          resolution = 1)

# try GEBCO, higher resolution

bathydata <- readGEBCO.bathy("analysis/gebco_2024_n53.8066_s20.4785_w-99.4219_e-44.4375.nc")


# Get depth of coordinates

depths = cbind(site = dat$Lab.ID,
               get.depth(bathydata, coords.gps, locator = FALSE))
hist(depths$depth, breaks=50)
hist(depths$depth, breaks=300, xlim=c(-100, (max(depths$depth)+10)))
max(depths$depth)

sum(depths$depth >= 0)
# 25

datmerge <- merge(depths, dat, by.x="site", by.y="Lab.ID")

dat_err <- datmerge[(depths$depth > 0),]

# figure out whats going on with these ones that are > 0
coords.gps2 = dplyr::select(dat_err, Long , Lat)


# map and check values:
#plot(bathydata, image = TRUE, deep = TRUE, shallow = TRUE)

#plot.bathy(bathydata, image= TRUE, land = TRUE, n = 0,
#           bpal = list(c(0, max(depths$depth), "grey"),
#                       c(min(depths$depth), 0, "royalblue")))
#points(coords.gps2$Long, coords.gps2$Lat, pch = 21, 
#       bg ="orange", col = "black", cex = 2)



#bathy_df <- fortify.bathy(bathydata)

#ggplot(bathy_df, aes(x=x, y=y)) + coord_quickmap() +
#  geom_raster(aes(fill=z), data=bathy_df[bathy_df$z <= 0,]) +
#  scale_x_continuous(expand=c(0,0)) +
#  scale_y_continuous(expand=c(0,0))

coords.gps2


library(raster)
library(terra)
clim <- rast("analysis/gebco_2024_n53.8066_s20.4785_w-99.4219_e-44.4375.nc")
plot(clim)

points(coords.gps, pch=16)
points(coords.gps2, pch=16, col="orange")

coords.gps

#make df to store output:
dfout <- as.data.frame(matrix(ncol=4, nrow=nrow(coords.gps)))
colnames(dfout) <- c("id", "lon", "lat", "depth_gebco")
dfout$id <- dat$Lab.ID
dfout$lon <- coords.gps$Long
dfout$lat <- coords.gps$Lat

for(i in 1:nrow(coords.gps)){
  val<-terra::extract(x=clim, y=coords.gps[i,])
  dfout$depth_gebco[i] <- val[1,2]
}


head(dfout)

hist(dfout$depth_gebco)

sum(dfout$depth_gebco >= 0)
#33

plot(dfout$depth_gebco, depths$depth)
# they agree

missing_index <- which(dfout$depth_gebco >= 0)
length(missing_index)
missingdat <- data.frame(lon = coords.gps$Long[missing_index],
                         lat = coords.gps$Lat[missing_index])

buffer <- 0.02
i=1
plot(clim, xlim=c((missingdat$lon[i] - buffer),
                  (missingdat$lon[i] + buffer)), 
            ylim=c((missingdat$lat[i] - buffer),
                  (missingdat$lat[i] + buffer)))
points(missingdat[i,], pch=21, lwd=2, col="red")


# they're all mistakes. right next to shore for some, so in wrong grid. 
# or others are in marshes/estuaries where the smaller channels aren't in the data 
clim_wgs84 <- project(clim, "EPSG:4326")

dfout_corrected <- dfout

skip_index <- c(28) 

for(i in 1:length(missing_index)){
  if (i %in% skip_index) {
    cat("Skipping index and doing manual", i, "\n")
    start_val <- 130 # for index 28, closest is 130. this is super slow, so I do it manually here
    sample_vect <- vect(missingdat[i,1:2], geom = c("lon", "lat"),
                        crs="EPSG:4326")
    sample_ext <- distance(x=clim_wgs84, y=sample_vect)
    df_ext <- values(sample_ext)
    minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[start_val],df_ext[,1])
    # get these coordinates:
    close_coords <- as.data.frame((crds(sample_ext)[minrast_index,]))
    close_coords_vect <- vect(t(close_coords), 
                              crs="EPSG:4326") 
    val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
    
    
    dfout_corrected$depth_gebco[missing_index[i]] <- val[1,2]
    dfout_corrected$lon[missing_index[i]] <- close_coords[1,1]
    dfout_corrected$lat[missing_index[i]] <- close_coords[2,1]
    next
  }
  
  cat("Starting index", i, "\n")
  # first, convert to SpatVector, then do the extraction
  sample_vect <- vect(missingdat[i,1:2], geom = c("lon", "lat"),
                      crs="EPSG:4326")
  sample_ext <- distance(x=clim_wgs84, y=sample_vect)
  df_ext <- values(sample_ext)
  
  # returns distance in meters
  # the problem is that it identifies those even with NA. so get next closest
  minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[2],df_ext[,1])
  # get these coordinates:
  close_coords <- as.data.frame((crds(sample_ext)[minrast_index,]))
  close_coords_vect <- vect(t(close_coords), 
                            crs="EPSG:4326") 
  # get the depth
  val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
  # some move to another empty cell. if this happens, go to 2nd match
  start_val <- 3
  while(val[1,2] >= 0){
    cat("Depth still 0, starting start_val:",(start_val), "\n")
    minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[start_val],df_ext[,1])
    # get these coordinates:
    close_coords <- as.data.frame((crds(sample_ext)[minrast_index,]))
    close_coords_vect <- vect(t(close_coords), 
                              crs="EPSG:4326") 
    val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
    start_val <- start_val + 1
  }
  
  dfout_corrected$depth_gebco[missing_index[i]] <- val[1,2]
  dfout_corrected$lon[missing_index[i]] <- close_coords[1,1]
  dfout_corrected$lat[missing_index[i]] <- close_coords[2,1]
  cat("done with index", i, "\n")
  cat("\n")
}

# for index 28, closest is 130. this is super slow, so I do it manually here


# 33 needs to move south to stay in the same river:
i = 33
sample_vect <- vect(data.frame(lon=c(-79.964),lat=c(32.787 )), geom = c("lon", "lat"),
                    crs="EPSG:4326")

sample_ext <- distance(x=clim_wgs84, y=sample_vect)
df_ext <- values(sample_ext)
minrast_index <- match(sort(df_ext[,1], decreasing=FALSE)[1],df_ext[,1])
# get these coordinates:
close_coords <- as.data.frame((crds(sample_ext)[minrast_index,]))
close_coords_vect <- vect(t(close_coords), 
                          crs="EPSG:4326") 
val<-terra::extract(x=clim_wgs84, y=close_coords_vect)
dfout_corrected$depth_gebco[missing_index[i]] <- val[1,2]
dfout_corrected$lon[missing_index[i]] <- close_coords[1,1]
dfout_corrected$lat[missing_index[i]] <- close_coords[2,1]

points(dfout_corrected[missing_index[i],c(2,3)], pch=21, lwd=2, col="purple")



#check all manually, because sometimes wacky things happen. 

buffer <- 0.05
i=33
plot(clim, xlim=c((dfout$lon[missing_index[i]] - buffer),
                  (dfout$lon[missing_index[i]] + buffer)), 
     ylim=c((dfout$lat[missing_index[i]] - buffer),
            (dfout$lat[missing_index[i]] + buffer)))
points(dfout[missing_index[i],c(2,3)], pch=21, lwd=2, col="red")
points(dfout_corrected[missing_index[i],c(2,3)], pch=21, lwd=2, col="orange")
dfout[missing_index[i],]
dfout_corrected[missing_index[i],]

#
#
colnames(dfout_corrected) <- c("id", "lon_corrected", "lat_corrected", "corrected_depth")
alldat <- merge(dfout, dfout_corrected, by="id")

sum(alldat$depth_gebco == alldat$corrected_depth)
nrow(alldat)

write.csv(alldat, file="depths.csv", row.names = F)
















# one is > 0, so move it slightly. This is a bug bc of the grid
# this also causes problems with distance calcs
coords.gps$long[depths$depth>0] <- -82.58 # from -82.56

depths = cbind(site = dat$Sample,
               get.depth(bathydata, coords.gps, locator = FALSE))
depths



# depths in meters

# Plot bathymetry data and coordinates
#coords <- coords %>%
#  mutate(colors = case_when(
#    Pop.Structure.Location == "Atlantic" ~ "#eac435",
#    Pop.Structure.Location == "Dry Tortuga" ~ "#557fc3",
#    Pop.Structure.Location == "NGOMex" ~ "#03cea4",
#    Pop.Structure.Location == "WGOMex" ~ "#fb4d3d",
#    TRUE ~ "#999999"  # Default color if none of the above match
#  ))

# Create coordinate grids from dimnames
x <- as.numeric(dimnames(bathydata)[[1]])
y <- as.numeric(dimnames(bathydata)[[2]])
coords <- expand.grid(x = x, y = y)

# Convert matrix to vector and combine with coordinates
bathy_df <- data.frame(
  coords,
  depth = as.vector(bathydata)
)





depthp <- ggplot() +
  # Bathymetry contours
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_contour(data = bathy_df, 
               aes(x = x, y = y, z = depth),
               breaks = seq(-500, 0, by = 20),
               color = "grey70") +
  #geom_point(data = dat, aes(x = long, y = lat, fill=Pop),
  #           shape= 21, color="black", size = 2.5,
  #           alpha=0.5) +
  theme_bw() +
  theme(
    #panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    #axis.text = element_blank(),
    #axis.ticks = element_blank()
    legend.position = "top",
    legend.title=element_blank()) +
  xlab("Longitude") + 
  ylab("Latitude") +   
  coord_sf(xlim = c(-100, -78), ylim = c(23, 32), expand = FALSE)+
  annotation_scale() +
  scale_fill_manual(values = Colorsdf$color,guide = guide_legend(override.aes = list(alpha = 1, size = 2.5))) 

ggsave("figures/map_depths.pdf", depthp, h=4, w=5)
ggsave("figures/map_depths.png",depthp, h=4, w=5)



# calculate the distance to nearest shoreline:

d <- dist2isobath(bathydata, x=dat$long, y=dat$lat, isobath = 0)
# isobath=0 means find coastline.
# output distance is in meters


# Visualize the great circle distances
blues <- c("lightsteelblue4","lightsteelblue3","lightsteelblue2","lightsteelblue1")
plot(bathydata, image=TRUE, lwd=0.1, land=TRUE, bpal = list(c(0,max(bathydata),"grey"), c(min(bathydata),0,blues)))
plot(bathydata, deep=-200, shallow=-200, step=0, lwd=0.6, add=TRUE)
points(dat$long,dat$lat, pch=21, col="orange4", bg="orange2", cex=.8)
linesGC(d[2:3],d[4:5])

head(d)

# make df of depth and distance to shore and write

dfout <- data.frame(
                    id = dat$Sample,
                    depth = depths$depth,
                    distance_to_shore = d$distance
)

write.csv(dfout, file="analysis/environmental_variables/depth_distance.csv", 
          row.names=F, quote=F)



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# get pairwise geographic distances between indivs
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# min depth 1 stops them from crossing land
# these are slow to run, especially lc.dist
trans1 <- trans.mat(bathydata, min.depth = -1, max.depth = NULL)
lc_paths <- lc.dist(trans1, coords.gps, res = "path")
# Compute least-cost distances (km) matrix
lc_dist <- lc.dist(trans1, coords.gps, res = "dist")
save(lc_paths, file = "least_cost_paths.RData")
save(lc_dist, file = "least_cost_paths_dist.RData")

load("least_cost_paths.RData")
load("least_cost_paths_dist.RData")

plot.bathy(bathydata, image= TRUE, land = TRUE, n = 0,
           bpal = list(c(0, max(bathydata), "grey"),
                       c(min(bathydata), 0, "royalblue")))
lapply(lc_paths, lines, col = "orange", lwd = 2, lty = 1)

# Compute least-cost distances (km) matrix

# Convert to matrix, rename columns and rows, and export as csv file
lc_mat = as.matrix(lc_dist)
colnames(lc_mat) = as.vector(location$Sample)
rownames(lc_mat) = as.vector(location$Sample)
lc_mat

# calculate geographic straight line distances directly. to resolve the close together samples that are called 0.

library(geodist)
distmat <- geodist (coords, measure = "geodesic")
distmat_km <- distmat/1000
                
lc_mat[lc_mat == 0] <- distmat_km[lc_mat == 0]

lc_mat[1:10,1:10]

str(lc_mat)                                        

write.csv(lc_mat, file="analysis/environmental_variables/lc_distances_km.csv")






