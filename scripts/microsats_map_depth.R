# microsat depths:


#library(ggplot2)
#library(sf)
#library(rnaturalearth)
#library(rnaturalearthdata)
#library(tidyverse)
#library(ggspatial)
library(marmap)
library(tidyverse)
#library(RColorBrewer)

micro_data <- read.csv("analysis/AllATL_metadata-FINAL_Ttru_only_All_samples_n4071.csv")
coords.gps = dplyr::select(micro_data, Long , Lat)

bathydata <- readGEBCO.bathy("analysis/gebco_2024_n53.8066_s20.4785_w-99.4219_e-44.4375.nc")

4Tt184=4Tt037=4Tt328
# Get depth of coordinates

depths = cbind(site = micro_data$Lab.ID,
               get.depth(bathydata, coords.gps, locator = FALSE))
hist(depths$depth, breaks=50)
hist(depths$depth, breaks=500, xlim=c(-100, (max(depths$depth)+10)))
max(depths$depth)

sum(depths$depth >= 0)
# 466

datmerge <- merge(depths, micro_data, by.x="site", by.y="Lab.ID")
head(datmerge)
dat_err <- datmerge[(depths$depth >= 0),]
dat_err <- dat_err[dat_err$Source != "stranding",]
nrow(dat_err)
#457

# figure out whats going on with these ones that are > 0
coords.gps2 = dplyr::select(dat_err, Long , Lat)

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
dfout$id <- micro_data$Lab.ID
dfout$lon <- coords.gps$Long
dfout$lat <- coords.gps$Lat

for(i in 1:nrow(coords.gps)){
  val<-terra::extract(x=clim, y=coords.gps[i,])
  dfout$depth_gebco[i] <- val[1,2]
}

missing_index <- which(dfout$depth_gebco >= 0)
length(missing_index)
# 37
missingdat <- data.frame(lon = coords.gps$Long[missing_index],
                         lat = coords.gps$Lat[missing_index])

buffer <- 0.005
i=1
plot(clim, xlim=c((missingdat$lon[i] - buffer),
                  (missingdat$lon[i] + buffer)), 
     ylim=c((missingdat$lat[i] - buffer),
            (missingdat$lat[i] + buffer)))
points(missingdat[i,], pch=21, lwd=2, col="red")

# these are all right by shore. just change to -1 for all. 
depths$depth[depths$depth >= 0] <- -1

nrow(depths)

hist(depths$depth, breaks=200)

colnames(depths) <- c("Lab.ID", "lon", "lat", "depth")

write.table(depths, file="analysis/microsat_depths.csv", row.names = F, quote = F,sep=",")

#--------------------------------------------------------
# plot, like fig 3 in ms:


dat1 <- read.csv("analysis/microsat_depths.csv")

# merge assignments from structure

all <- merge(dat1, micro_data, by="Lab.ID")

all$population <- NA
all$population[all$MicroStructure_Jul2022_n4068_ClumppK4.0.50.cutoff == "1"] <- "Offshore"
all$population[all$MicroStructure_Jul2022_n4068_ClumppK4.0.50.cutoff == "2"] <- "Intermediate"
all$population[all$MicroStructure_Jul2022_n4068_ClumppK4.0.50.cutoff == "3"] <- "Coastal_Gulf"
all$population[all$MicroStructure_Jul2022_n4068_ClumppK4.0.50.cutoff == "4"] <- "Coastal_Atlantic"
# assignments in col: MicroStructure_Jul2022_n4068_ClumppK4.0.50.cutoff
## 1 = offshore, 2 = intermediate, 3 = coastal gulf 4 = coastal atlantic
table(all$population)


df <- all %>%
  filter(depth != "stranding") %>%
  mutate(depth = as.numeric(depth),
         logdepth = log10(depth*-1)*-1) %>%
  mutate(logdepth = if_else(logdepth > 0, 0, logdepth))
df2 <- all %>%
  filter(depth != "stranding") %>%
  mutate(depth = if_else(depth > 0, 0, depth))

# need to split intermediat and offshore into atlantic and gulf:
# do this based on lat lon.

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(ggspatial)
library(marmap)
library(RColorBrewer)
library(ggpubr)

# Get the world map data

# Get the world map data
world <- ne_countries(scale = "medium", returnclass = "sf")


usa <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

df$sixpop <- df$population
df$sixpop[df$population == "Intermediate"] <- "Intermediate_Gulf"
df$sixpop[df$sixpop == "Intermediate_Gulf" & df$lon > -81.5] <- "Intermediate_Atlantic"
df$sixpop[df$sixpop == "Intermediate" & df$lon > -82 & df$lat < 26 ] <- "Intermediate_Atlantic"
df$sixpop[df$population == "Offshore"] <- "Offshore_Gulf"
df$sixpop[df$sixpop == "Offshore_Gulf" & df$lon > -81.5] <- "Offshore_Atlantic"

p1 <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_point(data = df, 
             aes(x = Long, y = Lat, fill=sixpop, shape=sixpop),
             size = 2.5,
             alpha=1,
             color="black") +
  coord_sf(xlim = c(-100, -64), ylim = c(23, 42), expand = FALSE)+
  scale_shape_manual(values=c(21,21,22,22,24,24))+
  annotation_scale()

p1



df2 %>%
  group_by(population) %>%
  summarise(
    n = n(),
    q0 = quantile(depth, 0, na.rm = TRUE),
    q01 = quantile(depth, 0.01, na.rm = TRUE),
    q05 = quantile(depth, 0.05, na.rm = TRUE),
    q50 = quantile(depth, 0.5, na.rm = TRUE),
    q95 = quantile(depth, 0.95, na.rm = TRUE),
    q99 = quantile(depth, 0.99, na.rm = TRUE),
    q100 = quantile(depth, 1, na.rm = TRUE),
    .groups = 'drop'
  )



summary_stats <- df %>%
  group_by(sixpop) %>%
  summarise(
    n = n(),
    mean_depth = mean(depth, na.rm = TRUE),
    min_depth = min(depth, na.rm = TRUE),
    #range_depth = max_depth - min_depth,
    #sd_depth = sd(depth, na.rm = TRUE),
    #se_depth = sd_depth / sqrt(n),
    quantile_01 = quantile(depth, 0.01, na.rm = TRUE),
    quantile_05 = quantile(depth, 0.05, na.rm = TRUE),
    quantile_50 = quantile(depth, 0.5, na.rm = TRUE),
    quantile_95 = quantile(depth, 0.95, na.rm = TRUE),
    quantile_99 = quantile(depth, 0.99, na.rm = TRUE),
    max_depth = max(depth, na.rm = TRUE),
    .groups = 'drop'
  ) 
  #dplyr::select(population, n, mean_depth,sd_depth, se_depth, min_depth, max_depth, range_depth, ci_95)
summary_stats

write.csv(summary_stats, file="depth_summary_micros.csv", row.names=F)

summary_stats$max_depth[2:6]

df%>% filter(sixpop == "Intermediate\nAtlantic") %>%
  filter(corrected_depth > -10)

y_breaks <- c(0, -log10(10), -log10(100), -log10(200), -log10(500), -log10(1000), -log10(2000), -log10(3000))
y_labels <- c("1", "10", "100", "200", "500", "1000", "2000", "3000")

df <- df %>%
  filter(!is.na(sixpop)) %>%
  mutate(sixpop = str_replace_all(sixpop, "_", "\n"))
  

fill_colors <- c(
  "Coastal\nAtlantic" = "#4782d4",
  "Coastal\nGulf" = "#e1526b",
  "Intermediate\nAtlantic" = "#B4ED50",
  "Intermediate\nGulf" = "#2E8B57", 
  "Offshore\nAtlantic" = "#FFDD33", 
  "Offshore\nGulf" = "#C49E45"
)


depth_plot <- ggplot(df, aes(x = sixpop, y = logdepth, fill=sixpop)) +
  geom_hline(yintercept = y_breaks, color = "grey80", linetype = "solid", linewidth = 0.3) +
  geom_violin() +
  #geom_swarm(overflow = "compress", color="black",size=0.5) +
  theme_classic(base_size = 14) +
  scale_y_continuous(
    breaks = y_breaks,
    labels = y_labels,
    limits = c(min(df$logdepth, na.rm = TRUE), max(df$logdepth, na.rm = TRUE))
  ) +
  labs(y = "Depth (m)",
       x = NULL) +
  scale_fill_manual(values=fill_colors) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
  )

depth_plot

ggsave(file="figures/depth_micros.png",depth_plot, h=5, w=7)
# look into the shallow ones. where are they? Highlight on the map.


df$sixpop <- df$population
df$sixpop[df$population == "Intermediate"] <- "Intermediate_Gulf"
df$sixpop[df$sixpop == "Intermediate_Gulf" & df$lon > -81.5] <- "Intermediate_Atlantic"
df$sixpop[df$sixpop == "Intermediate" & df$lon > -82 & df$lat < 26 ] <- "Intermediate_Atlantic"
df$sixpop[df$population == "Offshore"] <- "Offshore_Gulf"
df$sixpop[df$sixpop == "Offshore_Gulf" & df$lon > -81.5] <- "Offshore_Atlantic"


library(ggOceanMaps)
library(ggnewscale)

map_limits <- data.frame(lon = c(-100, -64), lat = c(23, 42))

breaks <- c(0, 200, 4000)

p1 <- basemap(data = map_limits, bathy.style = "rcb", grid.col = NA) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) colorRampPalette(c("#F7FBFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08306B"))(length(breaks)),
    limits = c(min(breaks), max(breaks)),
    breaks = breaks,
    show.limits = TRUE,
    guide = guide_coloursteps(title = "Depth")
  ) +
  new_scale_fill() +
  ggspatial::geom_spatial_point(data = df, 
                                aes(x = Long, y = Lat, fill = sixpop, shape = sixpop),
                                size = 2.5,
                                alpha = 1,
                                color = "black") +
  scale_fill_manual(values = fill_colors, name = "Population",
                    guide = guide_legend(
                      override.aes = list(fill_new = "white")
                    )) +
  scale_shape_manual(values = c(21,21,22,22,24,24), name = "Population") +
  guides(fill = guide_legend(override.aes = list(shape = c(21,21,22,22,24,24),
                                                 fill_new = "white"
                                                 ))) +
  annotation_scale() +
  theme(
    legend.key = element_blank()
  )

p1


# now zoom in:

map_limits <- data.frame(lon = c(-100, -78), lat = c(23, 32))

p2 <- basemap(data = map_limits, bathy.style = "rcb", grid.col = NA) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) colorRampPalette(c("#F7FBFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08306B"))(length(breaks)),
    limits = c(min(breaks), max(breaks)),
    breaks = breaks,
    show.limits = TRUE,
    guide = guide_coloursteps(title = "Depth")
  ) +
  new_scale_fill() +
  ggspatial::geom_spatial_point(data = df, 
                                aes(x = Long, y = Lat, fill = sixpop, shape = sixpop),
                                size = 2.5,
                                alpha = 1,
                                color = "black") +
  scale_fill_manual(values = fill_colors, name = "Population",
                    guide = guide_legend(
                      override.aes = list(fill_new = "white")
                    )) +
  scale_shape_manual(values = c(21,21,22,22,24,24), name = "Population") +
  guides(fill = guide_legend(override.aes = list(shape = c(21,21,22,22,24,24),
                                                 fill_new = "white"
  ))) +
  annotation_scale() +
  theme(
    legend.key = element_blank()
  )

p2




# atlantic

map_limits <- data.frame(lon = c(-82, -64), lat = c(23, 45))

p3 <- basemap(data = map_limits, bathy.style = "rcb", grid.col = NA) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) colorRampPalette(c("#F7FBFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08306B"))(length(breaks)),
    limits = c(min(breaks), max(breaks)),
    breaks = breaks,
    show.limits = TRUE,
    guide = guide_coloursteps(title = "Depth")
  ) +
  new_scale_fill() +
  ggspatial::geom_spatial_point(data = df, 
                                aes(x = Long, y = Lat, fill = sixpop, shape = sixpop),
                                size = 2.5,
                                alpha = 1,
                                color = "black") +
  scale_fill_manual(values = fill_colors, name = "Population",
                    guide = guide_legend(
                      override.aes = list(fill_new = "white")
                    )) +
  scale_shape_manual(values = c(21,21,22,22,24,24), name = "Population") +
  guides(fill = guide_legend(override.aes = list(shape = c(21,21,22,22,24,24),
                                                 fill_new = "white"
  ))) +
  annotation_scale() +
  theme(
    legend.key = element_blank()
  )

p3


# find the shallow intermediates. plot just those:


df_int <- df[df$sixpop == "Intermediate\nAtlantic" |
                   df$sixpop == "Intermediate\nGulf",]

nrow(df_int)
# 849

df_shallow <- df_int[df_int$depth > -10,]
nrow(df_shallow)
# 510
map_limits <- data.frame(lon = c(-100, -64), lat = c(23, 42))

p4 <- basemap(data = map_limits, bathy.style = "rcb", grid.col = NA) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) colorRampPalette(c("#F7FBFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08306B"))(length(breaks)),
    limits = c(min(breaks), max(breaks)),
    breaks = breaks,
    show.limits = TRUE,
    guide = guide_coloursteps(title = "Depth")
  ) +
  new_scale_fill() +
  ggspatial::geom_spatial_point(data = df_shallow, 
                                aes(x = Long, y = Lat, fill = sixpop, shape = sixpop),
                                size = 2.5,
                                alpha = 1,
                                color = "black") +
  scale_fill_manual(values = fill_colors, name = "Population",
                    guide = guide_legend(
                      override.aes = list(fill_new = "white")
                    )) +
  scale_shape_manual(values = c(22,22), name = "Population") +
  guides(fill = guide_legend(override.aes = list(shape = c(22,22),
                                                 fill_new = "white"
  ))) +
  annotation_scale() +
  theme(
    legend.key = element_blank()
  )

p4


# zoom in on S florida:

map_limits <- data.frame(lon = c(-83, -79), lat = c(24, 29))

p5 <- basemap(data = map_limits, bathy.style = "rcb", grid.col = NA) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) colorRampPalette(c("#F7FBFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08306B"))(length(breaks)),
    limits = c(min(breaks), max(breaks)),
    breaks = breaks,
    show.limits = TRUE,
    guide = guide_coloursteps(title = "Depth")
  ) +
  new_scale_fill() +
  ggspatial::geom_spatial_point(data = df_shallow, 
                                aes(x = Long, y = Lat, fill = sixpop, shape = sixpop),
                                size = 2.5,
                                alpha = 1,
                                color = "black") +
  scale_fill_manual(values = fill_colors, name = "Population",
                    guide = guide_legend(
                      override.aes = list(fill_new = "white")
                    )) +
  scale_shape_manual(values = c(22,22), name = "Population") +
  guides(fill = guide_legend(override.aes = list(shape = c(22,22),
                                                 fill_new = "white"
  ))) +
  annotation_scale() +
  theme(
    legend.key = element_blank()
  )

p5


map_limits <- data.frame(lon = c(-92, -82.5), lat = c(29, 31))

p6 <- basemap(data = map_limits, bathy.style = "rcb", grid.col = NA) +
  binned_scale(
    aesthetics = "fill",
    scale_name = "stepsn",
    palette = function(x) colorRampPalette(c("#F7FBFF", "#DEEBF7", "#9ECAE1", "#4292C6", "#08306B"))(length(breaks)),
    limits = c(min(breaks), max(breaks)),
    breaks = breaks,
    show.limits = TRUE,
    guide = guide_coloursteps(title = "Depth")
  ) +
  new_scale_fill() +
  ggspatial::geom_spatial_point(data = df_shallow, 
                                aes(x = Long, y = Lat, fill = sixpop, shape = sixpop),
                                size = 2.5,
                                alpha = 1,
                                color = "black") +
  scale_fill_manual(values = fill_colors, name = "Population",
                    guide = guide_legend(
                      override.aes = list(fill_new = "white")
                    )) +
  scale_shape_manual(values = c(22,22), name = "Population") +
  guides(fill = guide_legend(override.aes = list(shape = c(22,22),
                                                 fill_new = "white"
  ))) +
  annotation_scale(location="tl") +
  theme(
    legend.key = element_blank()
  )

p6


ggsave(file="figures/map_all_micros.png",p1, h=7, w=7)

bottomrow <- ggarrange(p2, p3, nrow = 1, ncol = 2, 
                         legend = "none")
ggsave(file="figures/map_subset-1_micros.png",p2, h=4, w=6)
ggsave(file="figures/map_subset-2_micros.png",p3, h=6, w=4)



left_column <- ggarrange(p4, p6, nrow = 2, ncol = 1, 
                         legend = "none")


pout_shallow <- ggarrange(left_column, p5, 
                          nrow = 1, ncol = 2, 
                          common.legend = TRUE)

ggsave(file="figures/map_shallow_micros.png",pout_shallow, h=6, w=8)
