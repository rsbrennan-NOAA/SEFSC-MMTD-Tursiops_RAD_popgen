#map with populations

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

pops <- read.table("analysis/population_assignments_summary.txt", header=T)

dat <- read.csv("Tursiops_RADseq_Metadata_new.csv")


out <- merge(dat, pops, by.x="Lab.ID", by.y="indiv")


#----------------------------------------------
# add pop colors

# coastal
#56B4E9
#004488

#Offshore:
#F0B800
#B65A00

#intermediate:
#1B9E77
#66A61E

p <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_point(data = out, 
             aes(x = Long, y = Lat, fill=fourpop, shape=fourpop),
             size = 2.5,
             alpha=1,
             color="black") +
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
  scale_shape_manual(values=c(21,22,23, 24)) +
  scale_fill_manual(values=c("#56B4E9","#004488","#66A61E","#F0B800")) +
  coord_sf(xlim = c(-100, -64), ylim = c(23, 42), expand = FALSE)+
  annotation_scale()

p

ggsave("figures/map_allpops.pdf", p, h=4, w=5)
ggsave("figures/map_allpops.png", p, h=4, w=5)


# gulf only
p <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_point(data = out, 
             aes(x = Long, y = Lat, fill=fourpop, shape=fourpop),
             size = 2,
             alpha=1,
             color="black") +
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
  scale_shape_manual(values=c(21,22,23, 24)) +
  coord_sf(xlim = c(-98, -79), ylim = c(23, 32), expand = FALSE)+
  annotation_scale()

p

ggsave("figures/map_gulf.pdf", p, h=4, w=5)
ggsave("figures/map_gulf.png", p, h=4, w=5)


# atlantic

p <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_point(data = out, 
             aes(x = Long, y = Lat, fill=fourpop, shape=fourpop),
             size = 2,
             alpha=1,
             color="black") +
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
  scale_shape_manual(values=c(21,22,23, 24)) +
  coord_sf(xlim = c(-83, -72), ylim = c(25, 40), expand = FALSE)+
  annotation_scale()

p

ggsave("figures/map_atlantic.pdf", p, h=5, w=5)
ggsave("figures/map_atlantic.png", p, h=5, w=5)



#--------------------------------------------------------------
# add the putative hybrids on to the map:

out_hyb <- out %>% filter(offshore_putative_hybrids == TRUE)

p <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_point(data = out, 
             aes(x = Long, y = Lat, fill=fourpop, shape=fourpop),
             size = 2,
             alpha=1,
             color="black") +
  geom_point(data = out_hyb, 
             aes(x = Long, y = Lat),
             fill="red", 
             shape=21,
             size = 4,
             alpha=1,
             color="black") +
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
  scale_shape_manual(values=c(21,22,23, 24)) +
  coord_sf(xlim = c(-100, -60), ylim = c(23, 47), expand = FALSE)+
  annotation_scale()

p

ggsave("figures/map_allpops_hybrids.pdf", p, h=4, w=5)
ggsave("figures/map_allpops_hybrids.png", p, h=4, w=5)




