#figure 1

#map with populations

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(ggspatial)
library(marmap)
#library(RColorBrewer)
library(ggpubr)


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

desired_order <- c("Coastal_Atlantic", "Offshore_Gulf", "Intermediate_Gulf", "Coastal_Gulf",
                   "Offshore_Atlantic", "Intermediate_Atlantic")

# Arrange data in this order
df_sorted <- out %>%
  arrange(factor(sixpop, levels = (desired_order)))


p1 <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_point(data = df_sorted, 
             aes(x = Long, y = Lat, fill=sixpop, shape=sixpop),
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
    legend.position = "left",
    legend.title=element_blank()) +
  xlab("Longitude")+
  ylab("Latitude") +
  #scale_shape_manual(values=c(21,22,23,24)) +
  #scale_fill_manual(values=c("#56B4E9","#004488","#66A61E","#F0B800")) +
  scale_shape_manual(values=c(21,21,22,22,24,24))+
  scale_fill_manual(values=c("#4782d4", "#e1526b","#B4ED50","#2E8B57", "#FFDD33", "#C49E45")) +
  coord_sf(xlim = c(-100, -64), ylim = c(23, 42), expand = FALSE)+
  annotation_scale()
  

p1

ggsave("figures/map_allpops.pdf", p1, h=4, w=5)
ggsave("figures/map_allpops.png", p1, h=4, w=5)



p1_facet <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_point(data = df_sorted, 
             aes(x = Long, y = Lat, fill=sixpop, shape=sixpop),
             size = 2.5,
             alpha=1,
             color="black") +
  facet_wrap(~sixpop, nrow=3) +
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
  #scale_shape_manual(values=c(21,22,23,24)) +
  #scale_fill_manual(values=c("#56B4E9","#004488","#66A61E","#F0B800")) +
  scale_shape_manual(values=c(21,21,22,22,24,24))+
  scale_fill_manual(values=c("#4782d4", "#e1526b","#B4ED50","#2E8B57", "#FFDD33", "#C49E45")) +
  coord_sf(xlim = c(-100, -64), ylim = c(23, 42), expand = FALSE)+
  annotation_scale()


p1_facet

ggsave("figures/manuscript/map_facet.pdf", p1_facet, h=9, w=7)
ggsave("figures/manuscript/map_facet.png", p1_facet, h=9, w=7)


# add the hybrids:

pops1 <- read.table("analysis/population_assignments_hybrids_summary.txt", header=T)
pops2 <- pops1[!is.na(pops1$newhybrids_category),]
out_hyb <- df_sorted[df_sorted$Lab.ID %in%pops2$indiv,]
nrow(out_hyb)

p1_facet <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_point(data = df_sorted, 
             aes(x = Long, y = Lat, fill=sixpop, shape=sixpop),
             size = 2.5,
             alpha=1,
             color="black") +
  geom_point(data = out_hyb, 
             aes(x = Long, y = Lat),
             color="purple", 
             fill=NA,
             shape=1,
             size = 3.5,
             stroke=1,
             alpha=1) +
  facet_wrap(~sixpop, nrow=3) +
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
  #scale_shape_manual(values=c(21,22,23,24)) +
  #scale_fill_manual(values=c("#56B4E9","#004488","#66A61E","#F0B800")) +
  scale_shape_manual(values=c(21,21,22,22,24,24))+
  scale_fill_manual(values=c("#4782d4", "#e1526b","#B4ED50","#2E8B57", "#FFDD33", "#C49E45")) +
  coord_sf(xlim = c(-100, -64), ylim = c(23, 42), expand = FALSE)+
  annotation_scale() 


p1_facet

ggsave("figures/manuscript/map_facet_hybrids.pdf", p1_facet, h=9, w=7)
ggsave("figures/manuscript/map_facet_hybrids.png", p1_facet, h=9, w=7)




#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#B: PCA
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#PCA

library(ggplot2)
library(dplyr)
library(SNPRelate)
library(SeqArray)
library(stringr)

pops <- read.table("analysis/population_assignments_summary.txt", header=T)

nrow(pops)
filename = "analysis/variants/filtered.final_ids_LDthin"
filename.gds = paste0(filename, ".gds")
filename.vcf.gz = paste0(filename, ".vcf.gz")
# Convert VCF to GDS
# run this the first time only
SeqArray::seqVCF2GDS(vcf.fn = filename.vcf.gz, out.fn = filename.gds, storage.option="ZIP_RA")

#gdsin = SeqArray::seqOpen(filename.gds)
#SeqArray::seqClose(filename.gds)

print(paste0("The number of SAMPLES in data: ", length(c(SeqArray::seqGetData(gdsin, "sample.id")))))
# [1] "The number of SAMPLES in data: 345"
print(paste0("The number of SNPs in data: ",  length(c(SeqArray::seqGetData(gdsin, "variant.id")))))
# [1] "The number of SNPs in data: 4496"
summary(m1 <- SeqArray::seqMissing(gdsin, per.variant=TRUE))
summary(m2 <- SeqArray::seqMissing(gdsin, per.variant=FALSE))

samples <- SeqArray::seqGetData(gdsin, "sample.id")

samples[!samples %in% pops$Lab.ID]

#pop.ids <- str_extract(samples, "^[^Tt]+")
#pop.ids <- sub("^(FB).*$", "\\1", pop.ids)

SeqArray::seqSummary(gdsin)

SeqArray::seqSetFilter(gdsin)
snp_in_all <- SeqArray::seqGetData(gdsin, "variant.id")
samp_in <- SeqArray::seqGetData(gdsin, "sample.id")

# remove sex chr
sexsnp <- snp_in_all[which(SeqArray::seqGetData(gdsin, "chromosome") == "NC_047055.1")]
length(sexsnp)

snp_in <- as.character(snp_in_all[which(SeqArray::seqGetData(gdsin, "chromosome") != "NC_047055.1")])
length(snp_in)


pca.out = SNPRelate::snpgdsPCA(autosome.only = F, gdsin, num.thread=2, 
                               remove.monosnp = F, maf = 0,
                               snp.id=snp_in,
                               missing.rate = NaN) # filtering for pruned SNPs


eig = pca.out$eigenval[!is.na(pca.out$eigenval)]
barplot(100*eig/sum(eig), main="PCA Eigenvalues")

# scale 
eig <- 100*eig/sum(eig)

dat <- as.data.frame(pca.out$eigenvect)

colnames(dat) <- c("PC1", "PC2", "PC3",colnames(dat)[4:ncol(dat)] )
dat$IDs <- pca.out$sample.id
#dat$Population <- pop.ids

# repeat with dapc populations
df <- merge(dat, pops, by.x="IDs", by.y="indiv")


# coastal
#56B4E9
#004488

#Offshore:
#F0B800
#B65A00

#intermediate:
#1B9E77
#66A61E


# with shapes
p2 <- ggplot(df, aes(PC1, PC2, fill=sixpop, shape=sixpop)) +
  geom_point(size=3, color="black") +
  xlab(paste0("PC1: ",round(eig[1], 0),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 0),"% variance")) +
  theme_classic() +
  #ggtitle('PCA: DAPC populations')+
  scale_shape_manual(values=c(21,21,22,22,24,24))+
  scale_fill_manual(values=c("#4782d4", "#e1526b","#B4ED50","#2E8B57", "#FFDD33", "#C49E45"))
  #guides(fill=guide_legend(override.aes=list(shape=21)))

p2

ggsave("figures/pca_fig1.pdf", p2, h=3.5, w=5.5)





#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# admixture plot

# admixture
library(tidyverse)
library(dplyr)
library(forcats)
library(patchwork)

pops <- read.table("analysis/population_assignments_summary.txt", header=T)

samplelist <- read_delim("analysis/variants/LDthin_noX_numCorrect.fam",
                         col_names = c("individual", "id2", "a", "b", "c", "d"),
                         delim=" ")

# read in all date, in a loop
## first create an empty dataframe
all_data <- tibble(individual=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

# then add all results to this
for (k in 2:4){
  data <- read_delim(paste0("analysis/pop_structure/LDthin_noX_numCorrect.",k,".Q"),
                     col_names = paste0("Q",seq(1:k)),
                     delim=" ")
  data$sample <- samplelist$individual
  data$k <- k
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k) -> data
  all_data <- rbind(all_data,data)
}

# order samples by population:

all_data <- merge(all_data, pops, by.x="sample", by.y="indiv")

#df$corrected_depth[df$Source == "stranding"] <- 0
#sampleorder <- df$Lab.ID[order(df$region, (as.numeric(df$corrected_depth)*-1))]
#df$corrected_depth[order((as.numeric(df$corrected_depth)*-1))]
#all_data$IDs <- as.factor(all_data$sample)
#all_data$IDs <- factor(all_data$IDs, levels=sampleorder)
#all_data$sample <- all_data$IDs

# within each population, order indivs by q val
new_dat<- all_data[all_data$k == 2 & all_data$Q == "Q1",]

new_dat2 <- new_dat %>%
  mutate(sample = fct_reorder2(sample, region, value)) %>%
  arrange(region, value)
head(new_dat2)
all_data$sample <- factor(all_data$sample, levels=c(new_dat2$sample))
all_data$k <- as.numeric(all_data$k)

all_data$Q[which(all_data$k == "2" & all_data$Q == "Q2")] <- "Coastal"
all_data$Q[which(all_data$k == "2" & all_data$Q == "Q1")] <- "Offshore"
all_data$Q[which(all_data$k == "3" & all_data$Q == "Q1")] <- "Offshore"
all_data$Q[which(all_data$k == "3" & all_data$Q == "Q2")] <- "Coastal_Gulf"
all_data$Q[which(all_data$k == "3" & all_data$Q == "Q3")] <- "Coastal"
all_data$Q[which(all_data$k == "4" & all_data$Q == "Q4")] <- "Coastal_Gulf"
all_data$Q[which(all_data$k == "4" & all_data$Q == "Q1")] <- "Intermediate"
all_data$Q[which(all_data$k == "4" & all_data$Q == "Q3")] <- "Offshore"
all_data$Q[which(all_data$k == "4" & all_data$Q == "Q2")] <- "Coastal"



p3 <-  
  all_data %>%
  filter(k < 8) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  geom_rug(aes(x=sample, y=value, color=sixpop),
           linewidth = 1,
           length= unit(0.06, "npc"),
           sides="b",
           outside = TRUE) +
  coord_cartesian(clip = "off")+
  xlab("") + ylab("Ancestry") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_blank(),     # Remove the facet label text
        strip.background = element_blank()) +
  scale_color_manual(values = c("black","grey55"),
                     name = "Population") +
  scale_color_manual(values=c("#4782d4", "#e1526b","#B4ED50","#2E8B57", "#FFDD33", "#C49E45")) +
  scale_fill_manual(values=c("#4782d4", "#e1526b","#61BA5C", "#E2BF3C"))+
  facet_wrap(~k,ncol=1)
p3


#combined_plot <- wrap_plots(p, p2, heights = c(0.15, 1), ncol=1)
#combined_plot

#ggsave("figures/Admixture_plot.pdf", combined_plot, width = 7, height = 8, units="in")
#ggsave("figures/Admixture_plot.png", combined_plot, width = 7, height = 8, units="in")

p_out1 <- ggarrange(p2, p3, ncol = 2, nrow = 1, labels=c("C", "D"), 
                    common.legend = T, legend = "none",
                    widths = c(0.5, 0.5))
# add in morphology:

merged_micro_morph_PCA_unique <- read.csv("analysis/micro_morph_PCA_plotting.csv")

none_data <- merged_micro_morph_PCA_unique[is.na(merged_micro_morph_PCA_unique$ClumppK4.0.50.cutoff), ]
other_data <- merged_micro_morph_PCA_unique[!(is.na(merged_micro_morph_PCA_unique$ClumppK4.0.50.cutoff)), ]
table(other_data$ClumppK4.0.50.cutoff)

fill_colors <- c(
  "Coastal\nAtlantic" = "#4782d4",
  "Coastal\nGulf" = "#e1526b",
  "Intermediate\nAtlantic" = "#B4ED50",
  "Intermediate\nGulf" = "#2E8B57", 
  "Offshore\nAtlantic" = "#FFDD33", 
  "Offshore\nGulf" = "#C49E45"
)
none_data$group <- "5"  
dat <- read.csv("metadata_FINAL.csv")


other_data$region <- "Atlantic"
other_data$region[other_data$Long < -81] <- "Gulf"
other_data$fourpop <- NA
other_data$fourpop[other_data$ClumppK4.0.50.cutoff =="1"] <- "Offshore"
other_data$fourpop[other_data$ClumppK4.0.50.cutoff =="2"] <- "Coastal_Gulf"
other_data$fourpop[other_data$ClumppK4.0.50.cutoff =="3"] <- "Coastal_Atlantic"
other_data$fourpop[other_data$ClumppK4.0.50.cutoff =="4"] <- "Intermediate"

other_data$sixpop <- paste(other_data$fourpop, other_data$region, sep= "_")

table(other_data$sixpop)

other_data$sixpop <- gsub("Coastal_Atlantic_Atlantic", "Coastal_Atlantic",other_data$sixpop)
other_data$sixpop <- gsub("Coastal_Gulf_Gulf", "Coastal_Gulf",other_data$sixpop)
other_data <- other_data[order(other_data$sixpop == "Coastal_Gulf"), ]
write.table(other_data, file="analysis/micro_morph_toShare.txt", quote = F, row.names = F, sep="\t")

other_data[other_data$sixpop == "Intermediate_Gulf",]

# flip the x axis so it corresponds to genetics pca

none_data$PC1 <- none_data$PC1*-1
other_data$PC1 <- other_data$PC1*-1

p_m_sixpop <- ggplot() +
  geom_point(data = none_data, aes(x=PC1, y=PC2, fill=group, shape=group), color="grey",size = 3) +
  geom_point(data = other_data, aes(x=PC1, y=PC2, fill=sixpop, shape=sixpop),color="black",size = 3) +
  #geom_point(data=pccoord_subset_CG, aes(x=PC1, y=PC2, shape=RAD_fourpop, fill=RAD_fourpop, size="RADseq"),
  #           color="black", shape=21, fill="#e1526b") +
  #geom_point(data=pccoord_subset_CA, aes(x=PC1, y=PC2, shape=RAD_fourpop, fill=RAD_fourpop, size="RADseq"),
  #           color="black", shape=21, fill="#4782d4") +
  scale_fill_manual(values = c( "Coastal_Atlantic" = "#4782d4", "Coastal_Gulf" = "#e1526b", 
                                "Intermediate_Atlantic" = "#B4ED50", "Intermediate_Gulf" = "#2E8B57", 
                                "Offshore_Atlantic" = "#FFDD33",
                                "5" = "grey"),
                    labels = c("Coastal_Atlantic" = "Coastal Atlantic", "Coastal_Gulf" = "Coastal Gulf", 
                               "Intermediate_Atlantic" = "Intermediate Atlantic", "Intermediate_Gulf" = "Intermediate Gulf", 
                               "Offshore_Atlantic" = "Offshore Atlantic",
                               "5" = "Not genotyped"),
                    breaks = c("Coastal_Atlantic", "Coastal_Gulf", "Intermediate_Atlantic", "Intermediate_Gulf",
                               "Offshore_Atlantic", "5"),
                    name = "Genotype group") +
  
  scale_shape_manual(values=c("Coastal_Atlantic"=21, "Coastal_Gulf"=21, 
                              "Intermediate_Atlantic" = 22, "Intermediate_Gulf" = 22, 
                              "Offshore_Atlantic" = 24,
                              "5"=16),
                     labels = c("Coastal_Atlantic" = "Coastal Atlantic", "Coastal_Gulf" = "Coastal Gulf", 
                                "Intermediate_Atlantic" = "Intermediate Atlantic", "Intermediate_Gulf" = "Intermediate Gulf", 
                                "Offshore_Atlantic" = "Offshore Atlantic",
                                "5" = "Not genotyped"),                    
                     breaks = c("Coastal_Atlantic", "Coastal_Gulf", "Intermediate_Atlantic", "Intermediate_Gulf",
                                "Offshore_Atlantic", "5"),
                     name = "Genotype group") +
  #scale_size_manual(values = c("microsatellite" = 3, "RADseq" = 6),
  #                  name = "Genotype method") +
  theme_classic() +
  guides(fill = guide_legend(override.aes = list(shape = c(21,21,22,22,24,16))),
         shape = guide_legend(override.aes = list(fill = c("#4782d4", "#e1526b","#B4ED50", "#2E8B57", "#FFDD33", "grey"),
                                                  size=4)),
         size = guide_legend(override.aes = list(fill = "black", color = "black", shape = 21)))+
  scale_y_continuous(labels = function(x) sprintf("%.2f", x))

p_m_sixpop


# now merge
p_no_legend <- p_m_sixpop + theme(legend.position = "none")

p1_no_legend <- p1 + theme(legend.position = "none")

#p_out1 <- ggarrange(p2, p3, ncol = 2, nrow = 1, labels=c("C", "D"), 
#                    common.legend = T, legend = "none",
#                    widths = c(0.5, 0.5))

p_out1 <- ggarrange(p3,p_no_legend, ncol = 2, nrow = 1, 
                    labels=c("C", "D"), 
                    legend = "none",
                    heights = c(0.5, 0.5),
                    widths = c(0.6, 0.4))

p_combine1 <- ggarrange(p1_no_legend, p2,
                     labels=c("A", "B"),
                     legend = "none",
                     heights = c(0.5, 0.5),
                     widths = c(0.6, 0.4))

p_combine <- ggarrange(p_combine1,p_out1, nrow=2, labels=c("", "", ""), common.legend=F,
                   heights = c(0.5, 0.5))
p_combine

ggsave(filename="figures/Fig_01-v4.pdf", p_combine, h=135, w=170, units= "mm")
ggsave(filename="figures/Fig_01-v4.png", p_combine, h=5, w=7)


