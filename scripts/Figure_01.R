#figure 1


#A

#map with populations

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

p1 <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_point(data = out, 
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
  scale_fill_manual(values=c("#A6DDF0", "#276FBF","#B4ED50","#2E8B57", "#FFDD33", "#C49E45")) +
  coord_sf(xlim = c(-100, -64), ylim = c(23, 42), expand = FALSE)+
  annotation_scale()
  

p1

ggsave("figures/map_allpops.pdf", p1, h=4, w=5)
ggsave("figures/map_allpops.png", p1, h=4, w=5)




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
SeqArray::seqVCF2GDS(vcf.fn = filename.vcf.gz, out.fn = filename.gds, storage.option="ZIP_RA")

gdsin = SeqArray::seqOpen(filename.gds)
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

snp_in <- snp_in_all[which(SeqArray::seqGetData(gdsin, "chromosome") != "NC_047055.1")]
length(snp_in)


pca.out = SNPRelate::snpgdsPCA(autosome.only = F, gdsin, num.thread=2, 
                               remove.monosnp = T, maf = 0.05,
                               snp.id=snp_in) # filtering for pruned SNPs


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
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_classic() +
  #ggtitle('PCA: DAPC populations')+
  scale_shape_manual(values=c(21,21,22,22,24,24))+
  scale_fill_manual(values=c("#A6DDF0", "#276FBF","#B4ED50","#2E8B57", "#FFDD33", "#C49E45"))
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
  scale_color_manual(values=c("#A6DDF0", "#276FBF","#B4ED50","#2E8B57", "#FFDD33", "#C49E45")) +
  scale_fill_manual(values=c("#A6DDF0", "#276FBF","#61BA5C", "#E2BF3C"))+
  facet_wrap(~k,ncol=1)
p3


#combined_plot <- wrap_plots(p, p2, heights = c(0.15, 1), ncol=1)
#combined_plot

#ggsave("figures/Admixture_plot.pdf", combined_plot, width = 7, height = 8, units="in")
#ggsave("figures/Admixture_plot.png", combined_plot, width = 7, height = 8, units="in")

p_out1 <- ggarrange(p2, p3, ncol = 2, nrow = 1, labels=c("B", "C"), 
                    common.legend = T, legend = "none")

p_out <- ggarrange(p1,p_out1, nrow=2, labels=c("A", "", ""), common.legend=F,
                   heights = c(0.5, 0.5))
p_out
#widths = c(1, 0.75)

ggsave(filename="figures/Fig_01.pdf", p_out, h=135, w=170, units= "mm")
ggsave(filename="figures/Fig_01.png", p_out, h=5, w=7)


