#PCA

library(ggplot2)
library(dplyr)
library(SNPRelate)
library(SeqArray)
library(stringr)

pops <- read.csv("Tursiops_RADseq_Metadata_new.csv")

filename = "analysis/filtered.final_ids_LDthin"
filename.gds = paste0(filename, ".gds")
filename.vcf.gz = paste0(filename, ".vcf.gz")
# Convert VCF to GDS
SeqArray::seqVCF2GDS(vcf.fn = filename.vcf.gz, out.fn = filename.gds, storage.option="ZIP_RA")

gdsin = SeqArray::seqOpen(filename.gds)

print(paste0("The number of SAMPLES in data: ", length(c(SeqArray::seqGetData(gdsin, "sample.id")))))
# [1] "The number of SAMPLES in data: 337"
print(paste0("The number of SNPs in data: ",  length(c(SeqArray::seqGetData(gdsin, "variant.id")))))
# [1] "The number of SNPs in data: 4356"
summary(m1 <- SeqArray::seqMissing(gdsin, per.variant=TRUE))
summary(m2 <- SeqArray::seqMissing(gdsin, per.variant=FALSE))

samples <- SeqArray::seqGetData(gdsin, "sample.id")

samples[!samples %in% pops$Lab.ID]

pop.ids <- str_extract(samples, "^[^Tt]+")
pop.ids <- sub("^(FB).*$", "\\1", pop.ids)

SeqArray::seqSummary(gdsin)

SeqArray::seqSetFilter(gdsin)
snp_in <- SeqArray::seqGetData(gdsin, "variant.id")
samp_in <- SeqArray::seqGetData(gdsin, "sample.id")



#--------------------------------------------
# pca

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
dat$Population <- pop.ids

# with ids
d <- ggplot(dat, aes(PC1, PC2, label=Population, color=Population)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('All individuals: PC1, PC2')
d

ggsave("figures/pca_popIds.png", d, h=4, w=5)


# color by gulf or atlantic

pops$region <- ifelse(pops$Long > -81.7, "Atlantic",
                     ifelse(pops$Long > 80, "Gulf", NA))
pops$region <- ifelse(pops$Long < -81.7, "Gulf", pops$region)

df <- merge(dat, pops, by.x="IDs", by.y="Lab.ID")


# with ids
d <- ggplot(df, aes(PC1, PC2, label=Population, color=region)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('All individuals: PC1, PC2')
d

ggsave("figures/pca_popIds_region.png", d, h=4, w=5)



# with shapes
d <- ggplot(df, aes(PC1, PC2, label=Population, fill=region, shape=region)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('PCA: PC1, PC2')+
  scale_shape_manual(values=c(21, 24))
d

ggsave("figures/pca_region.png", d, h=4, w=5)




# color by depth
depth <- read.csv("depths.csv", header=T)

df <- merge(df, depth, by.x="IDs", by.y="id")
df$log_corrected_depth <- log10(df$corrected_depth*-1)

d <- ggplot(df, aes(PC1, PC2, label=Population, fill=corrected_depth, shape=region)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('All individuals: PC1, PC2') +
  scale_fill_distiller(
    palette = "Blues",
    #breaks = c(seq(0, 500, by = 100), 4000)*-1,
    labels = function(x) ifelse(x == -200, " < -200", x),
    name = "Depth",
    limits = c(-200, 0),  
    oob = scales::squish, # anything below -800 will be the yellow
    direction = -1  # reverse direction if needed
  ) +
  scale_shape_manual(values=c(21, 24))

d

ggsave("figures/pca_depth.png", d, h=4, w=5)




# color by distance from shore
depth <- read.csv("analysis/depth_distance.csv", header=T)

df <- merge(df, depth, by.x="IDs", by.y="id")

d <- ggplot(df, aes(PC1, PC2, label=Population, fill=distance_to_shore/1000, shape=region)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('All individuals: PC1, PC2') +
  scale_fill_distiller(
    palette = "Blues",
    #breaks = c(seq(0, 500, by = 100), 4000)*-1,
    labels = function(x) ifelse(x == 200, " > 200", x),
    name = "Distance to \nshore (km)",
    limits = c(0,200),  
    oob = scales::squish,
    direction = 1
  ) +
  scale_shape_manual(values=c(21, 24))
d

ggsave("figures/pca_distance.png", d, h=4, w=5)



#-------------------------------------------------------------------------------
# look at loadings:
snp_ids <- SeqArray::seqGetData(gdsin, "annotation/id")
chromosomes <- sub("_[^_]*$", "", snp_ids)
table(chromosomes)
chromosomes <- sub(".1$", "", chromosomes)      # Remove .1 suffix
chromosome_ids <- sapply(strsplit(chromosomes, "_"), `[`, 2)
chromosomes <- as.numeric(chromosome_ids)

positions <- SeqArray::seqGetData(gdsin, "position")

SnpLoad <- snpgdsPCASNPLoading(pca.out, gdsin)
dim(SnpLoad$snploading)

pc<-1

# make df
plot_data <- data.frame(
  CHR = as.numeric(chromosomes),
  BP = positions,
  P = abs(SnpLoad$snploading[pc,])  # Using absolute values of loadings
)

# Order by chromosome and position
plot_data <- plot_data[order(plot_data$CHR, plot_data$BP),]

# Calculate cumulative position for x-axis
plot_data$cumpos <- NA
lastbase <- 0
nbp <- c()
for(i in unique(plot_data$CHR)){
  nbp[i] <- max(plot_data[plot_data$CHR == i,]$BP)
  plot_data[plot_data$CHR == i,"cumpos"] <- plot_data[plot_data$CHR == i,"BP"] + lastbase
  lastbase <- lastbase + nbp[i]
}


png(filename = paste0("figures/manhattan_plot_pc",pc,".png"), width=8, height=4, units="in", res=200)
par(mar=c(5,5,4,2))

plot(plot_data$cumpos, plot_data$P,
     xaxt="n",
     xlab="Chromosome",
     ylab="PC Loading",
     main=paste("Manhattan Plot of PC", pc, "Loadings"),
     pch=20,
     col=ifelse(plot_data$CHR %% 2 == 0, "navy", "grey45"),
     cex=0.8)

dev.off()



pc<-2

# make df
plot_data <- data.frame(
  CHR = as.numeric(chromosomes),
  BP = positions,
  P = abs(SnpLoad$snploading[pc,])  # Using absolute values of loadings
)

# Order by chromosome and position
plot_data <- plot_data[order(plot_data$CHR, plot_data$BP),]

# Calculate cumulative position for x-axis
plot_data$cumpos <- NA
lastbase <- 0
nbp <- c()
for(i in unique(plot_data$CHR)){
  nbp[i] <- max(plot_data[plot_data$CHR == i,]$BP)
  plot_data[plot_data$CHR == i,"cumpos"] <- plot_data[plot_data$CHR == i,"BP"] + lastbase
  lastbase <- lastbase + nbp[i]
}


png(filename = paste0("figures/manhattan_plot_pc",pc,".png"), width=8, height=4, units="in", res=200)
par(mar=c(5,5,4,2))

plot(plot_data$cumpos, plot_data$P,
     xaxt="n",
     xlab="Chromosome",
     ylab="PC Loading",
     main=paste("Manhattan Plot of PC", pc, "Loadings"),
     pch=20,
     col=ifelse(plot_data$CHR %% 2 == 0, "navy", "grey45"),
     cex=0.8)

dev.off()





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# color by admixture populations

admixPop <- read.table("analysis/populations_admixture.txt", header=T)
pops$region <- ifelse(pops$Long > -81.7, "Atlantic",
                      ifelse(pops$Long > 80, "Gulf", NA))
pops$region <- ifelse(pops$Long < -81.7, "Gulf", pops$region)

df <- merge(dat, pops, by.x="IDs", by.y="Lab.ID")

ndat <- merge(df, admixPop, by.x="IDs", by.y="Lab.ID")
# with shapes
d <- ggplot(ndat, aes(PC1, PC2, label=Population, fill=admixture_population, shape=region)) +
  geom_point(size=3, color="black", aes(fill=admixture_population)) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('PCA: admixture populations')+
  scale_shape_manual(values=c(21,24))+
  guides(fill=guide_legend(override.aes=list(shape=21)))

d

ggsave("figures/pca_admixturePops.png", d, h=4, w=5)

# repeat with dapc populations
admixPop <- read.table("analysis/populations_dapc.txt", header=T)
pops$region <- ifelse(pops$Long > -81.7, "Atlantic",
                      ifelse(pops$Long > 80, "Gulf", NA))
pops$region <- ifelse(pops$Long < -81.7, "Gulf", pops$region)

df <- merge(dat, pops, by.x="IDs", by.y="Lab.ID")

ndat_dapc <- merge(df, admixPop, by.x="IDs", by.y="Lab.ID")
# with shapes
d <- ggplot(ndat_dapc, aes(PC1, PC2, label=Population, fill=dapc_population, shape=region)) +
  geom_point(size=3, color="black", aes(fill=dapc_population)) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('PCA: DAPC populations')+
  scale_shape_manual(values=c(21,24))+
  guides(fill=guide_legend(override.aes=list(shape=21)))

d

ggsave("figures/pca_DAPCPops.png", d, h=4, w=5)

misid <- data.frame(IDs= c("37Tt023","13Tt073","8Tt144","10Tt007","7Tt312","10Tt059"))

misall <- merge(misid, ndat_dapc, by="IDs")

d <- ggplot(ndat_dapc, aes(PC1, PC2, label=Population, fill=dapc_population, shape=region)) +
  geom_point(size=3, color="black", aes(fill=dapc_population)) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('PCA: admixture vs dapc misassignments')+
  scale_shape_manual(values=c(21,24))+
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  geom_point(data=misall, size=6, fill=NA, color="black", shape=21, aes(PC1, PC2, fill=admixture_population, shape=region))

d


ggsave("figures/pca_DAPC_admixture_misAssignments.png", d, h=4, w=5)




# pull out the intermediate offshore cluster:
ndat[(ndat$PC1> 0.025 & ndat$PC1 < 0.1 & ndat$PC2 < 0.05 & ndat$PC2 > -0.01),]
ndat[(ndat$PC1> 0.025 & ndat$PC1 < 0.1 & ndat$PC2 < 0.045 & ndat$PC2 > -0.01),]

ndat$IDs[(ndat$PC1> 0.025 & ndat$PC1 < 0.1 & ndat$PC2 < 0.045 & ndat$PC2 > -0.01)]
# "8Tt114" "8Tt243" "8Tt322" "8Tt363" "8Tt365" "8Tt424" "8Tt428" "8Tt465" "8Tt495" "8Tt513" "8Tt521" "8Tt619" "8Tt654"
ndat_dapc$IDs == ndat$IDs

# table with dapc pops, low confidence indivs, and the intermediate offshore indivs
summary_dat <- data.frame(indiv = ndat$IDs,
           region = ndat$region,
           dapc_population = ndat_dapc$dapc_population,
           admixture_population = ndat$admixture_population,
           offshore_putative_hybrids = (ndat$PC1> 0.025 & ndat$PC1 < 0.1 & ndat$PC2 < 0.045 & ndat$PC2 > -0.01),
           low_confidence_population = ndat$IDs %in% misid$IDs)
head(summary_dat)

# make names better:
summary_dat$dapc_population[summary_dat$dapc_population == "Coastal_Atl"] <- "Coastal_Atlantic" 
summary_dat$admixture_population[summary_dat$admixture_population == "Coastal_Atl"] <- "Coastal_Atlantic" 

summary_dat$fourpop <- summary_dat$dapc_population
summary_dat$sixpop <- paste(summary_dat$region, summary_dat$fourpop, sep="_")

summary_dat2 <- as_tibble(summary_dat) %>%
  mutate(sixpop = str_replace(sixpop, "Atlantic_Coastal_Atlantic", "Coastal_Atlantic") %>%
           str_replace("Atlantic_Intermediate", "Intermediate_Atlantic") %>%
           str_replace("Atlantic_Offshore", "Offshore_Atlantic") %>%
           str_replace("Gulf_Coastal_Gulf", "Coastal_Gulf") %>%
           str_replace("Gulf_Intermediate", "Intermediate_Gulf") %>%
           str_replace("Gulf_Offshore", "Offshore_Gulf"))


summary_dat2 <- as_tibble(summary_dat) %>%
  mutate(across(sixpop, ~str_replace(., "Atlantic_Coastal_Atlantic", "Coastal_Atlantic") %>%
                  str_replace("Atlantic_Intermediate", "Intermediate_Atlantic") %>%
                  str_replace("Atlantic_Offshore", "Offshore_Atlantic") %>%
                  str_replace("Gulf_Coastal_Gulf", "Coastal_Gulf") %>%
                  str_replace("Gulf_Intermediate", "Intermediate_Gulf") %>%
                  str_replace("Gulf_Offshore", "Offshore_Gulf")))

head(as.data.frame(summary_dat2))

sum(summary_dat2$offshore_putative_hybrids == TRUE)

write.table(file="analysis/population_assignments_summary.txt", summary_dat2, sep="\t", 
            quote = F, row.names=F)








#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
# write out files for treemix:
# need to have relatively equal sample sizes

# 1: 4 pops from admixture
# 2: separate gulf and atlantic for intermediate and offshore.

# 1: 4 pop
table(summary_dat2$fourpop)
sum(table(ndat$admixture_population))

write.table(file="analysis/pop_structure/fourpop_all.clust", 
            data.frame(samp1 = summary_dat2$indiv,
                       group = summary_dat2$fourpop),
            sep="\t", quote=F, col.names=F, row.names=F)



# select 70 indivs from Coastal_Atl
# Subset Coastal_Atl IDs, keep all others
coastal_atl_subset_idx <- sample(summary_dat2$indiv[summary_dat2$fourpop == "Coastal_Atlantic"], 70)
coastal_atl_subset <-  summary_dat2[summary_dat2$indiv %in% coastal_atl_subset_idx,]
nrow(coastal_atl_subset)

other_dat <- summary_dat2[summary_dat2$fourpop != "Coastal_Atlantic",]

combodat <- rbind(coastal_atl_subset, other_dat)

# Verify new counts
table(combodat$fourpop)
#Coastal_Atl Coastal_Gulf Intermediate     Offshore 
#70           47           54           74 
sum(table(combodat$fourpop))
# 245

write.table(file="analysis/pop_structure/fourpop.clust", 
            data.frame(samp1 = combodat$indiv,
                       group = combodat$fourpop),
            sep="\t", quote=F, col.names=F, row.names=F)


###----------------------------------------------------------------------------
# 2: separate gulf and atlantic for intermeidate and offshore.

# add region to the pops

table(summary_dat2$sixpop)
#Coastal_Atlantic          Coastal_Gulf Intermediate_Atlantic     Intermediate_Gulf 
#162                    44                    44                    16 
#Offshore_Atlantic         Offshore_Gulf 
#33                    38 

# write the output, 

write.table(file="analysis/pop_structure/sixpop_all.clust", 
            data.frame(samp1 = summary_dat2$indiv,
                       group = summary_dat2$sixpop),
            sep="\t", quote=F, col.names=F, row.names=F)



# for treemix:
# Subset Coastal_Atl IDs, keep all others
coastal_atl_subset_idx <- sample(summary_dat2$indiv[summary_dat2$fourpop == "Coastal_Atlantic"], 44)
coastal_atl_subset <-  summary_dat2[summary_dat2$indiv %in% coastal_atl_subset_idx,]


other_dat <- summary_dat2[summary_dat2$sixpop != "Coastal_Atlantic",]

combodat <- rbind(coastal_atl_subset, other_dat)

# Verify new counts
table(combodat$sixpop)
#Coastal_Atlantic          Coastal_Gulf Intermediate_Atlantic     Intermediate_Gulf 
#44                    44                    44                    16 
#Offshore_Atlantic         Offshore_Gulf 
#33                    38 
sum(table(combodat$sixpop))
# 210

write.table(file="analysis/pop_structure/sixpop.clust", 
            data.frame(samp1 = combodat$indiv,
                       group = combodat$sixpop),
            sep="\t", quote=F, col.names=F, row.names=F)







