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

pc<-2
hist(sort(abs(SnpLoad$snploading[pc,]),decreasing=T,index.return=T)[[1]],breaks = 30,main="PC loadings: PC1",)


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


#png("manhattan_plot_pc2.png", width=1200, height=600)
par(mar=c(5,5,4,2))

plot(plot_data$cumpos, plot_data$P,
     xaxt="n",
     xlab="Chromosome",
     ylab=expression(-log[10](abs("PC Loading"))),
     main=paste("Manhattan Plot of PC", pc, "Loadings"),
     pch=20,
     col=ifelse(plot_data$CHR %% 2 == 0, "navy", "grey45"),
     cex=0.8)

#dev.off()