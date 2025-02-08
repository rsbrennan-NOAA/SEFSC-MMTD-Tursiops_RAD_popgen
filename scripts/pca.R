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


dat$region <- ifelse(dat$Long < -81.7, "Gulf", dat$region)
pops$Lat
pops$Long

# color by depth
depth <- read.csv("analysis/depth_distance.csv", header=T)

df <- merge(df, depth, by.x="IDs", by.y="id")

d <- ggplot(df, aes(PC1, PC2, label=Population, fill=distance_to_shore/1000, shape=region)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('All individuals: PC1, PC2') +
  scale_fill_gradient(low = "blue", high = "red", na.value = NA ,
                      breaks = seq(0, 300, by = 50),  # Now in km instead of meters
                      labels = scales::comma,  # This will add commas for thousands
                      name = "Distance to \nshore (km)"
      ) +
                      scale_shape_manual(values=c(21, 24))
d

ggsave("figures/pca_depth.png", d, h=4, w=5)

