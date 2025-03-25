

library(vcfR)
library(adegenet)
library(ggpubr)
library(dplyr)
library(poppr)

vcf <- read.vcfR( "analysis/filtered.final_ids_LDthin.vcf.gz", verbose = FALSE )

genl<-vcfR2genlight(vcf)

# add population ids
pops <- read.csv("Tursiops_RADseq_Metadata_new.csv")

#genl@ind.names <- gsub("b", "",genl@ind.names)

ids <- data.frame(IDs = genl@ind.names)
#result <- ids %>%
#  left_join(pops %>% select( Lab.ID, Pop.Structure.Location), by = c("IDs" = "Lab.ID"))

#genl@pop <- as.factor(result$Pop.Structure.Location)

## Full data set
# re-running pca, basically to make sure consistent again with plink, etc. 
pca1 <- glPca(genl,center = T, scale = F, nf = 5)

png("figures/dapc_eigen.png", h=4, w=4, units="in", res=300)
barplot(pca1$eig, col = heat.colors(50), main="PCA Eigenvalues",
        xlim=c(1,25)) # retain first 5 axes, incremental decrease after 2
title(ylab="value", line = 2)
title(xlab="Eigenvalues", line = 1)

dev.off()

#proportion of explained variance by first three axes
a1<-pca1$eig[1]/sum(pca1$eig) # proportion of variation explained by 1st axis
a2<-pca1$eig[2]/sum(pca1$eig) # proportion of variation explained by 2nd axis 
a3<-pca1$eig[3]/sum(pca1$eig) # proportion of variation explained by 3rd axis
pcvar <- data.frame(Axis = c(1:3), Proportion = c(a1,a2,a3))
pcvar


###------------------------------------------------------------------------
###------------------------------------------------------------------------
###------------------------------------------------------------------------
# run dapc
###------------------------------------------------------------------------
###------------------------------------------------------------------------
###------------------------------------------------------------------------

# keep ~200 pcs to start
grp <- find.clusters(genl, max.n.clust=40)

# 4 clusters

# do k-1 for pca
dapc1 <- dapc(genl, grp$grp, n.pca=3, n.da=4)

#temp <- optim.a.score(dapc1)
col.in <- c("#E69F00","#56B4E9", "#009E73", "#CC79A7")
scatter.dapc(dapc1)
scatter.dapc(dapc1, grp=grp$grp)

png(file="../figures/dapc_assignments.png", h=4, w=4, units="in", res=300)
scatter(dapc1, grp=genl@pop, 
        bg="white", pch=c(16,15,18,17), col=col.in,
        cex=2,
        cstar=0,clab=0,
        scree.pca=FALSE,scree.da=FALSE, leg=TRUE, posi.leg="bottomleft")

dev.off()


# Calculate optimal number of PCs 
dapcTemp <- dapc(genl, grp$grp, 
                 n.pca=70, n.da = 20)   # n.da = 4(species) - 1 = 3
ascores <- optim.a.score(dapcTemp, smart = FALSE, n.sim = 50)   
# Optimal number of PCs: 4

dapc4 <- dapc(genl, genl@pop, 
              n.pca=17, n.da = 3)

scatter(dapc4)

pramx <- xvalDapc(tab(genl, NA.method = "mean"), grp$grp)

names(pramx)
pramx[-1]
pramx[1]



######
# compare population assignments:

admix <- read.table("analysis/populations_admixture.txt", header=T)


df <- data.frame(
  Lab.ID = names(grp$grp),
  dapc_population = as.numeric(grp$grp)
)

head(df)

#


table(df$dapc_population)
table(admix$admixture_population)

df$dapc_population[df$dapc_population == 1] <- "Coastal_Atl"
df$dapc_population[df$dapc_population == 2] <- "Coastal_Gulf"
df$dapc_population[df$dapc_population == 3] <- "Intermediate"
df$dapc_population[df$dapc_population == 4] <- "Offshore"

df[which(df$dapc_population != admix$admixture_population),]
admix[which(admix$admixture_population != df$dapc_population),]

# 6 samples that don't agree. they're all highly admixed between population according to Admixture. 
# DAPC calls them all intermediate. which is reasonable honestly. But so are the admixture assigns. 

# admixture scores:
#IDs       Q1       Q2       Q3       Q4          highest_Q
#1 10Tt007 0.487264 0.000010 0.010431 0.502294 Coastal_Gulf
#2 10Tt059 0.455808 0.000010 0.000010 0.544172 Coastal_Gulf
#3 13Tt073 0.358837 0.148496 0.462016 0.030651     Offshore
#4 37Tt023 0.269821 0.000010 0.242777 0.487392 Coastal_Gulf
#5  7Tt312 0.413666 0.032105 0.458211 0.096018     Offshore
#6  8Tt144 0.440412 0.000010 0.463625 0.095952     Offshore


write.table(file="analysis/populations_dapc.txt", df, sep="\t", 
            quote = F, row.names=F)
