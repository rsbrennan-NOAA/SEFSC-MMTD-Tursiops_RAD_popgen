

library(vcfR)
library(adegenet)
library(ggpubr)
library(dplyr)
library(poppr)

vcf <- read.vcfR( "analysis/variants/filtered.final_ids_LDthin_noX.vcf.gz", verbose = FALSE )

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

# keep 200 pcs to start
grp <- find.clusters(genl, max.n.clust=30,n.pca=200)

# 4 clusters

# do k-1 for pca
dapc1 <- dapc(genl, grp$grp, n.pca=3, n.da=4)

#col.in <- c(c("#A6DDF0", "#276FBF","#B4ED50","#2E8B57", "#FFDD33", "#C49E45"))
col.in <- c(c("#A6DDF0", "#276FBF","#FFDD33","#B4ED50"))
#scatter.dapc(dapc1)
scatter.dapc(dapc1, grp=grp$grp)

pdf(file="figures/dapc_assignments.pdf", h=4, w=4)
scatter(dapc1, grp=grp$grp, 
        bg="white", pch=c(16,16,17,15), col=col.in,
        cex=2,
        cstar=0,clab=0,
        scree.pca=FALSE,scree.da=FALSE, leg=TRUE, posi.leg="topleft")

dev.off()


# Calculate optimal number of PCs 
dapcTemp <- dapc(genl, grp$grp, 
                 n.pca=70, n.da = 20)   # n.da = 4(species) - 1 = 3
ascores <- optim.a.score(dapcTemp, smart = FALSE, n.sim = 50)   
# Optimal number of PCs: 3

dapc4 <- dapc(genl, genl@pop, 
              n.pca=17, n.da = 3)

scatter(dapc4)

pramx <- xvalDapc(tab(genl, NA.method = "mean"), grp$grp)


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

df$dapc_population[df$dapc_population == 1] <- "Coastal_Atlantic"
df$dapc_population[df$dapc_population == 2] <- "Coastal_Gulf"
df$dapc_population[df$dapc_population == 3] <- "Offshore"
df$dapc_population[df$dapc_population == 4] <- "Intermediate"

df[which(df$dapc_population != admix$admixture_population),]
admix[which(admix$admixture_population != df$dapc_population),]

# 4 samples that don't agree. they're all highly admixed between population according to Admixture. 
# DAPC calls them all intermediate. which is reasonable honestly. But so are the admixture assigns. 

write.table(file="analysis/populations_dapc.txt", df, sep="\t", 
            quote = F, row.names=F)
