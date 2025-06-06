library(tidyverse)

dat <- read.csv("../Tursiops-NC-PopulationAssignment-RAD/analysis/out.imiss", 
                header=T, sep="\t")

hist(dat$F_MISS, breaks=40)

dat[grep("2Tt534",dat$INDV),]
tout <- table(gsub("b","",dat$INDV))

tout[order(tout)]

dat$INDV
dat[which(dat$F_MISS > 0.75),]
nrow(dat[which(dat$F_MISS > 0.75),])
# 26
# read in metadata. remove the stranded animals
datMeta <- read.csv("../Tursiops-NC-PopulationAssignment-RAD/Tursiops_RADseq_Metadata-original.csv")
stranded <- datMeta$Lab.ID[datMeta$Source == "stranding"]
# keep stranded in. 
length(stranded)# 11 indivs

# and add 7Tt252-rep, which is a replicate and of worse quality than the first one. 
# "8Tt449b-rep" "8Tt449b", but the rep is low cov and gets dropped. 
#rmindiv <- data.frame(d=unique(c(stranded, dat$INDV[which(dat$F_MISS > 0.75)], "7Tt252-rep")))
rmindiv <- data.frame(d=unique(c(dat$INDV[which(dat$F_MISS > 0.75)], "7Tt252-rep")))
table(c(stranded, dat$INDV[which(dat$F_MISS > 0.75)], "7Tt252-rep"))
nrow(rmindiv)
# 27

# three strandings are removed for data quality: 4Tt883 9Tt139 9Tt143 
# 2Tt534 is duplicated, but b gets removed for quality. 


write.table(file="scripts/rm_missing.txt",rmindiv, col.names=F, 
        row.names=F, quote=F)



###### post filtering

dat <- read.csv("analysis/filtered.6.depth.ldepth.mean", 
                header=T, sep="\t")
nrow(dat)
hist(dat$MEAN_DEPTH, breaks=40)
hist(dat$MEAN_DEPTH, breaks=500, xlim=c(0, 200))

mean(dat$MEAN_DEPTH)
# 46.2
mean(dat$MEAN_DEPTH)*3
# *3 = 139

quantile(dat$MEAN_DEPTH, probs=.975)
#136
sum(dat$MEAN_DEPTH > quantile(dat$MEAN_DEPTH, probs=.975))
#228
sum(dat$MEAN_DEPTH > quantile(dat$MEAN_DEPTH, probs=.975))
dat[which(dat$MEAN_DEPTH > quantile(dat$MEAN_DEPTH, probs=.975)),]

# HDplot
library(vcfR)

vcfInput<-read.vcfR("analysis/filtered.6.vcf.gz")
vcfInput

gt_modified <- vcfInput@gt

# Replace .:.:.:.:.:.:.:. with NA
gt_modified[gt_modified == ".:.:.:.:.:.:.:."] <- NA

# Create a new vcfR object with the modified data
vcfInput_modified <- vcfInput
vcfInput_modified@gt <- gt_modified

source("scripts/HDplot.R")

HDplotResults<-HDplot(vcfInput_modified)

head(HDplotResults)
mean(HDplotResults$H)

hist(HDplotResults$num_hets)

HDplotResults %>% ggplot()+geom_point(aes(x=H,y=abs(D)), alpha=0.05) + ylim(0,50)

#plot H and ratio
HDplotResults %>% ggplot()+geom_point(aes(x=H,y=ratio))


# sites beyond the distribution are likely paralogs
## D- > than 5 or 6 for sure, but maybe even 5?
## H - > than about 0.5


sum((HDplotResults$H > 0.5))
#33
sum((abs(HDplotResults$D) > 7), na.rm=T)
#1113

# positions to exclude:
datexclude <- HDplotResults[which(HDplotResults$H > 0.5 | abs(HDplotResults$D) > 7),]
posexclude <- datexclude[,1:2]
nrow(posexclude)
#1120
write.table(posexclude, file="scripts/HD_exclude.txt",
            quote=F, col.names = FALSE, row.names=FALSE,
            sep="\t")





#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# no maf filtee


# HDplot
library(vcfR)

vcfInput<-read.vcfR("analysis/variants/filtered.6.noMAF.vcf.gz")
vcfInput

gt_modified <- vcfInput@gt

# Replace .:.:.:.:.:.:.:. with NA
gt_modified[gt_modified == ".:.:.:.:.:.:.:."] <- NA

# Create a new vcfR object with the modified data
vcfInput_modified <- vcfInput
vcfInput_modified@gt <- gt_modified

source("scripts/HDplot.R")

HDplotResults<-HDplot(vcfInput_modified)

head(HDplotResults)
mean(HDplotResults$H)

hist(HDplotResults$num_hets)

HDplotResults %>% ggplot()+geom_point(aes(x=H,y=abs(D)), alpha=0.05) + ylim(0,20)

#plot H and ratio
HDplotResults %>% ggplot()+geom_point(aes(x=H,y=ratio))


# sites beyond the distribution are likely paralogs
## D- > than 5 or 6 for sure, but maybe even 5?
## H - > than about 0.5


sum((HDplotResults$H > 0.5))
#33
sum((abs(HDplotResults$D) > 7), na.rm=T)
#1113

# positions to exclude:
datexclude <- HDplotResults[which(HDplotResults$H > 0.5 | abs(HDplotResults$D) > 7),]
posexclude <- datexclude[,1:2]
nrow(posexclude)
#1120
write.table(posexclude, file="scripts/HD_exclude_noMAF.txt",
            quote=F, col.names = FALSE, row.names=FALSE,
            sep="\t")



