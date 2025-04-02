#bgc

library(bgchm)

library(triangulaR)
library(vcfR)

vcf <- read.vcfR( "analysis/filtered.final_ids.vcf.gz", verbose = FALSE )
vcf

pops1 <- read.table("analysis/population_assignments_summary.txt", header=T)

pops <- data.frame(id = pops1$indiv,
                   pop = pops1$fourpop)

# identify ancestry informative markers (AIMS)
# Create a new vcfR object composed only of sites above the given allele frequency difference threshold
diff1 <- alleleFreqDiff(vcfR = vcf, 
                        pm = pops, 
                        p1 = "Coastal_Gulf", 
                        p2 = "Offshore", 
                        difference = 0.7)

genos <- t(extract.gt(diff1))
genos <- gsub("0/0", 0,genos)
genos <- gsub("0/1", 1,genos)
genos <- gsub("1/1", 2,genos)

genos <- apply(genos, 2, as.numeric)
row.names(genos) <- colnames(extract.gt(diff1))
head(genos)

Coastal_Gulf <- pops$id[pops$pop == "Coastal_Gulf"]
Offshore <- pops$id[pops$pop == "Offshore"]
Intermediate <- pops$id[pops$pop == "Intermediate"]

GenP0 <- genos[row.names(genos) %in% Coastal_Gulf,]
GenP1 <- genos[row.names(genos) %in% Offshore,]
GenHybrids <- genos[row.names(genos) %in% Intermediate,]


#


# 0 (homozygote), 1 (heterozygote) and 2 (alternative homozygote).
# 1 row per indiv, 1 col per genotype
# set missing genotypes as NA
## estimate parental allele frequencies, uses analytical solution 
p_out<-est_p(G0=GenP0,
             G1=GenP1,
             model="genotype",ploidy="diploid",
             HMC=FALSE)



## estimate hybrid indexes, uses default HMC settings
## and uses point estimates (posterior medians) of allele frequencies
plin <- is.na(GenHybrids)  
plin[1:5,1:3]
plin <- matrix(ifelse(is.na(GenHybrids), 0, 2), 
               nrow=nrow(GenHybrids), ncol=ncol(GenHybrids))

plin[1:5,1:3]

h_out<-est_hi(Gx=GenHybrids,
              p0=p_out$p0[,1],
              p1=p_out$p1[,1],
              model="genotype",
              ploidy="mixed",
              pldat=plin)

## plot hybrid index estimates with 90% equal-tail probability intervals
## sorted by hybrid index, just a nice way to visualize that in this example we have
## few hybrids with intermediate hybrid indexes
plot(sort(h_out$hi[,1]),ylim=c(0,1),pch=19,xlab="Individual (sorted by HI)",ylab="Hybrid index (HI)")
segments(1:length(h_out$hi[,1]),
         h_out$hi[order(h_out$hi[,1]),3],
         1:60,
         h_out$hi[order(h_out$hi[,1]),4])


## fit a hierarchical genomic cline model for all 51 loci using the estimated
## hybrid indexes and parental allele frequencies (point estimates)
## use 4000 iterations and 2000 warmup to make sure we get a nice effective sample size
gc_out<-est_genocl(Gx=GenHybrids,p0=p_out$p0[,1],
                   p1=p_out$p1[,1],H=h_out$hi[,1],
                   model="genotype",ploidy="mixed",
                   pldat=plin,hier=TRUE,n_iters=4000)

## how variable is introgression among loci? Lets look at the cline SDs
## these are related to the degree of coupling among loci overall
gc_out$SDc
gc_out$SDv

## examine a plot of the joint posterior distribution for the SDs
pp_plot(objs=gc_out,param1="sdv",param2="sdc",probs=c(0.5,0.75,0.95),colors="black",addPoints=TRUE,palpha=0.1,pdf=FALSE,pch=19)

## impose sum-to-zero constraint on log/logit scale
## not totally necessary, but this is mostly a good idea
sz_out<-sum2zero(hmc=gc_out$gencline_hmc,transform=TRUE,ci=0.90)

## plot genomic clines for the 51 loci, first without the sum-to-zero constraint
## then with it... these differ more for some data sets than others
gencline_plot(center=gc_out$center[,1],v=gc_out$gradient,pdf=FALSE)
gencline_plot(center=sz_out$center[,1],v=sz_out$gradient,pdf=FALSE)

## summarize loci with credible deviations from genome-average gradients, here the focus is
## specifically on steep clines indicative of loci introgressing less than the average
which(sz_out$gradient[,2] > 1) ## index for loci with credibly steep clines
sum(sz_out$gradient[,2] > 1) ## number of loci with credibly steep clines

## last, lets look at interspecific ancestry for the same data set, this can
## be especially informative about the types of hybrids present
q_out<-est_Q(Gx=GenHybrids,p0=p_out$p0[,1],p1=p_out$p1[,1],model="genotype",ploidy="diploid")

## plot the results
tri_plot(hi=q_out$hi[,1],Q10=q_out$Q10[,1],pdf=FALSE,pch=19)
## note that some individuals appear to be likely backcrosses (close to the outer lines of the triangles)
## but the individals with intermediate hybrid indexes are clearly not F1s but rather late generation hybrids


