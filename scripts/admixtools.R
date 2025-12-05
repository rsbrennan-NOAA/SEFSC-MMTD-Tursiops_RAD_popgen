# admixtools

library(admixtools)
library(tidyverse)

# are the intermediate populations from admixture of the coastal and offshore populations?
# if f3 stat is negative, this suggests that A is admixed between a population 
  # related to B and one related to C

#------------------------------------------------------
#reformat the fam files to have the four and 6 pops
dat <- read.table("analysis/pop_structure/admixtools/tursiops_fourpop.fam")

pops <- read.table("analysis/population_assignments_summary.txt", header=T)

# make map
indiv_to_pop <- setNames(pops$fourpop, pops$indiv)
dat$V1 <- indiv_to_pop[dat$V2]

write.table(file="analysis/pop_structure/admixtools/tursiops_fourpop.ind", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)
write.table(file="analysis/pop_structure/admixtools/tursiops_fourpop.fam", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)


dat <- read.table("analysis/pop_structure/admixtools/tursiops_sixpop.fam")


# make map
indiv_to_pop <- setNames(pops$sixpop, pops$indiv)
dat$V1 <- indiv_to_pop[dat$V2]
write.table(file="analysis/pop_structure/admixtools/tursiops_sixpop.ind", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)
write.table(file="analysis/pop_structure/admixtools/tursiops_sixpop.fam", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)

#------------------------
# with aduncus
dat <- read.table("analysis/pop_structure/admixtools/tursiops_aduncus_fourpop.fam")

# make map
indiv_to_pop <- setNames(pops$fourpop, pops$indiv)
dat$V1 <- indiv_to_pop[dat$V2]

# add aduncus:
dat$V1[1:3] <- "Aduncus"


write.table(file="analysis/pop_structure/admixtools/tursiops_aduncus_fourpop.ind", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)
write.table(file="analysis/pop_structure/admixtools/tursiops_aduncus_fourpop.fam", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)

###################
# six pops

dat <- read.table("analysis/pop_structure/admixtools/tursiops_aduncus_sixpop.fam")


# make map
indiv_to_pop <- setNames(pops$sixpop, pops$indiv)
dat$V1 <- indiv_to_pop[dat$V2]
dat$V1[1:3] <- "Aduncus"

write.table(file="analysis/pop_structure/admixtools/tursiops_aduncus_sixpop.ind", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)
write.table(file="analysis/pop_structure/admixtools/tursiops_aduncus_sixpop.fam", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)


#------------------------------------------------------
# the subset pops
dat <- read.table("analysis/pop_structure/admixtools/tursiops_subset_fourpop.fam")


# make map
indiv_to_pop <- setNames(pops$fourpop, pops$indiv)
dat$V1 <- indiv_to_pop[dat$V2]

write.table(file="analysis/pop_structure/admixtools/tursiops_subset_fourpop.ind", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)
write.table(file="analysis/pop_structure/admixtools/tursiops_subset_fourpop.fam", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)


dat <- read.table("analysis/pop_structure/admixtools/tursiops_subset_sixpop.fam")

# make map
indiv_to_pop <- setNames(pops$sixpop, pops$indiv)
dat$V1 <- indiv_to_pop[dat$V2]
write.table(file="analysis/pop_structure/admixtools/tursiops_subset_sixpop.ind", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)
write.table(file="analysis/pop_structure/admixtools/tursiops_subset_sixpop.fam", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)



###----------------------------------------------------------------------------
# run admixtools

# first calc F2 stats, this is needed for downstream steps
prefix = 'analysis/pop_structure/admixtools/tursiops_fourpop'
my_f2_dir = 'analysis/pop_structure/admixtools/f2_tursiops_fourpop'
extract_f2(prefix, my_f2_dir,auto_only = FALSE,
           maxmiss =0.3, overwrite=T)

f2_blocks = f2_from_precomp(my_f2_dir)

dim(f2_blocks)
f2_blocks[,,1]                # f2-statistics of the 1st SNP block
apply(f2_blocks, 1:2, mean)   # average across all blocks

#dat <- read_plink("pop_structure/admixtools/fourpop",verbose = TRUE)
pop1 = "Intermediate"
pop2 = c("Coastal_Gulf", "Coastal_Atlantic")
pop3 = c('Offshore')

fout_foupop <- f3(f2_blocks, pop1, pop2, pop3)
write.table(fout_foupop, file = "analysis/pop_structure/admixtools/f3_stats_fourpop.txt",
            row.names=F, quote=F, sep="\t")


#---------------------------------------------------------
# six pop:

# first calc F2 stats, this is needed for downstream steps
prefix = 'analysis/pop_structure/admixtools/tursiops_sixpop'
my_f2_dir = 'analysis/pop_structure/admixtools/f2_tursiops_sixpop/'
extract_f2(prefix, my_f2_dir,auto_only = FALSE,
           maxmiss =0.3, overwrite=T)

f2_blocks = f2_from_precomp(my_f2_dir)

#dat <- read_plink("pop_structure/admixtools/fourpop",verbose = TRUE)
pop1 = c("Intermediate_Atlantic", "Intermediate_Gulf")
pop2 = c("Coastal_Gulf", "Coastal_Atlantic")
pop3 = c('Offshore_Atlantic', "Offshore_Gulf")

fout_sixpop <- f3(f2_blocks, pop1, pop2, pop3)

write.table(fout_sixpop, file = "analysis/pop_structure/admixtools/f3_stats_sixpop.txt",
            row.names=F, quote=F, sep="\t")


#----------------------------------------------
# four popsubset:

# first calc F2 stats, this is needed for downstream steps
prefix = 'analysis/pop_structure/admixtools/tursiops_subset_fourpop'
my_f2_dir = 'analysis/pop_structure/admixtools/f2_tursiops_subset_fourpop/'
extract_f2(prefix, my_f2_dir,auto_only = FALSE,
           maxmiss =0.3, overwrite=T)

f2_blocks = f2_from_precomp(my_f2_dir)

dim(f2_blocks)
f2_blocks[,,1]                # f2-statistics of the 1st SNP block
apply(f2_blocks, 1:2, mean)   # average across all blocks

pop1 = "Intermediate"
pop2 = c("Coastal_Gulf", "Coastal_Atlantic")
pop3 = c('Offshore')

f3(f2_blocks, pop1, pop2, pop3)
#pop1         pop2             pop3          est       se     z         p
#<chr>        <chr>            <chr>       <dbl>    <dbl> <dbl>     <dbl>
#  1 Intermediate Coastal_Atlantic Offshore  0.00203 0.000641  3.16 0.00156  
#2 Intermediate Coastal_Gulf     Offshore -0.00241 0.000581 -4.16 0.0000325   

#---------------------------------------------------------
# six pop subset

# first calc F2 stats, this is needed for downstream steps
prefix = 'analysis/pop_structure/admixtools/tursiops_subset_sixpop'
my_f2_dir = 'analysis/pop_structure/admixtools/f2_tursiops_subset_sixpop/'
extract_f2(prefix, my_f2_dir,auto_only = FALSE,
           maxmiss =0.3, overwrite=T)

f2_blocks = f2_from_precomp(my_f2_dir)

#dat <- read_plink("pop_structure/admixtools/fourpop",verbose = TRUE)
pop1 = c("Intermediate_Atlantic", "Intermediate_Gulf")
pop2 = c("Coastal_Gulf", "Coastal_Atlantic")
pop3 = c('Offshore_Atlantic', "Offshore_Gulf")

f3(f2_blocks, pop1, pop2, pop3)
#pop1                  pop2         pop3                    est       se      z        p
#<chr>                 <chr>        <chr>                 <dbl>    <dbl>  <dbl>    <dbl>
#1 Intermediate_Atlantic Coastal_Atlantic Offshore_Atlantic  0.00495  0.000777   6.37 1.84e-10
#2 Intermediate_Atlantic Coastal_Atlantic Offshore_Gulf      0.00492  0.000676   7.29 3.11e-13
#3 Intermediate_Atlantic Coastal_Gulf     Offshore_Atlantic -0.000804 0.000692  -1.16 2.45e- 1
#4 Intermediate_Atlantic Coastal_Gulf     Offshore_Gulf      0.000799 0.000602   1.33 1.85e- 1
#5 Intermediate_Gulf     Coastal_Atlantic Offshore_Atlantic -0.000972 0.000727  -1.34 1.81e- 1
#6 Intermediate_Gulf     Coastal_Atlantic Offshore_Gulf     -0.00102  0.000617  -1.66 9.75e- 2
#7 Intermediate_Gulf     Coastal_Gulf     Offshore_Atlantic -0.00699  0.000649 -10.8  4.47e-27
#8 Intermediate_Gulf     Coastal_Gulf     Offshore_Gulf     -0.00541  0.000559  -9.68 3.60e-22


# basic equation to figure out how we get various values
a=0.5
b=0.8
c=0.9 # potentially admixed population
(c-a)*(c-b)












#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#d stats
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------


# write the population tables for dstats:

# add aduncus:

four.out <- rbind(data.frame(indiv = c("SRR5357657", "SRR5357656", "SRR5357655"),
                 pop = c("Outgroup","Outgroup","Outgroup")),
      data.frame(indiv=pops$indiv, 
                 pop=pops$fourpop))
six.out <- rbind(data.frame(indiv = c("SRR5357657", "SRR5357656", "SRR5357655"),
                             pop = c("Outgroup","Outgroup","Outgroup")),
                  data.frame(indiv=pops$indiv, 
                             pop=pops$sixpop))

write.table(four.out, file="analysis/pop_structure/dstats/fourpop_all_dstats.pop", 
            quote=F, col.names=F, row.names = FALSE, sep="\t")
write.table(six.out, file="analysis/pop_structure/dstats/sixpop_all_dstats.pop", 
            quote=F, col.names=F, row.names = FALSE, sep="\t")



































#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#F4 stats
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# https://compvar-workshop.readthedocs.io/en/latest/contents/03_f3stats/f3stats.html
# use aduncus outgroup. 

# F4(A,B;C,D)=⟨(a−b)(c−d)⟩.

prefix = 'analysis/pop_structure/admixtools/tursiops_aduncus_fourpop'
my_f2_dir = 'analysis/pop_structure/admixtools/f2_tursiops_aduncus_fourpop/'
#extract_f2(prefix, my_f2_dir,auto_only = FALSE,
#           maxmiss =0.3, overwrite=T)

f2_blocks = f2_from_precomp(my_f2_dir,afprod = TRUE) #(afprod helps avoid bias)


f4_result <- f4(f2_blocks,
                pop1 = "Aduncus",  # Outgroup
                pop2 = "Coastal_Gulf",  # P1 in your DSuite output
                pop3 = "Coastal_Atl",  # P2 in your DSuite output
                pop4 = "Intermediate")  # P3 in your DSuite output

f4_result <- f4(f2_blocks,
                pop1 = "Aduncus",  # Outgroup
                pop2 = "Coastal_Gulf",  # P1 in  DSuite output
                pop3 = "Offshore",  # P2 in  DSuite output
                pop4 = "Intermediate")  # P3 in  DSuite output

f4_result
print(f4_result)


# This indicates excess allele sharing between Coastal_Gulf and Intermediate,
#For f4(A,B;C,D), it calculates the average of (a-b)*(c-d) 
# so (aduncus - coastal Gulf) * (coastal Atl - Intermediate)
# A negative f4 value means that when the allele frequencies 
   # of Aduncus and Coastal_Gulf differ in one direction, 
   # the frequencies of Coastal_Atl and Intermediate tend to differ in 
   # the opposite direction. This creates two possible explanations
# Excess allele sharing between Aduncus and Coastal_Atl
# Excess allele sharing between Coastal_Gulf and Intermediate
# but aduncus the outgroup, so it isn't this. 

# When the allele frequency difference (aduncus - coastal_gulf) is positive, 
     # the difference (coastal_atl - intermediate) tends to be negative
# OR when (aduncus - coastal_gulf) is negative, (coastal_atl - intermediate) 
    # tends to be positive
# This negative correlation suggests that coastal_gulf and intermediate share alleles 
   # that are not shared between aduncus and coastal_atl.


# Test 1: Testing gene flow between Intermediate and Coastal populations
f4_results <- f4(f2_blocks,
                 pop1 = "Aduncus",
                 pop2 = "Offshore", 
                 pop3 = c("Intermediate", "Coastal_Gulf", "Coastal_Atl"),
                 pop4 = c("Intermediate", "Coastal_Gulf", "Coastal_Atl"))
f4(f2_blocks, pop1, pop2, pop3, pop4)

f4_filtered <- f4_results[f4_results$pop3 != f4_results$pop4,]
f4_intermediate_coastal <- f4(f2_blocks,
                              pop1 = "Aduncus",  # outgroup
                              pop2 = "Coastal_Atl",  # reference population
                              pop3 = "Intermediate",
                              pop4 = "Coastal_Gulf")

# sig gene flow between the Intermediate population and the coastal population.

pop1 = 'Aduncus'
pop2 = 'Offshore'
pop3 = "Intermediate"
pop4 = c("Coastal_Gulf", "Coastal_Atl")
f4(f2_blocks, pop1, pop2, pop3, pop4)



# Test 2: Alternative arrangement to test gene flow between Intermediate and Offshore
pop1 = 'Aduncus'
pop2 = c("Coastal_Gulf", "Coastal_Atl")
pop3 = 'Offshore'
pop4 = "Intermediate"
f4(f2_blocks, pop1, pop2, pop3, pop4)

f4(f2_blocks, pop1, pop2, pop3, pop4, f4mode=FALSE)

f4_fourpops <- f4(f2_blocks)
#

#--------------------
# six pops


prefix = 'analysis/pop_structure/admixtools/tursiops_aduncus_sixpop'
my_f2_dir = 'analysis/pop_structure/admixtools/f2_tursiops_aduncus_sixpop/'
extract_f2(prefix, my_f2_dir,auto_only = FALSE,
           maxmiss =0.3, overwrite=T)

f2_blocks = f2_from_precomp(my_f2_dir,afprod = TRUE) #(afprod helps avoid bias)

# Test 1: Testing gene flow between Intermediate and Coastal populations
pop1 = 'Aduncus'
pop2 = c("Intermediate_Atlantic", "Intermediate_Gulf")
pop3 = c("Coastal_Gulf", "Coastal_Atl")
pop4 = c('Offshore_Gulf', "Offshore_Atlantic")
f4(f2_blocks, pop1, pop2, pop3, pop4)

# sig gene flow between the Intermediate population and the coastal population.

f4_sixpops <- f4(f2_blocks)

f4_sixpops <- f4_sixpops[which(f4_sixpops$pop1 == "Aduncus"),]

as.data.frame(f4_sixpops[order(f4_sixpops$pop2, f4_sixpops$pop3),])
as.data.frame((f4_sixpops))
#










#-----------------------------------------------------------------------
# admixture graph

# first calc F2 stats, this is needed for downstream steps
#prefix = 'analysis/pop_structure/admixtools/tursiops_aduncus_fourpop'
#my_f2_dir = 'analysis/pop_structure/admixtools/f2_tursiops_aduncus_fourpop/'
prefix = 'analysis/pop_structure/admixtools/tursiops_fourpop'
my_f2_dir = 'analysis/pop_structure/admixtools/f2_tursiops_fourpop/'
#extract_f2(prefix, my_f2_dir,auto_only = FALSE,
#           maxmiss =0.3, overwrite=T)

f2_blocks = f2_from_precomp(my_f2_dir)


# https://www.nature.com/articles/s41467-024-47015-y#Fig3
# https://onlinelibrary.wiley.com/doi/full/10.1111/mec.17539
# https://github.com/GenisGE/impalascripts/tree/d2679462e0091ecf529f48836da933097f719add/analyses/gene_flow/qpgraph


#The best-fit graph was estimated using five independent runs of command 
#find_graphs(stop_gen2 = 1000, numgraphs = 100, plusminus_generations = 20). 
#Confidence intervals were estimated using qpgraph_resample_snps with 100 bootstraps. The above methods used statistics (such as f2, f3, and D-statistics) estimated from individual SNPs.
outgraph <- find_graphs(f2_blocks,
                          stop_gen2 = 15, 
                        numgraphs = 10, 
                        plusminus_generations = 5,
                        numadmix=0,
                        max_admix = 1)
  
winner = outgraph %>% slice_min(score)
  
qpgraph(f2_blocks, winner, return_fstats=TRUE)
  
outgraph %>% slice_min(score)
plot_graph(winner$edges[[1]])
  
  write_tsv(winner$edges[[1]], 'fourpop_edges.tsv')
  write_lines(paste0('score\t', winner$score[[1]]),'fourpop_score.tsv')
# Confidence intervals were estimated using qpgraph_resample_snps with 100 bootstraps. 
# The above methods used statistics (such as f2, f3, and D-statistics) 
  # estimated from individual SNPs.

graph_in <- read_tsv("fourpop_edges.tsv")
qpg_results = qpgraph(f2_blocks, graph_in, return_fstats = T)
qpg_results$score
qpg_results$f3
qpg_results$edges
plot_graph(qpg_results$edges)
plotly_graph(qpg_results$edges, fixed=TRUE)
plotly_graph(qpg_results$edges, highlight_unidentifiable = TRUE)



# ------------------------------------
# ------------------------------------
# ------------------------------------
# ------------------------------------
# ------------------------------------

# other approach. this is what I actually use
library(foreach)
library(doParallel)
cl <- makeCluster(10)
#registerDoParallel(cl)

clusterExport(cl, c("find_graphs", "f2_blocks"))
winners <- list()
fits <- list()
nruns <- 1:10
all_pops <- foreach(x = nruns) %dopar% {
  find_graphs(f2_blocks, max_admix=0, stop_gen2 = 30)
}

stopCluster(cl)

winners <- lapply(all_pops, function(x){x %>% slice_min(score, with_ties = FALSE)}) 
fits <- lapply(1:length(winners), function(x){qpgraph(f2_blocks, winners[[x]]$graph[[1]], return_fstats=TRUE)})

winners <- do.call("rbind", winners)

winners$score

# select run with best likelihood
whichbest <- which.min(winners$score)
initial_graph <- fits[[1]]
initial_graph <- fits[[whichbest]]
#$worst_residual
#[1] 8.174426

p <- plot_graph(initial_graph$edges, 
           highlight_unidentifiable=FALSE, hide_weights=T)

ggsave(p, h=4, w=4, file="figures/admixture_graph_m0.pdf")
ggsave(p, h=4, w=4, file="figures/admixture_graph_m0.png")

initial_graph_out <- winners$graph[[whichbest]]




# now figure out if admixture improves the fit
# with admixture:

nruns <- 1:10

# Set up parallel cluster with 10 cores

cl <- makeCluster(10)
registerDoParallel(cl)

clusterExport(cl, c("find_graphs", "f2_blocks"))
winners <- list()
fits <- list()

all_pops <- foreach(x = nruns) %dopar% {
  find_graphs(f2_blocks, max_admix=3, stop_gen2 = 30)
}

stopCluster(cl)

winners <- lapply(all_pops, function(x){x %>% slice_min(score, with_ties = FALSE)}) 
fits <- lapply(1:length(winners), function(x){qpgraph(f2_blocks, winners[[x]]$graph[[1]], return_fstats=TRUE)})

winners <- do.call("rbind", winners)

winners$score

# select run with best likelihood
whichbest <- which.min(winners$score)
bestgraph <- fits[[6]]
bestgraph <- fits[[whichbest]]
#$worst_residual
#[1] 0.004591501
bestgraph_out_3 <- winners$graph[[whichbest]]


p2 <- plot_graph(bestgraph$edges, 
                 highlight_unidentifiable=FALSE, hide_weights=T)

ggsave(p2, h=4, w=4, file="figures/admixture_graph_m3.pdf")
ggsave(p2, h=4, w=4, file="figures/admixture_graph_m3.png")

# ------------------------------------

# 2 admix:

cl <- makeCluster(10)
registerDoParallel(cl)

clusterExport(cl, c("find_graphs", "f2_blocks"))
winners <- list()
fits <- list()

all_pops <- foreach(x = nruns) %dopar% {
  find_graphs(f2_blocks, max_admix=2, stop_gen2 = 30)
}

stopCluster(cl)

winners <- lapply(all_pops, function(x){x %>% slice_min(score, with_ties = FALSE)}) 
fits <- lapply(1:length(winners), function(x){qpgraph(f2_blocks, winners[[x]]$graph[[1]], return_fstats=TRUE)})

winners <- do.call("rbind", winners)

winners$score

# select run with best likelihood
whichbest <- which.min(winners$score)
bestgraph <- fits[[6]]
bestgraph <- fits[[whichbest]]
# 

p3 <- plot_graph(bestgraph$edges, 
                highlight_unidentifiable=FALSE, hide_weights=T)

ggsave(p3, h=4, w=4, file="figures/admixture_graph_m2.pdf")
ggsave(p3, h=4, w=4, file="figures/admixture_graph_m2.png")

bestgraph_out_2 <- winners$graph[[whichbest]]

# ------------------------------------

# 1 admix:

cl <- makeCluster(10)
registerDoParallel(cl)

clusterExport(cl, c("find_graphs", "f2_blocks"))
winners <- list()
fits <- list()

all_pops <- foreach(x = nruns) %dopar% {
  find_graphs(f2_blocks, numadmix =1, stop_gen2 = 30)
}

stopCluster(cl)

winners <- lapply(all_pops, function(x){x %>% slice_min(score, with_ties = FALSE)}) 
fits <- lapply(1:length(winners), function(x){qpgraph(f2_blocks, winners[[x]]$graph[[1]], return_fstats=TRUE)})

winners <- do.call("rbind", winners)

winners$score

# select run with best likelihood
whichbest <- which.min(winners$score)
bestgraph <- fits[[6]]
bestgraph <- fits[[whichbest]]
#$worst_residual
#[1] 0.009292572

p4 <- plot_graph(bestgraph$edges, 
                 highlight_unidentifiable=FALSE, hide_weights=T)

ggsave(p4, h=4, w=4, file="figures/admixture_graph_m1.pdf")
ggsave(p4, h=4, w=4, file="figures/admixture_graph_m1.png")


write_delim(bestgraph$edges, file="analysis/pop_structure/admixtools/graph_out.txt",delim="\t", quote="none")

bestgraph_out_1 <- winners$graph[[whichbest]]




# look at residuals of model
library(dplyr)



#getDZ <- function(pair, df) df[((df[,1] == pair[1] & df[,2]==pair[2]) | (df[,1] == pair[2] & df[,2] == pair[1])), c("diff", "z")]

#print(plotQPgraphRes(bestgraph,
#                     oma=c(1,3,1,0)))

library(ggpubr)

ggsave("figures/admixture_graphs_all.png",ggarrange(p,p4, p3, p2, labels="AUTO"),
       h=6, w=6)

# compare 3 vs 2 vs 1

fits = qpgraph_resample_multi(f2_blocks, list(initial_graph_out, bestgraph_out_1), nboot = 500)
compare_fits(fits[[1]]$score_test, fits[[2]]$score_test)
#$p_emp
#[1] 0.0239521
#$p_emp_nocorr
#[1] 0.02

hist(fits[[2]]$score_test, col="blue", breaks=30,xlim=c(0,200))
hist(fits[[1]]$score_test, breaks=30, add=T)

fits = qpgraph_resample_multi(f2_blocks, list(bestgraph_out_1, bestgraph_out_2), nboot = 500)
compare_fits(fits[[1]]$score_test, fits[[2]]$score_test)
#$p_emp
#[1] 0.3033932
#$p_emp_nocorr
#[1] 0.3
hist(fits[[2]]$score_test, col=scales::alpha("blue", 0.2), breaks=30,xlim=c(0,200))
hist(fits[[1]]$score_test, col=scales::alpha("red", 0.2),breaks=30, add=T)

fits = qpgraph_resample_multi(f2_blocks, list(bestgraph_out_1, bestgraph_out_3), nboot = 500)
compare_fits(fits[[1]]$score_test, fits[[2]]$score_test)
#$p
#[1] 0.9753403
#$p_emp
#[1] 0.9660679
hist(fits[[2]]$score_test, col=scales::alpha("blue", 0.2), breaks=30,xlim=c(0,200), ylim=c(0,20))
hist(fits[[1]]$score_test, breaks=30, add=T, col=scales::alpha("red", 0.2))


plot_graph(bestgraph_out_2, 
           highlight_unidentifiable=FALSE, 
           textsize=4,
           fix=F)


#https://link.springer.com/article/10.1186/s12862-023-02195-x#Fig3


# get bootstrap ci of best model:

fits = qpgraph_resample_snps(f2_blocks, bestgraph_out_1, boot = 100)
fits %>% summarize_fits() %>% plotly_graph(print_highlow = TRUE)

