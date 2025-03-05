# admixtools

library(admixtools)
library(tidyverse)

# are the intermediate populations from admixture of the coastal and offshore populations?
# if f3 stat is negative, this suggests that A is admixed between a population 
  # related to B and one related to C

#------------------------------------------------------
#reformat the fam files to have the four and 6 pops
dat <- read.table("analysis/pop_structure/admixtools/tursiops_fourpop.fam")

pops <- read.table("analysis/pop_structure/fourpop_all.clust")
colnames(pops) <- c("indiv", "pop")

# make map
indiv_to_pop <- setNames(pops$pop, pops$indiv)
dat$V1 <- indiv_to_pop[dat$V2]

write.table(file="analysis/pop_structure/admixtools/tursiops_fourpop.ind", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)
write.table(file="analysis/pop_structure/admixtools/tursiops_fourpop.fam", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)


dat <- read.table("analysis/pop_structure/admixtools/tursiops_sixpop.fam")

pops <- read.table("analysis/pop_structure/sixpop_all.clust")
colnames(pops) <- c("indiv", "pop")

# make map
indiv_to_pop <- setNames(pops$pop, pops$indiv)
dat$V1 <- indiv_to_pop[dat$V2]
write.table(file="analysis/pop_structure/admixtools/tursiops_sixpop.ind", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)
write.table(file="analysis/pop_structure/admixtools/tursiops_sixpop.fam", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)

#------------------------
# with aduncus
dat <- read.table("analysis/pop_structure/admixtools/tursiops_aduncus_fourpop.fam")

pops <- read.table("analysis/pop_structure/fourpop_all.clust")
colnames(pops) <- c("indiv", "pop")

# make map
indiv_to_pop <- setNames(pops$pop, pops$indiv)
dat$V1 <- indiv_to_pop[dat$V2]

# add aduncus:
dat$V1[1:3] <- "Aduncus"


write.table(file="analysis/pop_structure/admixtools/tursiops_aduncus_fourpop.ind", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)
write.table(file="analysis/pop_structure/admixtools/tursiops_aduncus_fourpop.fam", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)


dat <- read.table("analysis/pop_structure/admixtools/tursiops_aduncus_sixpop.fam")

pops <- read.table("analysis/pop_structure/sixpop_all.clust")
colnames(pops) <- c("indiv", "pop")

# make map
indiv_to_pop <- setNames(pops$pop, pops$indiv)
dat$V1 <- indiv_to_pop[dat$V2]
dat$V1[1:3] <- "Aduncus"

write.table(file="analysis/pop_structure/admixtools/tursiops_aduncus.ind", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)
write.table(file="analysis/pop_structure/admixtools/tursiops_aduncus.fam", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)


#------------------------------------------------------
# the subset pops
dat <- read.table("analysis/pop_structure/admixtools/tursiops_subset_fourpop.fam")

pops <- read.table("analysis/pop_structure/fourpop_all.clust")
colnames(pops) <- c("indiv", "pop")

# make map
indiv_to_pop <- setNames(pops$pop, pops$indiv)
dat$V1 <- indiv_to_pop[dat$V2]

write.table(file="analysis/pop_structure/admixtools/tursiops_subset_fourpop.ind", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)
write.table(file="analysis/pop_structure/admixtools/tursiops_subset_fourpop.fam", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)


dat <- read.table("analysis/pop_structure/admixtools/tursiops_subset_sixpop.fam")

pops <- read.table("analysis/pop_structure/sixpop_all.clust")
colnames(pops) <- c("indiv", "pop")

# make map
indiv_to_pop <- setNames(pops$pop, pops$indiv)
dat$V1 <- indiv_to_pop[dat$V2]
write.table(file="analysis/pop_structure/admixtools/tursiops_subset_sixpop.ind", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)
write.table(file="analysis/pop_structure/admixtools/tursiops_subset_sixpop.fam", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)





###----------------------------------------------------------------------------
# run admixtools

# first calc F2 stats, this is needed for downstream steps
prefix = 'analysis/pop_structure/admixtools/tursiops_fourpop'
my_f2_dir = 'analysis/pop_structure/admixtools/f2_tursiops_fourpop/'
extract_f2(prefix, my_f2_dir,auto_only = FALSE,
           maxmiss =0.3, overwrite=T)

f2_blocks = f2_from_precomp(my_f2_dir)

dim(f2_blocks)
f2_blocks[,,1]                # f2-statistics of the 1st SNP block
apply(f2_blocks, 1:2, mean)   # average across all blocks



#dat <- read_plink("pop_structure/admixtools/fourpop",verbose = TRUE)
pop1 = "Intermediate"
pop2 = c("Coastal_Gulf", "Coastal_Atl")
pop3 = c('Offshore')

f3(f2_blocks, pop1, pop2, pop3)

#pop1         pop2         pop3           est       se      z             p
#<chr>        <chr>        <chr>        <dbl>    <dbl>  <dbl>         <dbl>
#1 Intermediate Coastal_Atl  Offshore  0.00379  0.000649  5.85  0.00000000503
#2 Intermediate Coastal_Gulf Offshore -0.000549 0.000577 -0.951 0.342    

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
pop2 = c("Coastal_Gulf", "Coastal_Atl")
pop3 = c('Offshore_Atlantic', "Offshore_Gulf")

f3(f2_blocks, pop1, pop2, pop3)
#pop1                  pop2         pop3                    est       se      z        p
#<chr>                 <chr>        <chr>                 <dbl>    <dbl>  <dbl>    <dbl>
#1 Intermediate_Atlantic Coastal_Atl  Offshore_Atlantic  0.00580  0.000762  7.60  2.88e-14
#2 Intermediate_Atlantic Coastal_Atl  Offshore_Gulf      0.00576  0.000679  8.49  2.14e-17
#3 Intermediate_Atlantic Coastal_Gulf Offshore_Atlantic  0.000660 0.000667  0.989 3.23e- 1
#4 Intermediate_Atlantic Coastal_Gulf Offshore_Gulf      0.00200  0.000596  3.36  7.90e- 4
#5 Intermediate_Gulf     Coastal_Atl  Offshore_Atlantic  0.000401 0.000749  0.535 5.92e- 1
#6 Intermediate_Gulf     Coastal_Atl  Offshore_Gulf      0.000283 0.000653  0.433 6.65e- 1
#7 Intermediate_Gulf     Coastal_Gulf Offshore_Atlantic -0.00448  0.000695 -6.44  1.16e-10
#8 Intermediate_Gulf     Coastal_Gulf Offshore_Gulf     -0.00322  0.000609 -5.29  1.25e- 7

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
pop2 = c("Coastal_Gulf", "Coastal_Atl")
pop3 = c('Offshore')

f3(f2_blocks, pop1, pop2, pop3)
#pop1         pop2         pop3           est       se      z            p
#<chr>        <chr>        <chr>        <dbl>    <dbl>  <dbl>        <dbl>
#1 Intermediate Coastal_Atl  Offshore  0.00361  0.000662  5.45  0.0000000499
#2 Intermediate Coastal_Gulf Offshore -0.000549 0.000577 -0.951 0.342       

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
pop2 = c("Coastal_Gulf", "Coastal_Atl")
pop3 = c('Offshore_Atlantic', "Offshore_Gulf")

f3(f2_blocks, pop1, pop2, pop3)
#pop1                  pop2         pop3                    est       se      z        p
#<chr>                 <chr>        <chr>                 <dbl>    <dbl>  <dbl>    <dbl>
#1 Intermediate_Atlantic Coastal_Atl  Offshore_Atlantic  0.00617  0.000781  7.90  2.90e-15
#2 Intermediate_Atlantic Coastal_Atl  Offshore_Gulf      0.00601  0.000695  8.65  5.37e-18
#3 Intermediate_Atlantic Coastal_Gulf Offshore_Atlantic  0.000660 0.000667  0.989 3.23e- 1
#4 Intermediate_Atlantic Coastal_Gulf Offshore_Gulf      0.00200  0.000596  3.36  7.90e- 4
#5 Intermediate_Gulf     Coastal_Atl  Offshore_Atlantic  0.000753 0.000777  0.968 3.33e- 1
#6 Intermediate_Gulf     Coastal_Atl  Offshore_Gulf      0.000511 0.000675  0.756 4.49e- 1
#7 Intermediate_Gulf     Coastal_Gulf Offshore_Atlantic -0.00448  0.000695 -6.44  1.16e-10
#8 Intermediate_Gulf     Coastal_Gulf Offshore_Gulf     -0.00322  0.000609 -5.29  1.25e- 7










#-----------------------------------------------------------------------
# admixture graph

# first calc F2 stats, this is needed for downstream steps
prefix = 'analysis/pop_structure/admixtools/tursiops_sixpop'
my_f2_dir = 'analysis/pop_structure/admixtools/f2_tursiops_sixpop/'
#extract_f2(prefix, my_f2_dir,auto_only = FALSE,
#           maxmiss =0.3, overwrite=T)

f2_blocks = f2_from_precomp(my_f2_dir)

#The best-fit graph was estimated using five independent runs of command 
#find_graphs(stop_gen2 = 1000, numgraphs = 100, plusminus_generations = 20). 
#Confidence intervals were estimated using qpgraph_resample_snps with 100 bootstraps. The above methods used statistics (such as f2, f3, and D-statistics) estimated from individual SNPs.
  outgraph <- find_graphs(f2_blocks,stop_gen2 = 1000, 
                        numgraphs = 100, plusminus_generations = 20)
  outgraph %>% slice_min(score)
  winner = outgraph %>% slice_min(score)
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
plotly_graph(qpg_results$edges)



https://github.com/ekirving/qpbrute

https://link.springer.com/article/10.1186/s12862-023-02195-x#Fig3
