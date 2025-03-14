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

###################
# six pops

dat <- read.table("analysis/pop_structure/admixtools/tursiops_aduncus_sixpop.fam")

pops <- read.table("analysis/pop_structure/sixpop_all.clust")
colnames(pops) <- c("indiv", "pop")

# make map
indiv_to_pop <- setNames(pops$pop, pops$indiv)
dat$V1 <- indiv_to_pop[dat$V2]
dat$V1[1:3] <- "Aduncus"

write.table(file="analysis/pop_structure/admixtools/tursiops_aduncus_sixpop.ind", dat, 
            col.names=F, quote=FALSE, sep="\t", row.names=F)
write.table(file="analysis/pop_structure/admixtools/tursiops_aduncus_sixpop.fam", dat, 
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

fout_foupop <- f3(f2_blocks, pop1, pop2, pop3)
write.table(fout_foupop, file = "analysis/pop_structure/admixtools/f3_stats_fourpop.txt",
            row.names=F, quote=F, sep="\t")


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

fout_sixpop <- f3(f2_blocks, pop1, pop2, pop3)
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

dstat_four <- read.table("analysis/pop_structure/dstats/fourpop_all_BBAA.txt", header=T)
dstat_four

dstat_four$p.value_multTesting <- p.adjust(dstat_four$p.value,method="bonferroni")

write.table(dstat_four,"analysis/pop_structure/dstats/dstats_fourpop.txt",
            sep="\t", quote=F, row.names=F)
# remember, always arranged to be positive. 
# p1 and p2 can be flipped, would make negative.
# as is, indicates gene flow between p2 and p3. 
# so 
dstat_six <- read.table("analysis/pop_structure/dstats/sixpop_all_BBAA.txt", header=T)
dstat_six$p.value_multTesting <- p.adjust(dstat_six$p.value,method="bonferroni")

dstat_six[dstat_six$p.value_multTesting < 0.05,]
write.table(dstat_six,"analysis/pop_structure/dstats/dstats_sixpop.txt",
            sep="\t", quote=F, row.names=F)
#






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
                pop2 = "Coastal_Gulf",  # P1 in your DSuite output
                pop3 = "Offshore",  # P2 in your DSuite output
                pop4 = "Intermediate")  # P3 in your DSuite output

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
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# admixr
library(admixr)

snp_data <- eigenstrat(download_data())

result <- d(
  W = c("French", "Sardinian"), X = "Yoruba", Y = "Vindija", Z = "Chimp",
  data = snp_data
)

result
#













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
