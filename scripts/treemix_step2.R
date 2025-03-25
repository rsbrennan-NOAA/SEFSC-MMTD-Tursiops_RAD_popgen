
source("analysis/pop_structure/treemix/plotting_funcs.R")
setwd("~/projects/Tursiops_RAD_popgen/analysis/pop_structure/treemix/")
par(mfrow=c(1,1),   
    mar=c(2, 2, 1, 1),
    oma=c(0, 0, 0, 0))

pdf("../../../figures/TreeMix_noMigration_fourpop.pdf", h=4, w=5)                                          
plot_tree("fourpop_nomigration")
dev.off()

pdf("../../../figures/TreeMix_noMigration_sixpop.pdf", h=4, w=5)                                          
plot_tree("sixpop_nomigration")
dev.off()

plot_resid(stem="fourpop_noroot",pop_order="fourpop.order")
plot_resid(stem="sixpop_noroot",pop_order="sixpop.order")


# load devtools
library(dplyr)
library(data.table)
library(BITEV2)
library(OptM)
library(plyr)
source("~/projects/Tursiops_RAD_popgen/scripts/TreeMix_functions.R") #path to required functions for this analysis

########################################################
############# (A) Test migration events ################
########################################################
setwd("~/projects/Tursiops_RAD_popgen/analysis/pop_structure/treemix/fourpop/")

plot_tree("test_migrations/treemix.2.9")

folder <- file.path(path="test_migrations")                     #path to files of TreeMix replicates with different migration edges (m) to test
test.linear = optM(folder, method = "linear", tsv="linear.txt")   #test m: produces a list of which the $out dataframe suggests optimum m based on multiple linear models
plot_optM(test.linear, method = "linear")                         #shows changes in log likelihood for different m, and suggests optimum m as 'change points'

# 2 migreation events. 


#-------------------
##### six pop
setwd("~/projects/Tursiops_RAD_popgen/analysis/pop_structure/treemix/sixpop/")

plot_tree("test_migrations/treemix.4.4")
#plot_tree("bootstrap/sixpop.p.treemix.gz_constree_bootrep_3")


folder <- file.path(path="test_migrations")                     #path to files of TreeMix replicates with different migration edges (m) to test
test.linear = optM(folder, method = "linear", tsv="linear.txt")   #test m: produces a list of which the $out dataframe suggests optimum m based on multiple linear models
plot_optM(test.linear, method = "linear")                         #shows changes in log likelihood for different m, and suggests optimum m as 'change points'

# 4 migration events.









######################################################## --------------------------
################## (B) Plot Tree #######################
########################################################

## 1. From the final runs, compare tree likelihoods, 
   # select tree with highest likelihood, 
   # remove duplicates and retain tree(s) with unique topology. 
#Adapted from R functions written by Zecca, Labra and Grassi (2019).

setwd("~/projects/Tursiops_RAD_popgen/analysis/pop_structure/treemix/fourpop/bootstrap")
#folder with all TreeMix outputs from the final runs
source("~/projects/Tursiops_RAD_popgen/scripts/TreeMix_functions.R") #path to required functions for this analysis

maxLL("fourpop_treemix_bootrep_", nt=100, uel=FALSE)
#cfTrees("fourpop_treemix_bootrep_", nt=100, p=1, m='PH85')                        

## 2. Now plot and save unique tree with highest likelihood:

pdf("../../../../../figures/TreeMix_output_fourpop.pdf", h=4, w=5)                                          
treemix.bootstrap(in.file="fourpop_treemix_bootrep_7", #stem of TreeMix files (as above) + number of run with highest likelihood and unique topology
                  out.file = "tmp", 
                  phylip.file = "../fourpop_outtree.newick", 
                  #consensus tree in newick format (from the bootstrap procedure generated with PHYLIP)    
                  nboot = 100, #nboot is the number of bootstraps used
                  fill = TRUE,  
                  #pop.color.file = "col.txt",#specify colours with a tab delimited pop.color.file - first column is pop name, second the colour
                  boot.legend.location = "topleft")
dev.off()

treemix.drift(in.file = "fourpop_treemix_bootrep_7")
              #pairwise matrix for drift estimates with specified order of pops 

write.table(file="fourpop.order", data.frame(pops= c("Aduncus", "Off", "CAtl", "Int", "CGulf")),
            quote = F, row.names=F, col.names = F)

plot_resid("fourpop_treemix_bootrep_7",                                        #pairwise matrix for residuals
           pop_order = "fourpop.order") 

  #  title("Residuals")                         
#dev.off()




#-----------------------------
# six pop
setwd("~/projects/Tursiops_RAD_popgen/analysis/pop_structure/treemix/sixpop/bootstrap")
#folder with all TreeMix outputs from the final runs
source("~/projects/Tursiops_RAD_popgen/scripts/TreeMix_functions.R") #path to required functions for this analysis


result <- maxLL("sixpop_treemix_bootrep_", nt=100, uel=FALSE)
cfTrees("sixpop_treemix_bootrep_", nt=100, p=1, m='PH85')                        

## 2. Now plot and save unique tree with highest likelihood:

pdf("../../../../../figures/TreeMix_output_sixpop_maxLL.pdf", h=4, w=5)                                          
treemix.bootstrap(in.file="sixpop_treemix_bootrep_5", #stem of TreeMix files (as above) + number of run with highest likelihood and unique topology
                  out.file = "tmp", 
                  phylip.file = "../sixpop_outtree.newick", 
                  #consensus tree in newick format (from the bootstrap procedure generated with PHYLIP)    
                  nboot = 100, #nboot is the number of bootstraps used
                  fill = FALSE,  
                  #pop.color.file = "col.txt",#specify colours with a tab delimited pop.color.file - first column is pop name, second the colour
                  boot.legend.location = "topright")
dev.off()

for(i in 1:5){

pdf(paste0("../../../../../figures/TreeMix_output_sixpop_LL_",i,".pdf"), h=4, w=5)                                          
plot_tree(paste0("sixpop_treemix_bootrep_",result$top_likelihood_trees$treeN[i]))

dev.off()
}


pdf("../../../../../figures/TreeMix_output_sixpop_top6_LL.pdf", h=5, w=8)                                          
par(mfrow=c(3,3),   
    mar=c(2, 2, 1, 1),
    oma=c(0, 0, 0, 0))
for(i in 1:6){
  
  plot_tree(paste0("sixpop_treemix_bootrep_",result$top_likelihood_trees$treeN[i]))
  
}
dev.off()

png("../../../../../figures/TreeMix_output_sixpop_top6_LL.png", h=5, w=8, units="in", res=200)                                          
par(mfrow=c(3,3),   
    mar=c(2, 2, 1, 1),
    oma=c(0, 0, 0, 0))
for(i in 1:6){
  
  plot_tree(paste0("sixpop_treemix_bootrep_",result$top_likelihood_trees$treeN[i]))
  
}
dev.off()

































































treemix.drift(in.file = "fourpop_treemix_bootrep_5")
#pairwise matrix for drift estimates with specified order of pops 

write.table(file="fourpop.order", data.frame(pops= c("Aduncus", 
                                                     "OGulf", "OAtl",  
                                                     "CGulf", "CAtl", 
                                                     "IntGulf", "IntAtl")),
            quote = F, row.names=F, col.names = F)

plot_resid("sixpop_treemix_bootrep_42",                                        #pairwise matrix for residuals
           pop_order = "fourpop.order") 


########################################################
######### (C) Weights, Std. Err and p-values ###########
########################################################

#Output reports mean weight of edge, 
  # the jackknife estimate of the weight and standard error 
  #(averaged over N independent runs), 
#the least significant p-value recovered  over N runs for each migration event. 
  # Set working directory to the final_runs folder.

GetMigrStats(input_stem="fourpop_final__2mig_finalrun_", 
             nt=100) 
    #arguments as above, writes file "MS_and_stats.txt" into current directory

########################################################
#### Migration support and corrected MS (Optional) #####
########################################################

#From bootstrap replicates, few other support statistic might be calculated.
#The MS is the percentage of times each pair of label sets is present among n bootstrap replicates.
#Calculated as: (number of matches / number of bootstrap replicates)*100 from all independent runs in the current working directory.
#For the Extended MS (MSE), the number of counts is corrected for multiple matches to avoid over-counting.
#Based on R funcions written by Zecca, Labra and Grassi, 2019.

#Now set working directory to folder with all bootstrap replicates generated with optimum number of m in Step 3.
setwd("setwd("~/projects/Tursiops_RAD_popgen/analysis/pop_structure/treemix/final_runs_fourpop/bootstrap/")
")

#Copy PairsOfSets.txt into directory
GetMigrSupp(skipL=1)                                              #calculates MS over all bootstrap replicates, writes file "MigrSupp.txt" into current directory

GetMS_MSe(nmigr=2, min_n=2, fixed="To", skipL=1)                  #default input file is "MigrSupp.txt" created with GetMigrSupp(), writes file "MS_MSE.txt" into working directory
#nmigr = number of migrations, fixed = specifies which taxa/species label set of each pair is kept fixed
#fixed = "From" fixes the origin of m; fixed = "To" (default) fixes the destination of the same m 
#min_n = minimum number of taxa/species labels to be included within the unfixed set(s)

#Ouputs table with columns 'From' (subset of species below the origin of migration edges),
#'To' (the subset of species below the destination of migration edges), Migration Support (MS) and corrected MS with respect to bootstraps (MSE).

treemix.fit(in.file="fourpop_final__2mig_finalrun_2", 
            out.file="fourpop_final__2mig_finalrun_2", m.start=2, m.end=2)






#-------------------------------------------------------------------------------
# six pop

setwd("~/projects/Tursiops_RAD_popgen/analysis/pop_structure/treemix/final_runs_sixpop/final/")
#folder with all TreeMix outputs from the final runs
maxLL("sixpop_final_4mig_finalrun_", nt=100) 
#first argument is stem of TreeMix output files, 
# nt = number of runs

#If n of unique trees = 1, continue with step #2, 
# if n > 1, you might want to create a consensus tree. 
#Note that bootstrap and migration values will not be available for consensus tree, 
# thus you could also choose one ML tree 

cfTrees("sixpop_final_4mig_finalrun_", nt=100, p=1, m='PH85')                        
#m is the method to calculate pairwise distances between unique trees (default is Robinson-Foulds distance)
#p (number from 0.5 to 1)-proportion for a clade to be represented in the consensus tree. Default (1) is strict tree, 0.5 for majority-rule
#plots consensus tree and saves "Consensus.newick" into working directory

## 2. Now plot and save unique tree with highest likelihood:

pdf("../../../../../figures/TreeMix_output_fourpop.pdf", h=4, w=5)                                          
treemix.bootstrap(in.file="sixpop_final_4mig_finalrun_25", #stem of TreeMix files (as above) + number of run with highest likelihood and unique topology
                  out.file = "tmp", 
                  phylip.file = "../../sixpop_final__finalconstree.newick", 
                  #consensus tree in newick format (from the bootstrap procedure generated with PHYLIP)    
                  nboot = 100, #nboot is the number of bootstraps used
                  fill = TRUE, 
                  plotboot=T,
                  #pop.color.file = "col.txt",#specify colours with a tab delimited pop.color.file - first column is pop name, second the colour
                  boot.legend.location = "topright")
dev.off()

tree.text <- readLines("../../sixpop_final__finalconstree.newick")
print(tree.text)
result <- newick.split("../../sixpop_final__finalconstree.newick")


treemix.drift(in.file = "final_runs_fourpop/final/fourpop_final__2mig_finalrun_1")
#pairwise matrix for drift estimates with specified order of pops 

plot_resid("final_runs_fourpop/final/fourpop_final__2mig_finalrun_5",                                        #pairwise matrix for residuals
           pop_order = "fourpop.order") 

#  title("Residuals")                         
#dev.off()


########################################################
######### (C) Weights, Std. Err and p-values ###########
########################################################

#Output reports mean weight of edge, 
# the jackknife estimate of the weight and standard error 
#(averaged over N independent runs), 
#the least significant p-value recovered  over N runs for each migration event. 
# Set working directory to the final_runs folder.

GetMigrStats(input_stem="fourpop_final__2mig_finalrun_", 
             nt=100) 
#arguments as above, writes file "MS_and_stats.txt" into current directory

########################################################
#### Migration support and corrected MS (Optional) #####
########################################################

#From bootstrap replicates, few other support statistic might be calculated.
#The MS is the percentage of times each pair of label sets is present among n bootstrap replicates.
#Calculated as: (number of matches / number of bootstrap replicates)*100 from all independent runs in the current working directory.
#For the Extended MS (MSE), the number of counts is corrected for multiple matches to avoid over-counting.
#Based on R funcions written by Zecca, Labra and Grassi, 2019.

#Now set working directory to folder with all bootstrap replicates generated with optimum number of m in Step 3.
#setwd("final_runs_fourpop/bootstrap/")
setwd("~/projects/Tursiops_RAD_popgen/analysis/pop_structure/treemix/final_runs_fourpop/final/")

#Copy PairsOfSets.txt into directory
GetMigrSupp(skipL=1)      #calculates MS over all bootstrap replicates, writes file "MigrSupp.txt" into current directory

GetMS_MSe(nmigr=2, min_n=2, fixed="To", skipL=1)                  #default input file is "MigrSupp.txt" created with GetMigrSupp(), writes file "MS_MSE.txt" into working directory
#nmigr = number of migrations, fixed = specifies which taxa/species label set of each pair is kept fixed
#fixed = "From" fixes the origin of m; fixed = "To" (default) fixes the destination of the same m 
#min_n = minimum number of taxa/species labels to be included within the unfixed set(s)

#Ouputs table with columns 'From' (subset of species below the origin of migration edges),
#'To' (the subset of species below the destination of migration edges), Migration Support (MS) and corrected MS with respect to bootstraps (MSE).


treemix.fit(in.file="fourpop_final__2mig_finalrun_2", 
            out.file="fourpop_final__2mig_finalrun_2", m.start=2, m.end=2)
















