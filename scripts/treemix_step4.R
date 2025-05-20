
######################################################## --------------------------
################## (B) Plot Tree #######################
########################################################

library(BITEV2)

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
treemix.bootstrap(in.file="fourpop_treemix_bootrep_9", #stem of TreeMix files (as above) + number of run with highest likelihood and unique topology
                  out.file = "tmp", 
                  phylip.file = "../fourpop_outtree.newick", 
                  #consensus tree in newick format (from the bootstrap procedure generated with PHYLIP)    
                  nboot = 100, #nboot is the number of bootstraps used
                  fill = TRUE,  
                  #pop.color.file = "col.txt",#specify colours with a tab delimited pop.color.file - first column is pop name, second the colour
                  boot.legend.location = "topleft")
dev.off()

treemix.drift(in.file = "fourpop_treemix_bootrep_73")
#pairwise matrix for drift estimates with specified order of pops 

write.table(file="fourpop.order", data.frame(pops= c("Aduncus", "Offshore", "CoastAtl", "Intermediate", "CoastGulf")),
            quote = F, row.names=F, col.names = F)

plot_resid("fourpop_treemix_bootrep_7",                                        #pairwise matrix for residuals
           pop_order = "fourpop.order") 

#  title("Residuals")                         
#dev.off()




#-----------------------------
# six pop
setwd("~/projects/Tursiops_RAD_popgen/analysis/pop_structure/treemix/sixpop/bootstrap")
#folder with all TreeMix outputs from the final runs
#source("~/projects/Tursiops_RAD_popgen/scripts/TreeMix_functions.R") #path to required functions for this analysis


result <- maxLL("sixpop_treemix_bootrep_", nt=100, uel=FALSE)
result
cfTrees("sixpop_treemix_bootrep_", nt=100, p=1, m='PH85')                        

## 2. Now plot and save unique tree with highest likelihood:

pdf("../../../../../figures/TreeMix_output_sixpop_maxLL.pdf", h=4, w=5)                                          
treemix.bootstrap(in.file="sixpop_treemix_bootrep_60", #stem of TreeMix files (as above) + number of run with highest likelihood and unique topology
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
