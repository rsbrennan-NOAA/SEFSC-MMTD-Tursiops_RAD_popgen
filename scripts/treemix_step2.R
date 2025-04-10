
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

# 2 migration events. 


#-------------------
##### six pop
setwd("~/projects/Tursiops_RAD_popgen/analysis/pop_structure/treemix/sixpop/")

plot_tree("test_migrations/treemix.4.4")
#plot_tree("bootstrap/sixpop.p.treemix.gz_constree_bootrep_3")


folder <- file.path(path="test_migrations")                     #path to files of TreeMix replicates with different migration edges (m) to test
test.linear = optM(folder, method = "linear", tsv="linear.txt")   #test m: produces a list of which the $out dataframe suggests optimum m based on multiple linear models
plot_optM(test.linear, method = "linear")                         #shows changes in log likelihood for different m, and suggests optimum m as 'change points'

# 4 migration events.
