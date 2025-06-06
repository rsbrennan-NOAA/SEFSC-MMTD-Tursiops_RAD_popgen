
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(patchwork)

# fis,  obs and exp het, from stacks


fourpop <- read.csv("analysis/diversity/populations_fourpops/filtered.final.p.sumstats_summary.tsv", 
                    header=T, sep="\t", nrows=4, skip=1)
sixpop <- read.csv("analysis/diversity/populations_sixpops/filtered.final.p.sumstats_summary.tsv", 
                    header=T, sep="\t", nrows=6, skip=1)
colnames(fourpop)[1] <- "population"
colnames(sixpop)[1] <- "population"

fourpop_nohyb <- read.csv("analysis/diversity/populations_noHybs_fourpops/filtered.final_nohybrids.p.sumstats_summary.tsv", 
                    header=T, sep="\t", nrows=4, skip=1)
sixpop_nohyb <- read.csv("analysis/diversity/populations_noHybs_sixpops/filtered.final_nohybrids.p.sumstats_summary.tsv", 
                   header=T, sep="\t", nrows=6, skip=1)
colnames(fourpop)[1] <- "population"
colnames(sixpop)[1] <- "population"
colnames(fourpop_nohyb)[1] <- "population"
colnames(sixpop_nohyb)[1] <- "population"

head(dat)

obshet_plot_fourpop <- ggplot(fourpop, 
                      aes(x = population, 
                          y = Obs_Het, 
                          fill = population, shape=population)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = Obs_Het-StdErr.2, ymax = Obs_Het+StdErr.2), width = 0.2) +
  theme_classic(base_size = 12) +
  ylab("Observed\nheterozygosity") +
  scale_color_manual(values=c("#A6DDF0", "#276FBF","#61BA5C", "#E2BF3C"))+
  scale_fill_manual(values=c("#A6DDF0", "#276FBF","#61BA5C", "#E2BF3C")) +
  scale_shape_manual(values=c(21,21,22,24))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4)))+
  ylim(0, 0.3)

exphet_plot_fourpop <- ggplot(fourpop, 
                      aes(x = population, 
                          y = Exp_Het, 
                          fill = population, shape=population)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = Exp_Het-StdErr.3, ymax = Exp_Het+StdErr.3), width = 0.2) +
  scale_size_manual(values = c("All" = 8, "Population" = 4), guide = "none") +
  theme_classic(base_size = 12) +
  ylab("Expected\nheterozygosity") +
  scale_color_manual(values=c("#A6DDF0", "#276FBF","#61BA5C", "#E2BF3C"))+
  scale_fill_manual(values=c("#A6DDF0", "#276FBF","#61BA5C", "#E2BF3C")) +
  scale_shape_manual(values=c(21,21,22,24))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) +
  ylim(0, 0.3)



#FIS

fis_plot_fourpop <- ggplot(fourpop, 
                   aes(x = population, 
                       y = Fis, 
                       fill = population, shape=population)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = Fis-StdErr.7, ymax = Fis+StdErr.7), width = 0.2) +
  theme_classic(base_size = 12) +
  ylab("Inbreeding\ncoefficient") +
  scale_color_manual(values=c("#A6DDF0", "#276FBF","#61BA5C", "#E2BF3C"))+
  scale_fill_manual(values=c("#A6DDF0", "#276FBF","#61BA5C", "#E2BF3C")) +
  scale_shape_manual(values=c(21,21,22,24))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) +
  ylim(0, 0.3)


fourpop_out <- ggarrange(obshet_plot_fourpop, exphet_plot_fourpop,fis_plot_fourpop, 
          common.legend=T, nrow=1)


ggsave(filename = "figures/heterozygosity_fourpop.png", fourpop_out,
        h=3.5, w=7)

#------------ six pop

obshet_plot_sixpop <- ggplot(sixpop, 
                              aes(x = population, 
                                  y = Obs_Het, 
                                  fill = population, shape=population)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = Obs_Het-StdErr.2, ymax = Obs_Het+StdErr.2), width = 0.2) +
  #scale_size_manual(values = c("All" = 8, "Population" = 4), guide = "none") +
  theme_classic(base_size = 12) +
  ylab("Observed\nheterozygosity") +
  scale_color_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7", "firebrick3", "grey60"))+
  scale_fill_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7", "firebrick3", "grey60"))+
  scale_shape_manual(values=c(21,22,23,24, 25, 21))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) +
  ylim(0.1, 0.3)

exphet_plot_sixpop <- ggplot(sixpop, 
                              aes(x = population, 
                                  y = Exp_Het, 
                                  fill = population, shape=population)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = Exp_Het-StdErr.3, ymax = Exp_Het+StdErr.3), width = 0.2) +
  scale_size_manual(values = c("All" = 8, "Population" = 4), guide = "none") +
  theme_classic(base_size = 12) +
  ylab("Expected\nheterozygosity") +
  scale_color_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7", "firebrick3", "grey60"))+
  scale_fill_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7", "firebrick3", "grey60"))+
  scale_shape_manual(values=c(21,22,23,24, 25, 21))+  
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) +
  ylim(0.1, 0.3)



#FIS

fis_plot_sixpop <- ggplot(sixpop, 
                           aes(x = population, 
                               y = Fis, 
                               fill = population, shape=population)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = Fis-StdErr.7, ymax = Fis+StdErr.7), width = 0.2) +
  theme_classic(base_size = 12) +
  ylab("Inbreeding\ncoefficient") +
  scale_color_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7", "firebrick3", "grey60"))+
  scale_fill_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7", "firebrick3", "grey60"))+
  scale_shape_manual(values=c(21,22,23,24, 25, 21))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) +
  ylim(0, 0.3)


sixpop_out <- ggarrange(obshet_plot_sixpop, exphet_plot_sixpop,fis_plot_sixpop, 
                         common.legend=T, nrow=1)


ggsave(filename = "figures/heterozygosity_sixpop.png", sixpop_out,
       h=3.5, w=7)



#----------------------------------------------------------------
#----------------------------------------------------------------
# no hybrids

obshet_plot_fourpop <- ggplot(fourpop_nohyb, 
                              aes(x = population, 
                                  y = Obs_Het, 
                                  fill = population, shape=population)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = Obs_Het-StdErr.2, ymax = Obs_Het+StdErr.2), width = 0.2) +
  theme_classic(base_size = 12) +
  ylab("Observed\nheterozygosity") +
  scale_color_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_fill_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,22,23,24))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4)))+
  ylim(0.1, 0.3)

exphet_plot_fourpop <- ggplot(fourpop_nohyb, 
                              aes(x = population, 
                                  y = Exp_Het, 
                                  fill = population, shape=population)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = Exp_Het-StdErr.3, ymax = Exp_Het+StdErr.3), width = 0.2) +
  scale_size_manual(values = c("All" = 8, "Population" = 4), guide = "none") +
  theme_classic(base_size = 12) +
  ylab("Expected\nheterozygosity") +
  scale_color_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_fill_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,22,23,24))+  
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) +
  ylim(0.1, 0.3)



#FIS

fis_plot_fourpop <- ggplot(fourpop_nohyb, 
                           aes(x = population, 
                               y = Fis, 
                               fill = population, shape=population)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = Fis-StdErr.7, ymax = Fis+StdErr.7), width = 0.2) +
  theme_classic(base_size = 12) +
  ylab("Inbreeding\ncoefficient") +
  scale_color_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_fill_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,22,23,24))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) +
  ylim(0, 0.3)


fourpop_out <- ggarrange(obshet_plot_fourpop, exphet_plot_fourpop,fis_plot_fourpop, 
                         common.legend=T, nrow=1)


ggsave(filename = "figures/heterozygosity_fourpop_nohyb.png", fourpop_out,
       h=3.5, w=7)

#------------ six pop

obshet_plot_sixpop <- ggplot(sixpop_nohyb, 
                             aes(x = population, 
                                 y = Obs_Het, 
                                 fill = population, shape=population)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = Obs_Het-StdErr.2, ymax = Obs_Het+StdErr.2), width = 0.2) +
  #scale_size_manual(values = c("All" = 8, "Population" = 4), guide = "none") +
  theme_classic(base_size = 12) +
  ylab("Observed\nheterozygosity") +
  scale_color_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7", "firebrick3", "grey60"))+
  scale_fill_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7", "firebrick3", "grey60"))+
  scale_shape_manual(values=c(21,22,23,24, 25, 21))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) +
  ylim(0.1, 0.3)

exphet_plot_sixpop <- ggplot(sixpop_nohyb, 
                             aes(x = population, 
                                 y = Exp_Het, 
                                 fill = population, shape=population)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = Exp_Het-StdErr.3, ymax = Exp_Het+StdErr.3), width = 0.2) +
  scale_size_manual(values = c("All" = 8, "Population" = 4), guide = "none") +
  theme_classic(base_size = 12) +
  ylab("Expected\nheterozygosity") +
  scale_color_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7", "firebrick3", "grey60"))+
  scale_fill_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7", "firebrick3", "grey60"))+
  scale_shape_manual(values=c(21,22,23,24, 25, 21))+  
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) +
  ylim(0.1, 0.3)



#FIS

fis_plot_sixpop <- ggplot(sixpop_nohyb, 
                          aes(x = population, 
                              y = Fis, 
                              fill = population, shape=population)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = Fis-StdErr.7, ymax = Fis+StdErr.7), width = 0.2) +
  theme_classic(base_size = 12) +
  ylab("Inbreeding\ncoefficient") +
  scale_color_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7", "firebrick3", "grey60"))+
  scale_fill_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7", "firebrick3", "grey60"))+
  scale_shape_manual(values=c(21,22,23,24, 25, 21))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) +
  ylim(0, 0.3)


sixpop_out <- ggarrange(obshet_plot_sixpop, exphet_plot_sixpop,fis_plot_sixpop, 
                        common.legend=T, nrow=1)


ggsave(filename = "figures/heterozygosity_sixpop_nohyb.png", sixpop_out,
       h=3.5, w=7)




#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# pixy
### four pop
dat_pops <- read.csv("analysis/diversity/1kb_fourpop_pi.txt", header=T, sep="\t")
nrow(dat_pops)

dat.na <- dat_pops[!is.na(dat_pops$avg_pi),]
nrow(dat.na)
#183556

dat_pops.sub <- dat.na[grep("NC_",dat.na$chromosome),]
dat_pops.sub <- dat_pops.sub[grep("NC_047055.1",dat_pops.sub$chromosome, invert=T),]
nrow(dat_pops.sub)
# 177036
table(dat_pops.sub$chromosome)
table(dat_pops.sub$pop)
# 44259 for all

ggplot(dat_pops.sub, aes(x=pop, y=(avg_pi))) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 0.002))
sum(dat_pops.sub$avg_pi == 0)

sum(dat_pops.sub$count_diffs)/sum(dat_pops.sub$count_comparisons)
#[1] 0.0008097222


pi <- dat_pops.sub %>%
  group_by(pop) %>%
  summarise(overall_pi = sum(count_diffs) / sum(count_comparisons))

  # Load required libraries
  library(dplyr)
library(boot)
library(purrr)  
  

# function to calculate pi
calculate_pi <- function(data, indices) {
  d <- data[indices,] # need to specify indices here bc this is how boot samples. 
  pi <- sum(d$count_diffs) / sum(d$count_comparisons)
  return(pi)
}

# Function to run bootstrap 
bootstrap_population <- function(pop_data) {
  set.seed(123) 
  boot_results <- boot(data = pop_data, statistic = calculate_pi, R = 10000)
  ci <- boot.ci(boot_results, type = "perc")
  
  return(list(
    population = unique(pop_data$pop),
    estimate = boot_results$t0,
    ci_lower = ci$percent[4],
    ci_upper = ci$normal[5]
    ci = ci
  ))
}

# run bootstrap
results <- dat_pops.sub %>%
  group_by(pop) %>%
  group_split() %>% # create of df for each pop
  map_dfr(bootstrap_population) # apply function to each df in list

results
#population       estimate ci_lower ci_upper
#<chr>               <dbl>    <dbl>    <dbl>
#  1 Coastal_Atlantic 0.000696 0.000687 0.000705
#2 Coastal_Gulf     0.000641 0.000631 0.000651
#3 Intermediate     0.000940 0.000928 0.000951
#4 Offshore         0.00142  0.00141  0.00144 


results$population <- factor(results$population)
#legend_order <- c("All Populations", "Atlantic", "Dry Tortuga", "NGOMex", "WGOMex")
#results_combined$population_legend <- factor(results_combined$population, levels = legend_order)

results$population <- c("CoastAtl", "CoastGulf", "Intermediate", "Offshore")

piplot <- ggplot(results, 
                 aes(x = population, 
                     y = estimate, 
                     fill = population, shape=population)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  theme_classic(base_size = 12) +
  ylab("Genetic\ndiversity") +
  #scale_color_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  #scale_fill_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,21,22,24))+    
  xlab("") +
  ylim(0, 0.0015)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) +
  scale_fill_manual(values=c("#A6DDF0", "#276FBF","#61BA5C", "#E2BF3C"))
  


piplot

ggsave("figures/pi_fourpop.png", h=5, w=6)







### merge plots

fourpop_out <- ggarrange(obshet_plot_fourpop, exphet_plot_fourpop,
                         fis_plot_fourpop, piplot,
                         common.legend=T, nrow=2, ncol=2, labels="AUTO")
fourpop_out



ggsave(file="figures/diversity.pdf",fourpop_out,
       w=7, h=5)
ggsave(file="figures/diversity.png",fourpop_out,
       w=7, h=5)
















#-------------------------
# six pop
#-------------------------------------------------------------------------------


dat_pops <- read.csv("analysis/diversity/1kb_sixpop_pi.txt", header=T, sep="\t")
nrow(dat_pops)

dat.na <- dat_pops[!is.na(dat_pops$avg_pi),]
nrow(dat.na)

dat_pops.sub <- dat.na[grep("NC_",dat.na$chromosome),]
dat_pops.sub <- dat_pops.sub[grep("NC_047055.1",dat_pops.sub$chromosome, invert=T),]
nrow(dat_pops.sub)
table(dat_pops.sub$chromosome)
table(dat_pops.sub$pop)
# 44259 for all

ggplot(dat_pops.sub, aes(x=pop, y=(avg_pi))) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 0.002))
sum(dat_pops.sub$avg_pi == 0)

sum(dat_pops.sub$count_diffs)/sum(dat_pops.sub$count_comparisons)
#[1] 0.0007588548


pi <- dat_pops.sub %>%
  group_by(pop) %>%
  summarise(overall_pi = sum(count_diffs) / sum(count_comparisons))
pi
  
# Load required libraries
library(dplyr)
library(boot)
library(purrr)  


# function to calculate pi
calculate_pi <- function(data, indices) {
  d <- data[indices,] # need to specify indices here bc this is how boot samples. 
  pi <- sum(d$count_diffs) / sum(d$count_comparisons)
  return(pi)
}

# Function to run bootstrap 
bootstrap_population <- function(pop_data) {
  set.seed(123) 
  boot_results <- boot(data = pop_data, statistic = calculate_pi, R = 10000)
  ci <- boot.ci(boot_results, type = "perc")
  
  return(list(
    population = unique(pop_data$pop),
    estimate = boot_results$t0,
    ci_lower = ci$percent[4],
    ci_upper = ci$percent[5]
  ))
}

# run bootstrap
results <- dat_pops.sub %>%
  group_by(pop) %>%
  group_split() %>% # create of df for each pop
  map_dfr(bootstrap_population) # apply function to each df in list


results

#population       estimate ci_lower ci_upper
#<chr>               <dbl>    <dbl>    <dbl>
#1 Coastal_Atlantic      0.000696 0.000687 0.000705
#2 Coastal_Gulf          0.000641 0.000631 0.000651
#3 Intermediate_Atlantic 0.000935 0.000923 0.000947
#4 Intermediate_Gulf     0.000940 0.000927 0.000953
#5 Offshore_Atlantic     0.00142  0.00140  0.00144 
#6 Offshore_Gulf         0.00140  0.00138  0.00141 


results$population <- factor(results$population)
#legend_order <- c("All Populations", "Atlantic", "Dry Tortuga", "NGOMex", "WGOMex")
#results_combined$population_legend <- factor(results_combined$population, levels = legend_order)

piplot <- ggplot(results, 
                 aes(x = population, 
                     y = estimate, 
                     fill = population, shape=population)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  theme_classic(base_size = 12) +
  ylab("Genetic\ndiversity") +
  #scale_color_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  #scale_fill_manual(values=c("grey75","#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,22,23,24, 21, 22))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4))) 


piplot

ggsave("figures/pi_sixpop.png", h=5, w=6)

