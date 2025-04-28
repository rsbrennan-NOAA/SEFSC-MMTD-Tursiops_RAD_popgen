
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
library(ggplot2)
library(ggpubr)

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
  scale_color_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_fill_manual(values=c("#E69F00","#56B4E9", "#009E73", "#CC79A7"))+
  scale_shape_manual(values=c(21,22,23,24))+    
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4)),  
         shape = guide_legend(override.aes = list(size = 4)))+
  ylim(0.1, 0.3)

exphet_plot_fourpop <- ggplot(fourpop, 
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

fis_plot_fourpop <- ggplot(fourpop, 
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















# Modify each plot to remove x-axis labels except bottom plots
piplot_mod <- piplot + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

obshet_plot_mod <- obshet_plot + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())


# Combine the plots
grid_plots <- (guide_area() / (piplot_mod + obshet_plot_mod) / (exphet_plot + fis_plot)) +
  plot_layout(guides = "collect", heights = (c(0.1, 1,1))) +
  plot_annotation(tag_levels = 'A') &
  theme(legend.position = 'top',plot.tag = element_text(face = 'bold', size=16))

grid_plots

ggsave(file="../figures/fig4.pdf",grid_plots,
       w=7, h=5)
ggsave(file="../figures/fig4.png",grid_plots,
       w=7, h=5)