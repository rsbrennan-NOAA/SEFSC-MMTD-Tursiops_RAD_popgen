# mito tree

library(tidyverse)
library(ggtree)
library(treeio)
library(ggnewscale)
library(ggtreeExtra)
library(dplyr)


grouping <- read_csv("analysis/trees/RADseq_mtDNA_labels.csv")
tree <- read.tree(file="./analysis/trees/Run7.treefile")
#tree <- read.tree(file="./analysis/trees/Run8.treefile")
grouping$sixpop <- paste0(grouping$fourpop ,grouping$region...6)
haplotype_counts <- grouping %>%
  count(Haplotype, sixpop, name = "n")


# Create the base tree plot
p <- ggtree(tree)

plot_data <- p$data %>%
  filter(!is.na(label) & !isTip) %>%
  # split by /
  separate(label, into = c("alrt", "bootstrap"), sep = "/", remove = FALSE) %>%
  mutate(bootstrap = as.numeric(bootstrap))

# remove nodes of collapsed outgroups:
plot_data <- plot_data[!plot_data$node %in% seq(84,92,1),]


p2 <- p +
  geom_tiplab(
    aes(label = ""), # empty text
    align = TRUE,
    linetype = "dashed",
    linesize = 0.05,
    offset = 0.001,
    color="grey45"
  ) +
  
  # add barplot of haplotype freqs
  geom_fruit(
    data = haplotype_counts,
    geom = geom_col,
    mapping = aes(
      y = Haplotype,
      x = n,
      #fill = fourpop
      fill = sixpop
    ),
    offset = 0.01,
    pwidth = 0.4,
    axis.params = list(axis = "x", text.size = 3)
  )  +
  # set colors
    scale_fill_manual(
    name = "Population",
    values = c("#4782d4", "#e1526b", "#61BA5C", "#E2BF3C", "grey45", "black",
               "purple", "orange", "tan4")
  ) +
  #geom_nodepoint(data = plot_data, 
  #               aes(subset = bootstrap > 70, color = bootstrap), 
  #              #aes( color = bootstrap), 
  #               size = 3.5, alpha = 1) +
  geom_nodelab(data = plot_data, 
              # aes(subset = bootstrap > 70, label = bootstrap), 
               aes(label = bootstrap), 
               #aes( color = bootstrap), 
               size = 3, alpha = 1,
               nudge_x = -.0035, nudge_y = 1.1) +
  scale_color_gradient(
    name = "Bootstrap (%)",
    low = "red",
    high = "green"
  )+
  theme_tree() 
  #geom_text2(aes(label = node, subset = !isTip), hjust = -0.3, size = 3) 
  #geom_tiplab(size = 1)

p2
p3 <- collapse(p2,node = 84,mode = "max",fill="transparent",
           color="black",size=0.1)

p2
p3

ggsave(filename="figures/mtTree.pdf", p3, h=5, w=7)
ggsave(filename="figures/mtTree.png", p3, h=5, w=7)


p$data$label[p$data$isTip == TRUE][!p$data$label[p$data$isTip == TRUE] %in% haplotype_counts$Haplotype]
haplotype_counts$Haplotype[!haplotype_counts$Haplotype %in% p$data$label[p$data$isTip == TRUE]]

4Tt065



#---------------------------------------------------
# plot RADtree:



grouping <- read_csv("analysis/trees/RADseq_mtDNA_labels.csv")
tree <- read.tree(file="./analysis/trees/Run2.treefile")






grouping_for_plot <- grouping %>%
  mutate(label = indiv) %>%
  select(label, fourpop)

p <- ggtree(tree) %<+% grouping_for_plot

plot_data <- p$data %>%
  # Only work on internal nodes that have a label
  filter(!is.na(label) & !isTip) %>%
  # Split the 'label' column into two new columns at the "/"
  separate(label, into = c("alrt", "bootstrap"), sep = "/", remove = FALSE) %>%
  # Convert the new bootstrap column to a number
  mutate(bootstrap = as.numeric(bootstrap))


# Get unique population groups, excluding NA values
pop_groups <- grouping %>%
  filter(!is.na(fourpop) & fourpop != "NA") %>%
  pull(fourpop) %>%
  unique()

plot_data

p + 
  geom_text2(aes(label = node, subset = !isTip), hjust = -0.3, size = 3) +
  geom_tiplab(size = 1)

#collapse(p,node = 474,mode = "max",fill="transparent",
#         color="black",size=0.1)



pout <-p +
  # First color scale for the discrete 'fourpop' variable
  geom_tiplab(aes(color = fourpop), size = 3, align = TRUE, linetype = "dotted") +
  scale_color_manual(
    name = "Population", 
    values = c("#4782d4", "#e1526b", "#61BA5C", "#E2BF3C", "grey45", "black")
  ) +
  # Start a new color scale
  new_scale_color() +theme_tree()+
  
  # Second color scale for the continuous 'bootstrap' variable
  #geom_nodepoint(data = plot_data,
  #               aes(subset = bootstrap > 70, color = bootstrap),
  #               size = 2.5, alpha = 1) +
  scale_color_gradient(
    name = "Bootstrap (%)",
    low = "black",   # A light green for low values
    high = "green"  # A dark green for high values
  )+
  theme_tree()




ggsave(pout, filename="figures/rad_tree.pdf", h=6, w=7)
ggsave(pout, filename="figures/rad_tree.png", h=6, w=7)

#


## --------------------------------------------------------------------
## --------------------------------------------------------------------
## --------------------------------------------------------------------
## --------------------------------------------------------------------
# mito tree with all indivs:
grouping <- read_csv("analysis/trees/RADseq_mtDNA_labels.csv")
tree <- read.tree(file="./analysis/trees/Run9.treefile")

grouping_for_plot <- grouping %>%
  mutate(label = indiv) %>%
  select(label, fourpop)

p <- ggtree(tree) %<+% grouping_for_plot

plot_data <- p$data %>%
  # Only work on internal nodes that have a label
  filter(!is.na(label) & !isTip) %>%
  # Split the 'label' column into two new columns at the "/"
  separate(label, into = c("alrt", "bootstrap"), sep = "/", remove = FALSE) %>%
  # Convert the new bootstrap column to a number
  mutate(bootstrap = as.numeric(bootstrap))


# Get unique population groups, excluding NA values
pop_groups <- grouping %>%
  filter(!is.na(fourpop) & fourpop != "NA") %>%
  pull(fourpop) %>%
  unique()

plot_data

pout <-p +
  # First color scale for the discrete 'fourpop' variable
  #geom_tiplab(aes(color = fourpop), size = 3, align = TRUE, linetype = "dotted") +
  geom_tippoint(aes(color = fourpop), size=0.8) +
  
  scale_color_manual(
    name = "Population", 
    values = c("#4782d4", "#e1526b", "#61BA5C", "#E2BF3C", "grey45", "black")
  ) +
  # Start a new color scale
  new_scale_color() +theme_tree()+
  
  # Second color scale for the continuous 'bootstrap' variable
  #geom_nodepoint(data = plot_data,
  #               aes(subset = bootstrap > 70, color = bootstrap),
  #               size = 2.5, alpha = 1) +
  scale_color_gradient(
    name = "Bootstrap (%)",
    low = "black",   # A light green for low values
    high = "green"  # A dark green for high values
  )+
  theme_tree()

pout


ggsave(pout, filename="figures/mt_tree_allindiv_points.pdf", h=6, w=7)
ggsave(pout, filename="figures/mt_tree_allindiv_points.png", h=6, w=7)

































#ggsave(pout, filename="figures/rad_tree.pdf", h=6, w=7)
#ggsave(pout, filename="figures/rad_tree.png", h=6, w=7)







p +
  geom_tiplab(aes(color = fourpop), size = 3, align = TRUE) +
  geom_tippoint(aes(color = fourpop), size = 2) +
  theme_tree() +
  labs(color = "Population Group") + # Legend title
  scale_color_manual(values=c("#4782d4", "#e1526b","#61BA5C", "#E2BF3C", "grey45", "black"))+
  geom_nodelab(hjust = -0.2, size = 2.5, color = "firebrick")

p <- ggtree(tree) %<+% grouping

x <- as_tibble(tree)

head(x, n=15)

x$label

dm <- left_join(x, grouping, by="label")





p +
  # First color scale for the discrete 'fourpop' variable
  geom_tiplab(aes(color = fourpop), size = 3, align = TRUE, linetype = "dotted") +
  scale_color_manual(
    name = "Population", 
    values = c("#4782d4", "#e1526b", "#61BA5C", "#E2BF3C", "grey45", "black")
  ) +
  # Start a new color scale
  new_scale_color() +theme_tree()+
  
  # Second color scale for the continuous 'bootstrap' variable
  geom_nodepoint(data = plot_data,
                 aes(subset = bootstrap > 70, color = bootstrap),
                 size = 2.5, alpha = 1) +
  scale_color_gradient(
    name = "Bootstrap (%)",
    low = "#edf8e9",   # A light green for low values
    high = "forestgreen"  # A dark green for high values
  )+
  theme_tree()




















p + 
  geom_text2(aes(label = node, subset = !isTip), hjust = -0.3, size = 3) +
  geom_tiplab(size = 1)


    
    




table(grouping$Haplotype)

# Basic piechart
ggplot(data, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)








p +
  # Add tip labels colored by population
  geom_tiplab(aes(color = fourpop), size = 3, align = TRUE, linetype = "dotted") +
  
  # NEW: Add colored circles to nodes based on bootstrap support
  geom_nodepoint(data = plot_data, 
                 aes(subset = bootstrap > 70, color = bootstrap), 
                 size = 3.5, alpha = 0.8) +
  
  # Use a custom, intuitive color gradient for bootstrap values
  scale_color_gradientn(
    name = "Bootstrap (%)", 
    colors = c("blue", "firebrick", "red"), 
    limits = c(70, 100)
  ) +
  
  # A clean theme for the tree
  theme_tree() +
  
  # Adjust legend for the population colors
  guides(color = guide_legend(override.aes = list(label = ""))) +
  
  labs(fill = "Population Group") # Corrected labs call




p +
  geom_tiplab(aes(color = fourpop), size = 3, align = TRUE) +
  #geom_tippoint(aes(color = fourpop), size = 4) +
  # Add node labels using our new, clean data
  #geom_nodelab(data = plot_data, aes(label = bootstrap, subset = bootstrap > 70, fill=bootstrap), 
  #             size = 2.5, color = "black") +
  geom_nodepoint(data = plot_data, 
                 aes(subset = bootstrap > 70, color = bootstrap), 
                 size = 3.5, alpha = 0.8) +
  theme_tree() +
  labs(color = "Population Group")+
  geom_label(data=plot_data,aes(label = bootstrap, fill = bootstrap,subset = bootstrap > 70))





p +
  geom_tiplab(aes(color = fourpop), size = 3, align = TRUE) +
  geom_tippoint(aes(color = fourpop), size = 2) +
  theme_tree() +
  labs(color = "Population Group") + # Legend title
  scale_color_manual(values=c("#4782d4", "#e1526b","#61BA5C", "#E2BF3C", "grey45", "black"))+
  geom_nodelab(hjust = -0.2, size = 2.5, color = "firebrick")

p <- ggtree(tree) %<+% grouping

x <- as_tibble(tree)

head(x, n=15)

x$label

dm <- left_join(x, grouping, by="label")



p <- ggtree(tree) + theme_tree()

ggtree(tree) +
geom_tiplab() +
  theme_tree2() 






#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# cophylo tree
# https://github.com/rapidspeciation/mechanitis_melinaea/blob/main/scripts/cophylo%20plot.R
require("ape")
require(phytools)

#### Read in the data ####
# #note: the tip-labels have been renamed to include species and location of the individuals to aid the plotting
mech<-read.tree("~/Phylo.Mech.base.thin500.minGQ10.PhyloInd.ALL.wB.treefile.renamed")

# mt tree:
grouping <- read_csv("analysis/trees/RADseq_mtDNA_labels.csv")
mech_mt <- read.tree(file="./analysis/trees/Run7.treefile")

mech_mt <- read.tree(file="./analysis/trees/Run9.treefile")

grouping_for_plot <- grouping %>%
  mutate(label = indiv) %>%
  select(label, fourpop)

# rad tree:
grouping <- read_csv("analysis/trees/RADseq_mtDNA_labels.csv")
mech <- read.tree(file="./analysis/trees/Run2.treefile")

grouping_for_plot <- grouping %>%
  mutate(label = indiv) %>%
  select(label, fourpop)



# Root Mechanitis on Forbestra
mech <- midpoint.root(mech)
mech_mt <- midpoint.root(mech_mt)

# Ladderize the trees
mech<-ladderize(mech)
mech_mt<-ladderize(mech_mt)

mech$tip.label
mech_mt$tip.label

# Get a dataframe with the association (here just twice the same individual)
assoc<-cbind(mech$tip.label, mech$tip.label)
scale_color_manual(values=c("#4782d4", "#e1526b","#61BA5C", "#E2BF3C", "grey45", "black"))+
  
speciesCol<-as.data.frame(assoc[,1])
colnames(speciesCol) <- "label"
grouping_for_plot$col <- ifelse(grepl(grouping_for_plot$fourpop,pattern="Coastal_Atlantic"),"#4782d4",
                                ifelse(grepl(grouping_for_plot$fourpop,pattern="Coastal_Gulf"),"#e1526b",
                                       ifelse(grepl(grouping_for_plot$fourpop,pattern="Intermediate"),"#61BA5C",
                                              ifelse(grepl(grouping_for_plot$fourpop,pattern="Offshore"),"#E2BF3C",
                                                     ifelse(grepl(grouping_for_plot$fourpop,pattern="PutHyb_Intermediate"),"grey45",
                                                            ifelse(grepl(grouping_for_plot$fourpop,pattern="PutHyb_Offshore"),"black","white"
                                                              ))))))
colsin <- speciesCol %>%
  left_join(grouping_for_plot, by = "label")

obj<-cophylo(mech,mech_mt,assoc=assoc,cex=0.1)
plot.cophylo(obj,link.type="curved",link.lwd=1,fsize=0.3,pts=F,ftype="off",
             link.lty="solid",link.col=colsin$col)


