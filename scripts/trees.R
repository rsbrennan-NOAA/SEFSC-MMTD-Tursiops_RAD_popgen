# mito tree

library(tidyverse)
library(ggtree)
library(treeio)
library(ggnewscale)
library(ggtreeExtra)
library(dplyr)


grouping <- read_csv("analysis/trees/RADseq_mtDNA_labels.csv")
tree <- read.tree(file="./analysis/trees/Run7.treefile")
haplotype_counts <- grouping %>%
  count(Haplotype, fourpop, name = "n")


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
      fill = fourpop
    ),
    offset = 0.01,
    pwidth = 0.4,
    axis.params = list(axis = "x", text.size = 3)
  )  +
  # set colors
    scale_fill_manual(
    name = "Population",
    values = c("#4782d4", "#e1526b", "#61BA5C", "#E2BF3C", "grey45", "black")
  ) +
  geom_nodepoint(data = plot_data, 
                 aes(subset = bootstrap > 70, color = bootstrap), 
                 size = 3.5, alpha = 1) +
  scale_color_gradient(
    name = "Bootstrap (%)",
    low = "black",
    high = "green"
  )+
  theme_tree() 
  #geom_text2(aes(label = node, subset = !isTip), hjust = -0.3, size = 3) 
  #geom_tiplab(size = 1)


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


pout <- p %<+% grouping + 
  geom_tiplab(aes(color=group), size=0.9) +
  theme(legend.position="right")+ 
  #geom_text( show.legend  = F ) +
  geom_tippoint(aes(color=group), size=0.9) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = c("A_hudsonica" = "#1B9E77",
                                "A_lilljeborgi" = "#D95F02",
                                "F" = "#7570B3",
                                "IV" = "#E7298A",
                                "out_group" = "#666666",
                                "S" = "#66A61E",
                                "SB" = "#E6AB02",
                                "X" = "#A6761D"),
                     na.value = "black")

ggsave(pout, file="./output/tree_plot.pdf", h=13, w=10)
