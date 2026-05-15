# Figure 2
library(ggplot2)
library(gridExtra)
library(grid)

# scale_color_manual(values=c("#A6DDF0", "#276FBF","#88D03F","#1B9E77", "#FFDD33", "#C49E45"))+
# scale_fill_manual(values=c("#A6DDF0", "#276FBF","#61BA5C", "#E2BF3C"))+


# Treemix:


######################################################## --------------------------
################## (B) Plot Tree #######################
########################################################

## 1. From the final runs, compare tree likelihoods, 
# select tree with highest likelihood, 
# remove duplicates and retain tree(s) with unique topology. 
#Adapted from R functions written by Zecca, Labra and Grassi (2019).
library(BITEV2)

setwd("~/Documents/projects/Tursiops_RAD_popgen/analysis/pop_structure/treemix/fourpop/bootstrap")
#folder with all TreeMix outputs from the final runs
source("~/Documents/projects/Tursiops_RAD_popgen/scripts/TreeMix_functions.R") #path to required functions for this analysis

maxLL("fourpop_treemix_bootrep_", nt=100, uel=FALSE)
#cfTrees("fourpop_treemix_bootrep_", nt=100, p=1, m='PH85')                        

## 2. Now plot and save unique tree with highest likelihood:

#pdf("../../../../../figures/TreeMix_output_fourpop.pdf", h=4, w=5)                                          
treemix.bootstrap(in.file="fourpop_treemix_bootrep_9", #stem of TreeMix files (as above) + number of run with highest likelihood and unique topology
                  out.file = "tmp", 
                  phylip.file = "../fourpop_outtree.newick", 
                  #consensus tree in newick format (from the bootstrap procedure generated with PHYLIP)    
                  nboot = 100, #nboot is the number of bootstraps used
                  fill = TRUE,  
                  #pop.color.file = "col.txt",#specify colours with a tab delimited pop.color.file - first column is pop name, second the colour
                  boot.legend.location = "topleft")


# Figure 3
library(admixtools)
library(tidyverse)
library(igraph)
library(plotly)

graph <- read_delim("~/Documents/projects/Tursiops_RAD_popgen/analysis/pop_structure/admixtools/graph_out.txt")

# manually plot so I have control:
library(magrittr)

edges = graph %>% as_tibble
names(edges)[1:2] = c('V1', 'V2')
graph = igraph::graph_from_edgelist(as.matrix(graph)[,1:2])  
leaves = get_leafnames(graph)
#edges = as_tibble(edges, .name_repair = ~c('V1', 'V2'))
admixnodes = unique(edges[[2]][edges[[2]] %in% names(which(table(edges[[2]]) > 1))])

pos = data.frame(names(V(graph)), igraph::layout_as_tree(graph), stringsAsFactors = F) %>%
  set_colnames(c('node', 'x', 'y'))
eg = graph %>% igraph::as_edgelist() %>%
  as_tibble(.name_repair = ~c('V1', 'V2')) %>%
  left_join(pos, by=c('V1'='node')) %>%
  left_join(pos %>% transmute(V2=node, xend=x, yend=y), by='V2') %>%
  mutate(type = ifelse(V2 %in% admixnodes, 'admix', 'normal')) %>% rename(name=V1, to=V2)

shift_subgraph_down = function(graph, coord, node, by) {
  # returns new coords with some nodes shifted down
  if(by <= 0) return(coord)
  nodes = graph %>% subcomponent(node, mode = 'out') %>% names
  coord %>% mutate(shift = name %in% nodes,
                   shift1 = to %in% nodes,
                   y = ifelse(shift, y-by, y),
                   yend = ifelse(shift | shift1, yend-by, yend)) %>%
    select(-shift, -shift1)
}
fix_shiftdown = function(coord, graph) {
  adm = find_admixedges(graph) %>%
    left_join(coord, by = c('from'='name', 'to'='to')) %>%
    filter(yend >= y) %>% mutate(by = yend - y + 1)
  for(i in seq_len(nrow(adm))) {
    by = coord %>% filter(name == adm$from[i], to == adm$to[i]) %$% {yend - y + 1}
    coord = shift_subgraph_down(graph, coord, adm$to[i], by = by)
  }
  coord
}
eg = fix_shiftdown(eg, graph)

lab = ifelse(edges$type == 'admix', paste0(round(edges$weight*100), '%'), round(edges$weight*1000))
edges %<>% mutate(label = lab)
eg %<>% left_join(edges %>% transmute(name=V1, to=V2, label), by=c('name', 'to')) 



internal = NULL
nodes = eg %>% filter(to %in% leaves) %>%
  rename(x=xend, y=yend, xend=x, yend=y) %>% transmute(name = to, x, y, xend, yend)
pdat <- namedList(eg, nodes, internal)


pdat[["comb"]] = bind_rows(select(pdat$nodes, name, x, y), 
                           select(pdat$eg, name, x, y)) %>% distinct()
layers = vector("list")

pdat$eg %<>% mutate(label = ifelse(type == "normal", 
                                   "", as.character(label)))

# manuall adjust the offshore to make the intermediate physically intermediate between coastal and offshore:

pdat$eg$xend[1] <- -0.5
#pdat$eg$y[1] <- 2
pdat$eg$yend[1] <- 1
pdat$eg$xend[4] <- -0.5
pdat$eg$yend[4] <- 1
pdat$eg$xend[5] <- -1.5
pdat$eg$x[3] <- -0.5
pdat$eg$xend[3] <- -0.5

pdat$eg$y[3] <- 1
pdat$eg$yend[3] <- 0

pdat$nodes$x[3] <- -1.5
pdat$nodes$xend[3] <- -1.5
pdat$nodes$x[2] <- -0.5
pdat$nodes$xend[2] <- -0.5
pdat$nodes$y[2] <- 0
pdat$nodes$yend[2] <- 0

pdat$nodes$name <- gsub("_", " ", pdat$nodes$name)

admix_plot <- ggplot(pdat$eg, aes(x = x, xend = xend, y = y, yend = yend)) + 
  ylab("")+
  theme(panel.background = element_blank(), 
        axis.line = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), legend.position = "none") +
  #scale_linetype_manual(values = c("admix" = "dashed", "normal" = "solid")) +
  
  scale_fill_manual(values=c("#4782d4", "#e1526b","#61BA5C", "#E2BF3C")) +
  
  geom_label(data = pdat$nodes, 
             aes_string(label = "name", 
                        fill = "as.factor(name)"), 
             size = 2.5,
             nudge_y=-0.1,label.size = NA) +
  geom_text(aes(x = (x + xend)/2, 
                y = (y + yend)/2, 
                label = label), size = 3, nudge_y = -0.15)+ 
  xlab("")+
  scale_linetype_manual(values = c(admix = 3, 
                                   normal = 1))+
  scale_x_continuous(expand = c(0.15,  0.15)) +
  geom_segment(aes(linetype = as.factor(type)),
               
               linewidth = 0.5) +
  geom_segment(aes(x = xend - 0.01*(xend-x), 
                   xend = xend, 
                   y = yend - 0.01*(yend-y), 
                   yend = yend),
               linetype = 1, # Force solid line for arrowheads
               arrow = arrow(type = "closed", 
                             length = unit(0.12, "inches"),
                             angle = 20),
               linewidth = 0.5)

admix_plot 

ggsave("figures/admix_graph.pdf", h=3, w=4)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# triangle plot

setwd("~/Documents/projects/Tursiops_RAD_popgen/")

##-------------------------------------------------------------------------------
##-------------------------------------------------------------------------------
##-------------------------------------------------------------------------------
# triangle plot
library(triangulaR)
library(vcfR)

vcf <- read.vcfR( "analysis/variants/filtered.final_ids.vcf.gz", verbose = FALSE )
vcf

pops1 <- read.table("analysis/population_assignments_hybrids_summary.txt", header=T)
head(pops1)

table(pops1$newhybrids_category,pops1$sixpop)

pops <- data.frame(id = pops1$indiv,
                   pop = pops1$fourpop)

# identify ancestry informative markers (AIMS)
# Create a new vcfR object composed only of sites above the given allele frequency difference threshold

# [1] "438 sites passed allele frequency difference threshold"


diff2 <- alleleFreqDiff(vcfR = vcf, 
                        pm = pops, 
                        p1 = "Coastal_Atlantic", 
                        p2 = "Offshore", 
                        difference = 0.7)
#[1] "450 sites passed allele frequency difference threshold"

# Calculate hybrid index and heterozygosity for each sample. 
# Values are returned in a data.frame


hi.het1 <- hybridIndex(vcfR = diff2, 
                       pm = pops, 
                       p1 = "Coastal_Atlantic", p2 = "Offshore")



pops <- read.table("analysis/population_assignments_hybrids_summary.txt", header=T)
#hybs <- pops$indiv[pops$offshore_putative_hybrids == TRUE]


#pops2 <- read.table("analysis/pop_structure/newhybrids/newhybrids_posteriors.txt", header=T)

#hybs <- pops$indiv[pops$offshore_putative_hybrids == TRUE]
#newhybrids <- pops2[pops2$sig == "High_Confidence" & pops2$highest_category != "Pure1" & pops2$highest_category != "Pure2" ,]

F2 <- pops$indiv[pops$newhybrids_category == "F2"]
Backcross2  <- pops$indiv[pops$newhybrids_category == "Backcross2"]
Backcross1 <- pops$indiv[pops$newhybrids_category == "Backcross1"]

hi.het1$pop2 <- hi.het1$pop
hi.het1$pop2[hi.het1$id %in% F2] <- "F2"
hi.het1$pop2[hi.het1$id %in% Backcross2] <- "Backcross2"
hi.het1$pop2[hi.het1$id %in% Backcross1] <- "Backcross1"

df <- merge(hi.het1, pops,  by.x="id", by.y="indiv")
df$pop <- as.factor(df$sixpop)
hybsOnly <- df[which(df$pop2 == "F2" | df$pop2 == "Backcross1" | df$pop2 == "Backcross2" ),]
#hybsOnly$pop2 <- as.factor(hybsOnly$pop2)



fill_colors <- c(
  "Coastal_Atlantic" = "#4782d4",
  "Coastal_Gulf" = "#e1526b",
  "Intermediate_Atlantic" = "#B4ED50",
  "Intermediate_Gulf" = "#2E8B57", 
  "Offshore_Atlantic" = "#FFDD33", 
  "Offshore_Gulf" = "#C49E45"
)


t1 <- ggplot(df, aes(x = hybrid.index, y = heterozygosity)) + 
  annotate("segment", x = 0.5, xend = 1, y = 1, yend = 0, color = "black") +
  annotate("segment", x = 0, xend = 0.5, y = 0, yend = 1, color = "black") +
  annotate("segment", x = 0, xend = 1, y = 0, yend = 0, color = "black") +
  stat_function(fun = function(hi) 2 * hi * (1 - hi), xlim = c(0, 1), 
                color = "black", linetype = "dashed") +
  geom_point(aes(fill = pop, shape = pop), size = 2, alpha=1) + 
  # Hybrid  with different color 
  geom_point(data = hybsOnly, 
             aes(x = hybrid.index, y = heterozygosity, 
                 fill = pop, 
                 shape = pop,
                 color= pop2,
                 stroke=1.25), 
             size = 3) + 
  # Scales
  scale_fill_manual(values = fill_colors) +
  scale_shape_manual(values = c(21, 21, 22, 22, 24, 24, 21, 21)) +
  scale_color_manual(values = c("purple", "darkorange1", "firebrick2")) +
  # Guides/legends
  guides(
    shape = guide_legend(title = "", 
                         override.aes = list(size = 2, stroke = 0.5), 
                         order = 1),
    fill = guide_legend(title = "", 
                        override.aes = list(size = 2, stroke = 0.5), 
                        order = 1),
    color = guide_legend(title = "Hybrids", override.aes = list(size = 2),
                         order = 2)
  ) + 
  # Labels and limits
  xlab("Hybrid Index") + 
  ylab("Interclass Heterozygosity") + 
  ylim(c(-0.0, 1.0)) + 
  xlim(c(-0.0, 1.0)) + 
  theme_classic(base_size = 12)

t1
t2 <- t1 + theme(legend.position = "none")


ggsave(filename="figures/triangle_hybrids_legend.pdf",t1, w=8.0, h=3.5)
ggsave(filename="figures/triangle_hybrids.pdf",t2, w=4, h=3.5)



#plt1<- ggarrange(p,p, t2, common.legend=T, nrow=1)
#ggsave(filename="figures/bgc_combined.pdf",
#       ggarrange(plt1,p2, widths=c(0.33, 0.33,1), common.legend=T, nrow=2),
#       h=6, w=9)

library(ggpubr)

plt1 <- ggarrange(admix_plot,t2,t2,t2, common.legend=TRUE, nrow=2, ncol=2, labels="AUTO")
plt1
#final_plot <- ggarrange(plt1, p2, heights=c(1.7, 1), common.legend=TRUE, nrow=2, labels=c("", "D"))
#final_plot

# Save the final combined plot
ggsave(filename="figures/manuscript/fig2_v5.pdf",
       plt1,
       height=5, width=7)





