#ibd
# pca genetic distance matrix, based on tutorial pipeline https://github.com/TheWangLab/algatr/blob/558519684f5d346ca295f68b184c7ed99cf9c97e/R/Gen_dist.R#L17
library(vcfR)
library(algatr)
library(ade4)
library(dplyr)

vcf <- read.vcfR( "analysis/filtered.final_ids_LDthin.vcf.gz", verbose = FALSE )
#genin <- vcfR2genind(vcf)


# compare imputation methods
simple_dos  <- simple_impute(gen_to_geno(vcf))
str_dos <- str_impute(gen = vcf, K = 1:7, 
                      entropy = TRUE, repetitions = 5, 
                      quiet = FALSE, save_output = FALSE)

str_dos <- str_impute(gen = vcf, K = 4)

str_dist <- as.matrix(ecodist::distance(str_dos, method = "euclidean"))
#simple_dist <- as.matrix(ecodist::distance(simple_dos, method = "euclidean"))


pc_dists <- gen_dist((str_dos), dist_type = "pc", 
                     npc_selection = "manual") # use 3pcs, following k-1
pc_dists_auto <- gen_dist((str_dos), dist_type = "pc", 
                     npc_selection = "auto",criticalpoint = 2.0234)

pc_plink <- gen_dist(
  dist_type = "plink",
  plink_file = "analysis/pop_structure/LDthin.dist",
  plink_id_file = "analysis/pop_structure/LDthin.dist.id",
  npc_selection = "manual"
)

gen_dist_corr(dist_x = pc_dists, dist_y = pc_plink, 
              metric_name_x = "SNMF_PCA", metric_name_y = "plink")
gen_dist_corr(dist_x = pc_dists, dist_y = str_dist, 
              metric_name_x = "SNMF_PCA", metric_name_y = "euclidian")
gen_dist_corr(dist_x = str_dist, dist_y = pc_plink, 
              metric_name_x = "euclidian", metric_name_y = "plink")
gen_dist_corr(dist_x = pc_plink, dist_y = pc_dists_auto, 
              metric_name_x = "pc_plink", metric_name_y = "pc_dists_auto")
gen_dist_corr(dist_x = pc_dists, dist_y = pc_dists_auto, 
              metric_name_x = "SNMF_PCA", metric_name_y = "pc_dists_auto")

gen_dist_corr(dist_x = pc_plink, dist_y = str_dist, 
              metric_name_x = "pc_plink", metric_name_y = "euclidian")

# Gather data
plot <- gen_dist_corr(dist_x = pc_plink, dist_y = str_dist, 
                      metric_name_x = "pc_plink", metric_name_y = "euclidian")

# Now, build plot
ggplot(plot$data, aes(x = pc_plink, y = str_dist) ) +
  geom_hex() +
  theme_bw() +
  scale_fill_continuous(type = "viridis")



# TREE --------------------
library(ggtree)
library(ggplot2)
library(ape)
library(dplyr)
library(tidyr)
out_tree <- nj(as.matrix(pc_plink))
out_tree <- nj(as.matrix(pc_dists))
out <- as_tibble(nj(as.matrix(pc_dists)))

pops <- read.table("analysis/pop_structure/sixpop_all.clust", header=F)
colnames(pops) <-c("ids", "population")
#out$label <- gsub("X", "", out$label)
dm <- left_join(out, pops, by=c("label" ="ids"))
dm <- dm %>%
  separate(population, into = c("Population", "region"), sep = "_", remove = FALSE) %>%
  mutate(region = ifelse(region == "Atl", "Atlantic", region))


my_colors <- c("#3f88c5", "#ef476f", "#ffd166")
p <- ggtree(out_tree, layout="ape") + theme_tree()

pout <- 
  p %<+% dm + 
  #geom_tiplab(aes(color=newpop), size=0.9) +
  theme(legend.position="right")+ 
  #geom_text( show.legend  = F ) +
  geom_tippoint(aes(fill=Population, shape=region), size=4, alpha=1) +
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=my_colors)+
  guides(fill=guide_legend(override.aes=list(shape=21))) 
pout

ggsave("figures/nj_tree.pdf", pout, h=5, w=7)
ggsave("figures/nj_tree.png", pout, h=5, w=7)
# shirk paper says use 64 pc model when we have sub structure
# only a need a few for initial structure, but these others increase ability
 # to detect futher structure. 
# I think we have substantial sub structure, so we want to use 64

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# pull in pairwise distances

lc_dist <- read.csv("analysis/lc_distances.csv", header=T)

as.matrix(pc_dists)
head(dm)

# need to convert the pairwise genetic distances to long format:
pc_dists_long <- as.data.frame(pc_dists)
pc_dists_long$indiv_1 <- rownames(pc_dists)

library(tidyr)
pc_long <- pc_dists_long %>%
  pivot_longer(cols = -indiv_1, 
               names_to = "indiv_2", 
               values_to = "pc_distance")

lc_long <- lc_dist %>%
  select(-X) 

# Merge the datasets
merged_distances <- merge(pc_long, lc_long, 
                          by = c("indiv_1", "indiv_2"), 
                          all = FALSE)

# add the region and pop:
pops <- read.table("analysis/pop_structure/sixpop_all.clust", header=F)
colnames(pops) <-c("ids", "population")
#out$label <- gsub("X", "", out$label)
pops <- pops %>%
  mutate(population = ifelse(population == "Coastal_Atl", "Coastal_Atlantic", population)) %>%
  separate(population, into = c("Population", "region"), sep = "_", remove = FALSE) %>%
  mutate(region = ifelse(region == "Atl", "Atlantic", region))



id_to_pop <- setNames(pops$population, pops$ids)
merged_distances$population_1 <- id_to_pop[merged_distances$indiv_1]
merged_distances$population_2 <- id_to_pop[merged_distances$indiv_2]

id_to_region <- setNames(pops$region, pops$ids)
merged_distances$region_1 <- id_to_region[merged_distances$indiv_1]
merged_distances$region_2 <- id_to_region[merged_distances$indiv_2]


head(merged_distances)

merged_distances <- merged_distances %>%
  mutate(comparison = case_when(
    population_1 == "Coastal_Atlantic" & population_2 == "Coastal_Atlantic" ~ "Coastal_Atlantic",
    population_1 == "Coastal_Gulf" & population_2 == "Coastal_Gulf" ~ "Coastal_Gulf",
    population_1 == "Intermediate_Atlantic" & population_2 == "Intermediate_Atlantic" ~ "Intermediate_Atlantic",
    population_1 == "Intermediate_Gulf" & population_2 == "Intermediate_Gulf" ~ "Intermediate_Gulf",
    population_1 == "Offshore_Atlantic" & population_2 == "Offshore_Atlantic" ~ "Offshore_Atlantic",
    population_1 == "Offshore_Gulf" & population_2 == "Offshore_Gulf" ~ "Offshore_Gulf",
    TRUE ~ "inter-population"
  ))


p1 <- ggplot(merged_distances, aes(x = leastcost_distance, y = pc_distance)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw(base_size = 14) +
  labs(y = "Genetic Distance", 
       x = "Geographic Distance",
       title = "Isolation by distance: all populations")
p1

ggsave("figures/ibd_allpops.png", p1, h=5, w=7)


p2 <- ggplot(merged_distances, aes(x = leastcost_distance, y = pc_distance)) +
  geom_point(aes(color = comparison),alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE, linewidth=2) +
  theme_bw(base_size = 14) +
  labs(y = "Genetic Distance", 
       x = "Geographic Distance",
       title = "Isolation by distance: within population") +
  facet_wrap(~comparison, ncol=3) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size=5)))

p2
ggsave("figures/ibd_withinpops.png", p2, h=7, w=9)

merged_distances[which(merged_distances$comparison == "Offshore_Atlantic"),]

# 7Tt312 and 13Tt073 are divergent and causing that weird pattern in offshore atlantic


#-----------------------------------------
# split by Gulf and Atlantic:
merged_distances <- merged_distances %>%
  mutate(comparison_ocean = case_when(
    str_detect(population_1, "Gulf") & str_detect(population_2, "Gulf") ~ "Gulf",
    str_detect(population_1, "Atlantic") & str_detect(population_2, "Atlantic") ~ "Atlantic",
    TRUE ~ "inter-ocean"
  ))

p2 <- ggplot(merged_distances, aes(x = leastcost_distance, y = pc_distance)) +
  geom_point(aes(color = comparison),alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE, linewidth=2) +
  theme_bw(base_size = 14) +
  labs(y = "Genetic Distance", 
       x = "Geographic Distance",
       title = "Isolation by distance: within region") +
  facet_wrap(~comparison_ocean, ncol=2) +
guides(color = guide_legend(override.aes = list(alpha = 1, size=5)))

p2
ggsave("figures/ibd_withinregion.png", p2, h=5, w=9)


#-------------------------
# split by coastal, intermediate, offshore

merged_distances <- merged_distances %>%
  mutate(comparison_depth = case_when(
    str_detect(population_1, "Coastal") & str_detect(population_2, "Coastal") ~ "Coastal",
    str_detect(population_1, "Intermediate") & str_detect(population_2, "Intermediate") ~ "Intermediate",
    str_detect(population_1, "Offshore") & str_detect(population_2, "Offshore") ~ "Offshore",
    TRUE ~ "inter-depth"
  ))

p3 <- ggplot(merged_distances, aes(x = leastcost_distance, y = pc_distance)) +
  geom_point(aes(color = comparison),alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE, linewidth=2) +
  theme_bw(base_size = 14) +
  labs(y = "Genetic Distance", 
       x = "Geographic Distance",
       title = "Isolation by distance: within depth") +
  facet_wrap(~comparison_depth, ncol=2)+
  guides(color = guide_legend(override.aes = list(alpha = 1, size=5)))

p3
ggsave("figures/ibd_withindepth.png", p3, h=5, w=7)


#-------------------------

# only keep inter pop comparisons:

merged_distances <- merged_distances %>%
  mutate(comparison_depth = case_when(
    str_detect(population_1, "Coastal") & str_detect(population_2, "Coastal") ~ "Coastal",
    str_detect(population_1, "Intermediate") & str_detect(population_2, "Intermediate") ~ "Intermediate",
    str_detect(population_1, "Offshore") & str_detect(population_2, "Offshore") ~ "Offshore",
    TRUE ~ "inter-depth"
  ))

merged_distances_dropped <- merged_distances[which(merged_distances$comparison == "inter-population"),]

p3 <- ggplot(merged_distances_dropped, aes(x = leastcost_distance, y = pc_distance)) +
  geom_point(aes(color = comparison),alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE, linewidth=2) +
  theme_bw(base_size = 14) +
  labs(y = "Genetic Distance", 
       x = "Geographic Distance",
       title = "Isolation by distance: Only inter-population") +
  facet_wrap(~comparison_depth, ncol=2)+
  guides(color = guide_legend(override.aes = list(alpha = 1, size=5)))

p3
ggsave("figures/ibd_withindepth_interOnly.png", p3, h=5, w=7)
















#------------------------------------------------------------------------
# plot just the simple method, nj tree. 

## Load packages
library('vcfR')
library('adegenet')

#read the data
vcf <- read.vcfR("Tt_mappedTt_SP1060_relax_all_7000_minQ25_GQ20_miss_pop5indMAF_paral_nomissSP1060.recode.vcf")

vcf

x <- vcfR2genlight(vcf)


# Check the data
gt <- extract.gt(vcf, element = "GT")
gt[c(2,6,18), 1:3]

# transpose the data
t(as.matrix(x))[c(2,6,18), 1:3]

# remove individual 7Tt182 as low coveerage
x1 <- x[indNames(x) != "7Tt182"]

# add population to samples
pop(x1) <-as.factor(c(rep("ancient", 1),rep("NEAp", 10), rep("NWAp", 9), rep("NEPp", 11),rep("NEAc", 10), rep("NWAc", 7), rep("NEPc", 9)))

popNames(x1)

## Heat map of genotype
pdf('paralle_all_ancient_noMiss_no7Tt182.pdf',width=17,height=6)

glPlot(x1,col=c("blue","green","red"),legend=F,yaxt='n',las=1)

# add population labels on the y axis
axis(2, at=9, labels="NEPc", pos=1,las=2)
axis(2, at=16, labels="NWAc", pos=1,las=2)
axis(2, at=26, labels="NEAc", pos=1,las=2)
axis(2, at=37, labels="NEPp", pos=1,las=2)
axis(2, at=46, labels="NWAp", pos=1,las=2)
axis(2, at=56, labels="NEAp", pos=1,las=2)
axis(2, at=57, labels="ancient", pos=1,las=2)

dev.off()

# Stats

# Plot a neighbour-joining (NJ) tree
library(ape)

tree <- nj(dist(as.matrix(x1)))
plot(tre, typ="fan", cex=0.7)

# run 100 bootstrap of the NJ tree
tree_boot <- boot.phylo(tree,B=100, x1, FUN = function(xx) nj(dist(xx)))

# Plot the tree with boostrap values
plot(tre, typ="fan", show.tip=FALSE)
nodelabels(tree_boot,frame = "n", cex = 0.8)

library(scales)
mycol<-c("black","#F4A582","#92C5DE","#B2182B","#2166AC","#D6604D","#4393C3")
tiplabels(pch=20,col=alpha(mycol[pop(x1)],0.9), cex=4)




