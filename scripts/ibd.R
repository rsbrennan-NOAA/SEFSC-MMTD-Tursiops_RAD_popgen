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
                     npc_selection = "manual") 
          # shirk suggests using 64 axes. do this

# https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12684


pc_dists_auto <- gen_dist((str_dos), dist_type = "pc", 
                     npc_selection = "auto",criticalpoint = 2.0234)

pc_plink <- gen_dist(
  dist_type = "plink",
  plink_file = "analysis/pop_structure/LDthin.dist",
  plink_id_file = "analysis/pop_structure/LDthin.dist.id",
  npc_selection = "manual"
)

gen_dist_corr(dist_x = pc_dists, dist_y = pc_plink, 
              metric_name_x = "pc_dists", metric_name_y = "plink")
gen_dist_corr(dist_x = pc_dists, dist_y = str_dist, 
              metric_name_x = "pc_dists", metric_name_y = "euclidian")
gen_dist_corr(dist_x = str_dist, dist_y = pc_plink, 
              metric_name_x = "euclidian", metric_name_y = "plink")
gen_dist_corr(dist_x = pc_plink, dist_y = pc_dists_auto, 
              metric_name_x = "pc_plink", metric_name_y = "pc_dists_auto")
gen_dist_corr(dist_x = pc_dists, dist_y = pc_dists_auto, 
              metric_name_x = "pc_dists", metric_name_y = "pc_dists_auto")

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

pops <- read.table("analysis/population_assignments_summary.txt", header=T)
dm <- left_join(out, pops, by=c("label" ="indiv"))
dm <- dm %>%
  separate(sixpop, into = c("Population", "region"), sep = "_", remove = FALSE) %>%
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
#pops <- read.table("analysis/pop_structure/sixpop_all.clust", header=F)
#colnames(pops) <-c("ids", "population")
#out$label <- gsub("X", "", out$label)

# this is changing coastal_Atl to Coastal_Atlantic
#pops <- pops %>%
#  mutate(population = ifelse(population == "Coastal_Atl", "Coastal_Atlantic", population)) %>%
#  separate(population, into = c("Population", "region"), sep = "_", remove = FALSE) %>%
#  mutate(region = ifelse(region == "Atl", "Atlantic", region))


#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
# idb four pop
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------


# add the region and pop:
#pops <- read.table("analysis/pop_structure/fourpop_all.clust", header=F)
#colnames(pops) <-c("ids", "population")
#out$label <- gsub("X", "", out$label)

merged_distances <- merge(pc_long, lc_long, 
                          by = c("indiv_1", "indiv_2"), 
                          all = FALSE)

id_to_pop <- setNames(pops$fourpop, pops$indiv)
merged_distances$population_1 <- id_to_pop[merged_distances$indiv_1]
merged_distances$population_2 <- id_to_pop[merged_distances$indiv_2]

head(merged_distances)

merged_distances <- merged_distances %>%
  mutate(comparison = case_when(
    population_1 == "Coastal_Atlantic" & population_2 == "Coastal_Atlantic" ~ "Coastal_Atlantic",
    population_1 == "Coastal_Gulf" & population_2 == "Coastal_Gulf" ~ "Coastal_Gulf",
    population_1 == "Intermediate" & population_2 == "Intermediate" ~ "Intermediate",
    population_1 == "Offshore" & population_2 == "Offshore" ~ "Offshore",
    TRUE ~ "inter-population"
  ))
head(merged_distances)

p2 <- ggplot(merged_distances, aes(x = leastcost_distance, y = pc_distance)) +
  geom_point(aes(color = comparison),alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE, linewidth=2) +
  theme_bw(base_size = 14) +
  labs(y = "Genetic Distance", 
       x = "Geographic Distance",
       title = "Isolation by distance: four populations") +
  facet_wrap(~comparison, ncol=3) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size=5)))

p2
ggsave("figures/ibd_fourpops.png", p2, h=7, w=9)



#########-------------------------------------------
#########-------------------------------------------
#########-------------------------------------------
#########-------------------------------------------
# six pops

id_to_pop <- setNames(pops$sixpop, pops$indiv)
merged_distances$population_1 <- id_to_pop[merged_distances$indiv_1]
merged_distances$population_2 <- id_to_pop[merged_distances$indiv_2]

id_to_region <- setNames(pops$region, pops$indiv)
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
       title = "Isolation by distance: Six populations") +
  facet_wrap(~comparison, ncol=3) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size=5)))

p2
ggsave("figures/ibd_sixpops.png", p2, h=7, w=9)


#-----------------------------------------
#-----------------------------------------
#-----------------------------------------
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

