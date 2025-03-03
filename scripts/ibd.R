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
out_tree <- nj(as.matrix(pc_plink))
out_tree <- nj(as.matrix(pc_dists))
out <- as_tibble(nj(as.matrix(pc_dists)))
#out$label <- gsub("X", "", out$label)
ggtree(out_tree) + 
  theme_tree()

#p <- ggtree(out_tree,layout="daylight") + theme_tree()
#p <- ggtree(out_tree, layout="unrooted") + theme_tree()
p <- ggtree(out_tree, layout="ape") + theme_tree()
p

# shirk paper says use 64 pc model when we have sub structure
# only a need a few for initial structure, but these others increase ability
 # to detect futher structure. 
# I think we have substantial sub structure, so we want to use 64

# color the plot










