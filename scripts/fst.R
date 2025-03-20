# pairwise fst


library(snpR)
library(vcfR)
library(dplyr)

dat <- snpR::read_vcf("analysis/filtered.final_ids.vcf.gz")

pops <- read.table("analysis/pop_structure/sixpop_all.clust", header=F)

colnames(pops) <- c("IDs", "population")

# add meta data information:
## population
ids <- data.frame(IDs = colnames(dat))
result <- ids %>%
  left_join(pops, by = c("IDs"))

sample_meta <- data.frame(pop = result$population)
## order the population
#sample_meta$pop <- factor(sample_meta$pop, levels=c("GA", "HP", "BC", "PC", "TR")) 

# assign meta data to dat
sample.meta(dat) <- sample_meta

# calculate fst between the populations
my.dat <- calc_pairwise_fst(dat, facets="pop", method = "WC", boot = 500)
# the bootstrapping here is actually permutation, mixing pop assignments
fst_pvals <- get.snpR.stats(my.dat, facets = "pop", stats = "fst")$fst.matrix$pop$p






# convert p-values to long:

dat_fst_pvals2 <- as.matrix(fst_pvals[,2:ncol(fst_pvals)])
rownames(dat_fst_pvals2) <- fst_pvals$p1

melt_fst_pval<- as_tibble(dat_fst_pvals2, rownames = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value")%>%
  mutate(across(c(Var1, Var2), ~str_replace(., "Dry_Tortuga", "Dry Tortuga")))

# remove na that isn't diag
melt_fst_pval <- subset(melt_fst_pval, !(Var1 != Var2 & is.na(value)))


# save p-values
write.csv(melt_fst_pval, file="fst_pvalues.csv", quote=F, row.names = F)


# convert fst to long
dat_fst <- get.snpR.stats(my.dat, facets = "pop", stats = "fst")$fst.matrix$pop$fst

dat_fst2 <- as.matrix(dat_fst[,2:ncol(dat_fst)])
rownames(dat_fst2) <- dat_fst$p1

melt_fst <- as_tibble(dat_fst2, rownames = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value")%>%
  mutate(across(c(Var1, Var2), ~str_replace(., "Dry_Tortuga", "Dry Tortuga")))

melt_fst <- subset(melt_fst, !(Var1 != Var2 & is.na(value)))

# save fst
write.csv(melt_fst, file="fst.csv", quote=F, row.names = F)
