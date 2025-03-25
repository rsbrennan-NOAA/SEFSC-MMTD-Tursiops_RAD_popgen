# pairwise fst


library(snpR)
library(vcfR)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

dat <- snpR::read_vcf("analysis/filtered.final_ids.vcf.gz")

pops <- read.table("analysis/population_assignments_summary.txt", header=T)

# add meta data information:
## population
ids <- data.frame(IDs = names(dat))
result <- ids %>%
  left_join(pops, by = c("IDs" = "indiv"))


sample_meta <- data.frame(pop = result$sixpop)
## order the population
#sample_meta$pop <- factor(sample_meta$pop, levels=c("GA", "HP", "BC", "PC", "TR")) 

# assign meta data to dat
sample.meta(dat) <- sample_meta

# calculate fst between the populations
my.dat <- calc_pairwise_fst(dat, facets="pop", method = "WC", boot = 500)
# the bootstrapping here is actually permutation, mixing pop assignments
fst_pvals <- get.snpR.stats(my.dat, facets = "pop", stats = "fst")$fst.matrix$pop$p
dat_fst <- get.snpR.stats(my.dat, facets = "pop", stats = "fst")$fst.matrix$pop$fst
dat_fst

# convert p-values to long:

dat_fst_pvals2 <- as.matrix(fst_pvals[,2:ncol(fst_pvals)])
rownames(dat_fst_pvals2) <- fst_pvals$p1

melt_fst_pval<- as_tibble(dat_fst_pvals2, rownames = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value") 

# remove na that isn't diag
melt_fst_pval <- subset(melt_fst_pval, !(Var1 != Var2 & is.na(value)))


# save p-values
write.csv(melt_fst_pval, file="analysis/fst_pvalues.csv", quote=F, row.names = F)


# convert fst to long
dat_fst2 <- as.matrix(dat_fst[,2:ncol(dat_fst)])
rownames(dat_fst2) <- dat_fst$p1

melt_fst <- as_tibble(dat_fst2, rownames = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value")

melt_fst <- subset(melt_fst, !(Var1 != Var2 & is.na(value)))

# save fst
write.csv(melt_fst, file="analysis/fst.csv", quote=F, row.names = F)


##### figure ------------------------------------------------
# flip order for pvals, to put on other triangle:
melt_fst_pval2 <- melt_fst_pval[, c("Var2", "Var1", "value")]
colnames(melt_fst_pval2) <- c("Var1", "Var2", "value")

melt_fst$group <- c("fst")
melt_fst_pval2$group <- c("pval")
fst_plot_dat <- rbind(melt_fst, melt_fst_pval2)

fst_plot_dat$group[melt_fst$Var1 == melt_fst$Var2] <- "diagonal"

fst_plot_dat$Var1 <- as.character(fst_plot_dat$Var1 )
fst_plot_dat$Var2 <- as.character(fst_plot_dat$Var2 )

as.data.frame(fst_plot_dat)

#fst_plot_dat$Var1 <- factor(fst_plot_dat$Var1, c("Atlantic", "Dry Tortuga", "N. Gulf", "W. Gulf"))
#fst_plot_dat$Var2 <- factor(fst_plot_dat$Var2, c("Atlantic", "Dry Tortuga", "N. Gulf", "W. Gulf"))


ggplot(fst_plot_dat, aes()) +
  geom_tile(aes(x = Var1, y = Var2, fill = value),
            size = 0.1
  ) 

fst_plot_dat2 <- fst_plot_dat

# add other diags
fst_plot_dat2 <- rbind(fst_plot_dat2, data.frame(
  Var1 = c("Coastal_Atlantic", "Offshore_Gulf"),
  Var2 = c("Coastal_Atlantic", "Offshore_Gulf"),
  value = c(NA, NA),
  group= c("diagonal", "diagonal")
)
)


# drop p-vals
fst_plot_dat3 <- subset(fst_plot_dat2, !group %in% c( "pval"))
fst_plot_dat4 <- subset(fst_plot_dat2, group %in% c( "pval"))
fst_plot_diags <- subset(fst_plot_dat2, group %in% c( "diagonal"))


p_fst <- ggplot(fst_plot_dat3, aes(y = Var2, x = Var1, fill = value)) +
  geom_tile(size=0.2, color="black") +
  geom_text(aes(label = round(value, 3)), na.rm = TRUE, size=3) +
  scale_fill_gradient(limits = range(fst_plot_dat3$value, na.rm = TRUE),
                      low = "#FFCCCC", 
                      high = "firebrick3", 
                      na.value = "white") +  
  theme_bw(base_size = 13) +
  labs(x = NULL, y = NULL, fill =  bquote(F[ST])) +
  geom_tile(data = fst_plot_diags, 
            aes(x = Var1, y = Var2),
            fill="grey",
            size=0.2, color="black") +
  geom_tile(data = fst_plot_dat4, 
            aes(x = Var1, y = Var2),
            fill="white",
            size=0.2, color="black") +
  geom_text(data = fst_plot_dat4, 
            aes(label = round(value, 3)),
            size=3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.position = c(0.95, 0.05),
    #legend.justification = c(1, 0),
    legend.background = element_blank(),  # Remove legend background
    legend.margin = margin(0, 0, 0, 0),  # Remove legend margin
    legend.box.margin = margin(0, 0, 0, -10),  # negative value moves legend left
    legend.key.size = unit(0.5, "lines"),  #  reduced legend key size
    legend.title = element_text(size = 9, margin = margin(b = 0)),  # Reduce bottom margin of title
    legend.text = element_blank(),
    #legend.text = element_text(size = 7, margin = margin(l = 1, r = 0)),  # Reduced legend text size
    legend.spacing.y = unit(0.05, "cm"),  # Reduced spacing between legend elements
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p_fst

ggsave(file="figures/fst.png", h=4, w=4)
ggsave(file="figures/fst.pdf", h=4, w=4)








#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# fourpop

sample_meta <- data.frame(pop = result$fourpop)
## order the population
#sample_meta$pop <- factor(sample_meta$pop, levels=c("GA", "HP", "BC", "PC", "TR")) 

# assign meta data to dat
sample.meta(dat) <- sample_meta

# calculate fst between the populations
my.dat <- calc_pairwise_fst(dat, facets="pop", method = "WC", boot = 500)
# the bootstrapping here is actually permutation, mixing pop assignments
fst_pvals <- get.snpR.stats(my.dat, facets = "pop", stats = "fst")$fst.matrix$pop$p
dat_fst <- get.snpR.stats(my.dat, facets = "pop", stats = "fst")$fst.matrix$pop$fst
dat_fst

# convert p-values to long:

dat_fst_pvals2 <- as.matrix(fst_pvals[,2:ncol(fst_pvals)])
rownames(dat_fst_pvals2) <- fst_pvals$p1

melt_fst_pval<- as_tibble(dat_fst_pvals2, rownames = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value") 

# remove na that isn't diag
melt_fst_pval <- subset(melt_fst_pval, !(Var1 != Var2 & is.na(value)))


# save p-values
write.csv(melt_fst_pval, file="analysis/fst_fourpops_pvalues.csv", quote=F, row.names = F)


# convert fst to long
dat_fst2 <- as.matrix(dat_fst[,2:ncol(dat_fst)])
rownames(dat_fst2) <- dat_fst$p1

melt_fst <- as_tibble(dat_fst2, rownames = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value")

melt_fst <- subset(melt_fst, !(Var1 != Var2 & is.na(value)))

# save fst
write.csv(melt_fst, file="analysis/fst_fourpops.csv", quote=F, row.names = F)


##### figure ------------------------------------------------
# flip order for pvals, to put on other triangle:
melt_fst_pval2 <- melt_fst_pval[, c("Var2", "Var1", "value")]
colnames(melt_fst_pval2) <- c("Var1", "Var2", "value")

melt_fst$group <- c("fst")
melt_fst_pval2$group <- c("pval")
fst_plot_dat <- rbind(melt_fst, melt_fst_pval2)

fst_plot_dat$group[melt_fst$Var1 == melt_fst$Var2] <- "diagonal"

fst_plot_dat$Var1 <- as.character(fst_plot_dat$Var1 )
fst_plot_dat$Var2 <- as.character(fst_plot_dat$Var2 )

as.data.frame(fst_plot_dat)

#fst_plot_dat$Var1 <- factor(fst_plot_dat$Var1, c("Atlantic", "Dry Tortuga", "N. Gulf", "W. Gulf"))
#fst_plot_dat$Var2 <- factor(fst_plot_dat$Var2, c("Atlantic", "Dry Tortuga", "N. Gulf", "W. Gulf"))


ggplot(fst_plot_dat, aes()) +
  geom_tile(aes(x = Var1, y = Var2, fill = value),
            size = 0.1
  ) 

fst_plot_dat2 <- fst_plot_dat

# add other diags
fst_plot_dat2 <- rbind(fst_plot_dat2, data.frame(
  Var1 = c("Coastal_Atlantic", "Offshore"),
  Var2 = c("Coastal_Atlantic", "Offshore"),
  value = c(NA, NA),
  group= c("diagonal", "diagonal")
)
)


# drop p-vals
fst_plot_dat3 <- subset(fst_plot_dat2, !group %in% c( "pval"))
fst_plot_dat4 <- subset(fst_plot_dat2, group %in% c( "pval"))
fst_plot_diags <- subset(fst_plot_dat2, group %in% c( "diagonal"))


p_fst <- ggplot(fst_plot_dat3, aes(y = Var2, x = Var1, fill = value)) +
  geom_tile(size=0.2, color="black") +
  geom_text(aes(label = round(value, 3)), na.rm = TRUE, size=3) +
  scale_fill_gradient(limits = range(fst_plot_dat3$value, na.rm = TRUE),
                      low = "#FFCCCC", 
                      high = "firebrick3", 
                      na.value = "white") +  
  theme_bw(base_size = 13) +
  labs(x = NULL, y = NULL, fill =  bquote(F[ST])) +
  geom_tile(data = fst_plot_diags, 
            aes(x = Var1, y = Var2),
            fill="grey",
            size=0.2, color="black") +
  geom_tile(data = fst_plot_dat4, 
            aes(x = Var1, y = Var2),
            fill="white",
            size=0.2, color="black") +
  geom_text(data = fst_plot_dat4, 
            aes(label = round(value, 3)),
            size=3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.position = c(0.95, 0.05),
    #legend.justification = c(1, 0),
    legend.background = element_blank(),  # Remove legend background
    legend.margin = margin(0, 0, 0, 0),  # Remove legend margin
    legend.box.margin = margin(0, 0, 0, -10),  # negative value moves legend left
    legend.key.size = unit(0.5, "lines"),  #  reduced legend key size
    legend.title = element_text(size = 9, margin = margin(b = 0)),  # Reduce bottom margin of title
    legend.text = element_blank(),
    #legend.text = element_text(size = 7, margin = margin(l = 1, r = 0)),  # Reduced legend text size
    legend.spacing.y = unit(0.05, "cm"),  # Reduced spacing between legend elements
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p_fst

ggsave(file="figures/fst_fourpop.png", h=4, w=4)
ggsave(file="figures/fst_fourpop.pdf", h=4, w=4)

