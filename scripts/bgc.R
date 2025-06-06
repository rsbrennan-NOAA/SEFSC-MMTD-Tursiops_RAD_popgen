#bgc

library(bgchm)

library(triangulaR)
library(vcfR)
library(dplyr)
library(ggplot2)



# also useful here: https://github.com/zgompert/bgchm_test/blob/main/FitClineModel_Lycaeides.R

vcf <- read.vcfR( "analysis/variants/filtered.final_ids.vcf.gz", verbose = FALSE )
vcf

pops1 <- read.table("analysis/population_assignments_summary.txt", header=T)

pops <- data.frame(id = pops1$indiv,
                   pop = pops1$fourpop)

# identify ancestry informative markers (AIMS)
# Create a new vcfR object composed only of sites above the given allele frequency difference threshold
diff1 <- alleleFreqDiff(vcfR = vcf, 
                        pm = pops, 
                        p1 = "Coastal_Gulf", 
                        p2 = "Offshore", 
                        difference = 0.7)

# [1] "438 sites passed allele frequency difference threshold"

genos <- t(extract.gt(diff1))
genos <- gsub("0/0", 0,genos)
genos <- gsub("0/1", 1,genos)
genos <- gsub("1/1", 2,genos)

genos <- apply(genos, 2, as.numeric)
row.names(genos) <- colnames(extract.gt(diff1))
head(genos)

Coastal_Gulf <- pops$id[pops$pop == "Coastal_Gulf"]
Offshore <- pops$id[pops$pop == "Offshore"]
Intermediate <- pops$id[pops$pop == "Intermediate"]

GenP0 <- genos[row.names(genos) %in% Coastal_Gulf,]
GenP1 <- genos[row.names(genos) %in% Offshore,]
GenHybrids <- genos[row.names(genos) %in% Intermediate,]


#


# 0 (homozygote), 1 (heterozygote) and 2 (alternative homozygote).
# 1 row per indiv, 1 col per genotype
# set missing genotypes as NA
## estimate parental allele frequencies, uses analytical solution 
p_out<-est_p(G0=GenP0,
             G1=GenP1,
             model="genotype",ploidy="diploid",
             HMC=FALSE)



## estimate hybrid indexes, uses default HMC settings
## and uses point estimates (posterior medians) of allele frequencies
plin <- is.na(GenHybrids)  
plin[1:5,1:3]
plin <- matrix(ifelse(is.na(GenHybrids), 0, 2), 
               nrow=nrow(GenHybrids), ncol=ncol(GenHybrids))

plin[1:5,1:3]

h_out<-est_hi(Gx=GenHybrids,
              p0=p_out$p0[,1],
              p1=p_out$p1[,1],
              model="genotype",
              ploidy="mixed",
              pldat=plin)

## plot hybrid index estimates with 90% equal-tail probability intervals
## sorted by hybrid index, just a nice way to visualize that in this example we have
## few hybrids with intermediate hybrid indexes
# save these:

pdf(file="figures/bgc_hybrid_index.pdf", h=5, w=5)
plot(sort(h_out$hi[,1]),ylim=c(0,1),pch=19,xlab="Individual (sorted by HI)",ylab="Hybrid index (HI)")
segments(1:length(h_out$hi[,1]),
         h_out$hi[order(h_out$hi[,1]),3],
         1:60,
         h_out$hi[order(h_out$hi[,1]),4])
dev.off()

pdf(file="figures/bgc_hybrid_index_hmc.pdf", h=4, w=10)

#rstan::traceplot(h_out$hi_hmc)
rstan::traceplot(h_out$hi_hmc,pars=c("H[1]","H[10]","H[20]"))
dev.off()

## view the summary from rstan
h_out$hi_hmc

## fit a hierarchical genomic cline model for all  loci using the estimated
## hybrid indexes and parental allele frequencies (point estimates)
## use 4000 iterations and 2000 warmup to make sure we get a nice effective sample size
## note that this takes a while to run- 60 min or so.
gc_out<-est_genocl(Gx=GenHybrids,p0=p_out$p0[,1],
                   p1=p_out$p1[,1],H=h_out$hi[,1],
                   model="genotype",ploidy="mixed",
                   pldat=plin,hier=TRUE,
                   p_warmup = 0.5, n_iters=4000)

## how variable is introgression among loci? Lets look at the cline SDs
## these are related to the degree of coupling among loci overall
gc_out$SDc
gc_out$SDv

## examine a plot of the joint posterior distribution for the SDs
pp_plot(objs=gc_out,param1="sdv",param2="sdc",probs=c(0.5,0.75,0.95),colors="black",
        addPoints=TRUE,palpha=0.1,pdf=FALSE,pch=19)

## impose sum-to-zero constraint on log/logit scale
## not totally necessary, but this is mostly a good idea
sz_out<-sum2zero(hmc=gc_out$gencline_hmc,transform=TRUE,ci=0.90)

## plot genomic clines for the 51 loci, first without the sum-to-zero constraint
## then with it... these differ more for some data sets than others
gencline_plot(center=gc_out$center[,1],v=gc_out$gradient,pdf=FALSE)
gencline_plot(center=sz_out$center[,1],v=sz_out$gradient,pdf=FALSE)
# they look to be very similar

# loci with credible deviations from genome-average ancestry (90% CIs for cline gradient of center not overlapping null expectations).

## summarize loci with credible deviations from genome-average gradients, 
  # here the focus is
  # specifically on steep clines indicative of loci introgressing less than the average
# from the paper:
# Genomic cline slopes greater than 1 indicate a steeper cline than the admixture gradient,
# clines less than 1 indicate a shallower cline.
# centers greater than 0.5 indicate an overall excess of source 0 ancestry, 
# whereas centers less than 0.5 indicate an excess of source 1 ancestry.

# the cols of gradient are 50%, 5%, and 95% CI. 
  # So check that the 5% is greater than 1 for loci that are steeper than average
    # these are introgressing less than expected.
  # less than 1 is shallower than expected, higher introgression.
which(sz_out$gradient[,2] > 1) ## index for loci with credibly steep clines
sum(sz_out$gradient[,2] > 1) ## number of loci with credibly steep clines
# 26
## number of loci with credibly shallow clines, now 95% instead of 5%
sum(sz_out$gradient[,3] < 1) 
# 39


# what about the center:
which(sz_out$center[,2] > 1) ## index for loci with credibly steep clines

which(sz_out$center[,2] > .5)
which(sz_out$center[,3] < .5)
sum(sz_out$center[,2] > .5)
#33
sum(sz_out$center[,3] < .5)
# 50

## last, lets look at interspecific ancestry for the same data set, this can
## be especially informative about the types of hybrids present
q_out<-est_Q(Gx=GenHybrids,
             p0=p_out$p0[,1],p1=p_out$p1[,1],
             model="genotype",ploidy="mixed",
             pldat=plin,
             n_chains = 4,
             n_iters = 3000,
             p_warmup = 0.5)
q_out$Q_hmc
rstan::traceplot(q_out$Q_hmc)


## plot the results
#pdf(file="figures/bgc_trianglePlot.pdf", h=4, w=4)

#tri_plot(hi=q_out$hi[,1],Q10=q_out$Q10[,1],pdf=FALSE,pch=19)
## note that some individuals appear to be likely backcrosses (close to the outer lines of the triangles)
## but the individals with intermediate hybrid indexes are clearly not F1s but rather late generation hybrids
dev.off()


# bc the above are slow to run, save and load for later
save.image(file = "analysis/bgc.RData")

# load environment
load(file = "analysis/bgc.RData")








#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# analyze the actual results
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#--------------
## cline plots
#--------------

o<-sz_out

library(ggplot2)
library(scales) # For alpha 

# makehybrid index values
h <- seq(0, 1, 0.01)

# make empty df to store curves
plot_data <- data.frame()

# non-significant 
nonsig <- (o$gradient[,2] <= 1)
nonsig_center <- o$cent[nonsig,1]
nonsig_v <- o$gradient[nonsig,1]
nonsig_u <- log(nonsig_center/(1 - nonsig_center)) * nonsig_v

# create non-sig curves df
for (i in 1:length(nonsig_v)) {
  phi <- (h^nonsig_v[i])/(h^nonsig_v[i] + (1 - h)^nonsig_v[i] * exp(nonsig_u[i]))
  temp_df <- data.frame(h = h, phi = phi, 
                        group = paste0("nonsig_", i), 
                        type = "nonsig")
  plot_data <- rbind(plot_data, temp_df)
}

#---------
# shallow curves
shallow <- (o$gradient[,3] < 1)
shallow_center <- o$cent[shallow,1]
shallow_v <- o$gradient[shallow,1]
shallow_u <- log(shallow_center/(1 - shallow_center)) * shallow_v

for (i in 1:length(shallow_v)) {
  phi <- (h^shallow_v[i])/(h^shallow_v[i] + (1 - h)^shallow_v[i] * exp(shallow_u[i]))
  temp_df <- data.frame(h = h, phi = phi, 
                        group = paste0("shallow_", i), 
                        type = "shallow")
  plot_data <- rbind(plot_data, temp_df)
}

#---------------
#  steep curves
sig <- (o$gradient[,2] > 1)
sig_center <- o$cent[sig,1]
sig_v <- o$gradient[sig,1]
sig_u <- log(sig_center/(1 - sig_center)) * sig_v

for (i in 1:length(sig_v)) {
  phi <- (h^sig_v[i])/(h^sig_v[i] + (1 - h)^sig_v[i] * exp(sig_u[i]))
  temp_df <- data.frame(h = h, phi = phi, 
                        group = paste0("sig_", i), 
                        type = "sig")
  plot_data <- rbind(plot_data, temp_df)
}

# Add diagonal line data, just for plotting
diag_line <- data.frame(h = c(0, 1), phi = c(0, 1), group = "diagonal", 
                        type = "diagonal")
plot_data <- rbind(plot_data, diag_line)

head(plot_data)
# Create the plot
p <- ggplot() +
  # Plot the curves by type with appropriate colors and line widths
  # We need to use aes(color = ...) to get the legend to work
  geom_line(data = subset(plot_data, type == "nonsig"), 
            aes(x = h, y = phi, group = group, color = "Neutral"), 
            linewidth = 0.5, alpha = 0.5) +
  geom_line(data = subset(plot_data, type == "shallow"), 
            aes(x = h, y = phi, group = group, color = "Shallow"), 
            linewidth = 1, alpha = 0.5) +
  geom_line(data = subset(plot_data, type == "sig"), 
            aes(x = h, y = phi, group = group, color = "Steep"), 
            linewidth = 1, alpha = 0.5) +
  # Add the diagonal line
  geom_line(data = subset(plot_data, type == "diagonal"), 
            aes(x = h, y = phi), 
            linetype = 2, linewidth = 2, color = "black") +
  # Set custom colors for the legend
  scale_color_manual(name = NULL,
                     values = c("Neutral" = "grey45", 
                                "Shallow" = "skyblue3", 
                                "Steep" = "firebrick3")) +
  guides(color = guide_legend(override.aes = list(linewidth = 2, alpha=1))) +
  # Set plot labels and limits
  labs(x = "Hybrid index", 
       y = "Ancestry probability") +
  # Match the original plot limits
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  # Control theme elements to match original appearance
  theme_classic(base_size = 14) +
  theme(

    legend.position = "bottom",

  )

p
#
ggsave(file="figures/genomic_clines.png", p, 
       h=6.5, w=7)


#---------------------------------------------------
#---------------------------------------------------
#---------------------------------------------------
#---------------------------------------------------
# plot gradients (V) CI and result:

# gompert was ordering by linkage group
#sig<-1 + (o$gradient[order(snpLgSz[,1]),2] > 1)

split_cols <- strsplit(colnames(genos), "_")
# Extract the accession numbers (first element of each split)
chroms <- sapply(split_cols, function(x) paste0(x[1], "_", x[2]))
SNP <- colnames(genos)
# ggplot
plot_data <- data.frame(
  SNP_id <- SNP,
  CHR <- chroms,
  SNP_number = 1:nrow(o$gradient),
  log_gradient = log(o$gradient[,1]),
  lower_ci = log(o$gradient[,2]),
  upper_ci = log(o$gradient[,3])
)

colnames(plot_data) <- c("SNP", "CHR", "SNP_number", "log_gradient", "lower_ci", "upper_ci")

plot_data$category <- "Neutral"
plot_data$category[o$gradient[,2] > 1] <- "Steep"  # significant
plot_data$category[o$gradient[,3] < 1] <- "Shallow" # shallow

# make a factor for CHR 
plot_data$CHR <- factor(plot_data$CHR, levels = unique(plot_data$CHR))

# colors for chromosomes
chr_colors <- rep(c("black", "grey1"), length.out = length(levels(plot_data$CHR)))
names(chr_colors) <- levels(plot_data$CHR)

# order by chr then position:


# Split data by category
neutral_data <- subset(plot_data, category == "Neutral")
shallow_data <- subset(plot_data, category == "Shallow")
steep_data <- subset(plot_data, category == "Steep")

chr_bounds <- plot_data %>%
  # We need to join with the chr_colors information
  mutate(color = chr_colors[as.character(CHR)]) %>%
  # Only keep chromosomes with black color
  filter(color == "black") %>%
  group_by(CHR) %>%
  summarize(min_snp = min(SNP_number),
            max_snp = max(SNP_number)) %>%
  mutate(ymin = -Inf, ymax = Inf)

# for chr labels
chr_levels <- levels(plot_data$CHR)
n_chrs <- length(chr_levels)

# Create new chromosome labels (1:21 plus X)
new_chr_labels <- c(as.character(1:(n_chrs-1)), "X")

# Calculate midpoints for each chromosome for label placement
chr_midpoints <- plot_data %>%
  group_by(CHR) %>%
  summarize(mid_pos = mean(c(min(SNP_number), max(SNP_number)))) %>%
  arrange(match(CHR, chr_levels))  # Maintain original order

# Assign new labels to the midpoints dataframe
chr_midpoints$new_label <- new_chr_labels

# Create the plot with 3 geom_point layers - one for each category
p2 <- ggplot() +
  geom_rect(data = chr_bounds,
            aes(xmin = min_snp, xmax = max_snp,
                ymin = ymin, ymax = ymax,
                fill = CHR),
            alpha = 1, fill="grey80") +
  # Add points for each category with alternating chromosome colors
  # Add confidence interval segments colored by category
  geom_segment(data = neutral_data,
               aes(x = SNP_number, y = lower_ci, 
                   xend = SNP_number, yend = upper_ci),
               color = "grey30", linewidth = 0.3) +
  geom_segment(data = shallow_data,
               aes(x = SNP_number, y = lower_ci, 
                   xend = SNP_number, yend = upper_ci),
               color = "skyblue3", linewidth = 0.3) +
  geom_segment(data = steep_data,
               aes(x = SNP_number, y = lower_ci, 
                   xend = SNP_number, yend = upper_ci),
               color = "firebrick3", linewidth = 0.3) +
  geom_point(data = neutral_data, 
             aes(x = SNP_number, y = log_gradient, color = CHR),
             size = 3, shape = 19) +
  geom_point(data = shallow_data, 
             aes(x = SNP_number, y = log_gradient),
             size = 3.5, shape = 21, fill = "skyblue3", color="black") +
  geom_point(data = steep_data, 
             aes(x = SNP_number, y = log_gradient, color = CHR),
             size = 3.5, shape = 21,fill = "firebrick3", color="black") +
#  geom_point(data = shallow_data, 
#             aes(x = SNP_number, y = log_gradient),
#             color = "skyblue3",
#             size = 3, shape = 21, stroke =2) +
#  geom_point(data = steep_data, 
#             aes(x = SNP_number, y = log_gradient),
#             color = "firebrick3",
#             size = 3, shape = 21, stroke =2) +
  # Add horizontal line at y=0
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
  
  # Set point colors for chromosomes
  scale_color_manual(values = chr_colors, guide = "none") +  # Hide chromosome colors from legend
  
  # Add chromosome labels
  scale_x_continuous(
    breaks = chr_midpoints$mid_pos,
    labels = chr_midpoints$new_label,
    expand = c(0.01, 0.01)
  ) +
  
  # Create a dummy aesthetic for the legend
  #annotate("point", x = -Inf, y = -Inf, color = "grey30", size = 3, shape = 19) +
  #annotate("point", x = -Inf, y = -Inf, color = "skyblue3", size = 3, shape = 19) +
  #annotate("point", x = -Inf, y = -Inf, color = "firebrick3", size = 3, shape = 19) +
  
  labs(x = "Chromosome", 
       y = "Log gradient (v)") +
  theme_classic(base_size = 14) +
  theme(
    #legend.position = "bottom",
  )
p2

# Add manual legend
#p2 <- p2 + guides(color = guide_legend(
#  override.aes = list(
#    color = c("grey30", "skyblue3", "firebrick3"),
#    shape = c(19, 19, 19),
#    size = c(3, 3, 3)
#  ),
#  labels = c("Neutral", "Shallow", "Steep")
#))

#p2


ggsave(filename="figures/bgc_p1.pdf",p, h=4, w=3.5)
ggsave(filename="figures/bgc_p2.pdf",p2, h=4, w=7)

ggsave(filename="figures/bgc_combined.pdf",
       ggarrange(p,p2, widths=c(0.33,1), common.legend=T),
       h=3.5, w=11)


#---------------------------------------------------------
#---------------------------------------------------------
#---------------------------------------------------------
# analyze the outliers, their locations, expectations
#---------------------------------------------------------
#---------------------------------------------------------
#---------------------------------------------------------

head(plot_data)

# add snp id:
snp_ids <- row.names(extract.gt(diff1))
CHR <- sapply(strsplit(snp_ids, "_"), function(x) paste(x[1], x[2], sep="_"))
SNP <- sapply(strsplit(snp_ids, "_"), function(x) x[3])
nrow(plot_data)
snp_df <- data.frame(CHR2 = CHR, POS = SNP, snp_id = snp_ids)
alldat <- cbind(snp_df, plot_data)

head(alldat)

library(dplyr)
library(tidyr)

result <- alldat %>%
  group_by(CHR, category) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = category, values_from = count, values_fill = 0) %>%
  mutate(Total = rowSums(across(where(is.numeric))))

result$prop_steep <- round(result$Steep/result$Total, 2)
result$prop_shallow <- round(result$Shallow/result$Total,2)

as.data.frame(result)

# 8/21 loci on NC_047055.1  are steep outliers. 
# NC_047055.1 is the x chromosome. this is similar to pattern in gompert paper and many others
sum(result$Steep)

# want to look at chromosome legnth and number of outliers
# 

chrlen <- read.table("analysis/chr_length.txt", header=F)
colnames(chrlen) <- c("CHR", "Length")
chrlen <- chrlen[1:nrow(chrlen),1:2]


result <- merge(chrlen, result)

df_long <- pivot_longer(
  result,
  cols = c(Neutral, Steep, Shallow, Total),
  names_to = "Variable",
  values_to = "Value"
)

# Create the faceted plot
ggplot(df_long, aes(x = (Length)/1000000, y = Value)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ Variable, scales = "free_y") +
  labs(
    title = "CHR length vs # loci",
    x = "Chromosome Length (mb)",
    y = "Count"
  ) +
  theme_bw() 

ggsave("figures/bgc_length_vs_loci.png", h=4,w=5)


plt1 <- ggplot(result, aes(x = Total, y = Steep)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Total number of loci",
    y = "Number of steep outliers",
    title="Steep outliers"
  ) +
  theme_bw() 


plt2 <- ggplot(result, aes(x = Total, y = Shallow)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Total number of loci",
    y = "Number of shallow outliers",
    title="Shallow outliers"
  ) +
  theme_bw() 

ggsave("figures/bgc_sig_vs_total.png", ggarrange(plt1, plt2), h=3,w=6)



#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# need to run randomization test for enrichment of outliers on certain chromosomes
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

sig_v<-o$gradient[,2] > 1

result$null_steep <- NA
result$pval_steep <- NA

for(i in 1:nrow(result)){
  tmp_chr <- result$CHR[i]
  tmp_total <- result$Total[i]
  tmp_steep_obs <- result$Steep[i]
  null_v<-rep(NA,1000)
  for(j in 1:1000){
    null_v[j]<-sum(sample(as.numeric(sig_v),tmp_total,replace=FALSE))
  }
  result$pval_steep[i] <- mean(null_v >= tmp_steep_obs)
  result$null_steep[i] <- mean(null_v)
}


# shallow outliers:


sig_v<-o$gradient[,3] < 1

result$null_shallow <- NA
result$pval_shallow <- NA

for(i in 1:nrow(result)){
  tmp_chr <- result$CHR[i]
  tmp_total <- result$Total[i]
  tmp_steep_obs <- result$Shallow [i]
  null_v<-rep(NA,1000)
  for(j in 1:1000){
    null_v[j]<-sum(sample(as.numeric(sig_v),tmp_total,replace=FALSE))
  }
  result$pval_shallow[i] <- mean(null_v >= tmp_steep_obs)
  result$null_shallow[i] <- mean(null_v)
}

as.data.frame(result)

library(googlesheets4)
library(googledrive)

# need to only do this authenticatino once
#gs4_auth()
#drive_auth()

# Write  dataframe to a new google sjeet
sheet_url <- gs4_create("bgc_summary", sheets = list(data = result))

sheet_url

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# cline centers

# https://github.com/zgompert/bgchm_test/blob/18ba0853ce36dcc9cf8e69da0f7a8f3b3d33f4ef/FitClineModel_Lycaeides.R



# think about cline center vs slope. 


