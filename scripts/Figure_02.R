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

setwd("~/projects/Tursiops_RAD_popgen/analysis/pop_structure/treemix/fourpop/bootstrap")
#folder with all TreeMix outputs from the final runs
source("~/projects/Tursiops_RAD_popgen/scripts/TreeMix_functions.R") #path to required functions for this analysis

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

dev.off()


#----------------------------------------------------
# bgc:

setwd("~/projects/Tursiops_RAD_popgen/")


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
            linetype = 2, linewidth = 1, color = "black") +
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
  theme_classic(base_size = 12) +
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
             size = 2, shape = 21, fill="grey50", color="black") +
  geom_point(data = shallow_data, 
             aes(x = SNP_number, y = log_gradient),
             size = 2.5, shape = 21, fill = "skyblue3", color="black") +
  geom_point(data = steep_data, 
             aes(x = SNP_number, y = log_gradient, color = CHR),
             size = 2.5, shape = 21,fill = "firebrick3", color="black") +
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
  theme_classic(base_size = 12) +
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


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# triangle plot

setwd("~/projects/Tursiops_RAD_popgen/")

##-------------------------------------------------------------------------------
##-------------------------------------------------------------------------------
##-------------------------------------------------------------------------------
# triangle plot
library(triangulaR)
library(vcfR)

vcf <- read.vcfR( "analysis/variants/filtered.final_ids.vcf.gz", verbose = FALSE )
vcf

pops1 <- read.table("analysis/population_assignments_summary.txt", header=T)

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



pops <- read.table("analysis/population_assignments_summary.txt", header=T)
#hybs <- pops$indiv[pops$offshore_putative_hybrids == TRUE]


pops <- read.table("analysis/population_assignments_summary.txt", header=T)

pops2 <- read.table("analysis/pop_structure/newhybrids/newhybrids_posteriors.txt", header=T)

#hybs <- pops$indiv[pops$offshore_putative_hybrids == TRUE]
newhybrids <- pops2[pops2$sig == "High_Confidence" & pops2$highest_category != "Pure1" & pops2$highest_category != "Pure2" ,]

F2 <- newhybrids$indiv[newhybrids$highest_category == "F2"]
Backcross2  <- newhybrids$indiv[newhybrids$highest_category == "Backcross2"]
Backcross1 <- newhybrids$indiv[newhybrids$highest_category == "Backcross1"]

hi.het1$pop2 <- hi.het1$pop
hi.het1$pop2[hi.het1$id %in% F2] <- "F2"
hi.het1$pop2[hi.het1$id %in% Backcross2] <- "Backcross2"
hi.het1$pop2[hi.het1$id %in% Backcross1] <- "Backcross1"

df <- merge(hi.het1, pops,  by.x="id", by.y="indiv")
df$pop <- as.factor(df$sixpop)
hybsOnly <- df[which(df$pop2 == "F2" | df$pop2 == "Backcross1" | df$pop2 == "Backcross2" ),]
#hybsOnly$pop2 <- as.factor(hybsOnly$pop2)



fill_colors <- c(
  "Coastal_Atlantic" = "#A6DDF0",
  "Coastal_Gulf" = "#276FBF",
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
                 color= pop2), 
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
ggsave(filename="figures/triangle_hybrids.pdf",t2, w=3.6, h=3.5)



#plt1<- ggarrange(p,p, t2, common.legend=T, nrow=1)
#ggsave(filename="figures/bgc_combined.pdf",
#       ggarrange(plt1,p2, widths=c(0.33, 0.33,1), common.legend=T, nrow=2),
#       h=6, w=9)


plt1 <- ggarrange(t2,t2,p, common.legend=TRUE, nrow=1, labels="AUTO")
plt1
final_plot <- ggarrange(plt1, p2, heights=c(1.7, 1), common.legend=TRUE, nrow=2, labels=c("", "D"))
final_plot

# Save the final combined plot
ggsave(filename="figures/bgc_combined.pdf",
       final_plot,
       height=5, width=7)





