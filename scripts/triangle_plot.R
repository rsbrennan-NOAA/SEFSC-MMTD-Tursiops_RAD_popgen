# triangle plots

# https://omys-omics.github.io/triangulaR/index.html

library(triangulaR)
library(vcfR)

vcf <- read.vcfR( "analysis/variants/filtered.final_ids.vcf.gz", verbose = FALSE )
vcf

#***** Object of Class vcfR
#345 samples
#23 CHROMs
#7,819 variants
#Object size: 90.3 Mb
#0 percent missing data

# make a popmap
#>     id pop
#> 1  i55  P1
#> 2 i159  P1
#> 3 i245  P1
#> 4 i246  P1
#> 5 i264  P1
#> 6 i526  P1
#> 
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
#[1] "438 sites passed allele frequency difference threshold"
write.table(data.frame(snps=row.names(extract.gt(diff1))),
            quote = F,
            col.names = F,
            file="analysis/snp_set_gulfVSoffshore.txt",
            sep="\t",
            row.names=F
            )


diff2 <- alleleFreqDiff(vcfR = vcf, 
                        pm = pops, 
                        p1 = "Coastal_Atlantic", 
                        p2 = "Offshore", 
                        difference = 0.7)
#[1] "450 sites passed allele frequency difference threshold"

# Calculate hybrid index and heterozygosity for each sample. 
  # Values are returned in a data.frame
hi.het1 <- hybridIndex(vcfR = diff1, 
                      pm = pops, 
                      p1 = "Coastal_Gulf", p2 = "Offshore")

hi.het2 <- hybridIndex(vcfR = diff2, 
                      pm = pops, 
                      p1 = "Coastal_Atlantic", p2 = "Offshore")

plot(hi.het1$hybrid.index, hi.het2$hybrid.index)
plot(hi.het1$heterozygosity, hi.het2$heterozygosity)
# they are largely consistent

# Generate colors (or leave blank to use default)
#cols <- c("#af8dc3", "#7fbf7b", "#bababa", "#878787", "#762a83", "#1b7837")
# View triangle plot
#triangle.plot(hi.het, colors = cols)
cols <- c("#1B9E77FF", "#D95F02FF", "#7570B3FF", "#E7298AFF", "#66A61EFF", "#E6AB02FF", "#A6761DFF", "#666666FF")

t1 <- triangle.plot(hi.het1) + ggtitle("Coastal Gulf vs. Offshore")
t2 <- triangle.plot(hi.het2) + ggtitle("Coastal Atlantic  vs. Offshore")

ggsave(file="figures/triangle_plot_fourpop.png",
       ggpubr::ggarrange(t1, t2, common.legend = T),
       h=4, w=7)
ggsave(file="figures/triangle_plot_fourpop.pdf",
       ggpubr::ggarrange(t1, t2, common.legend = T),
       h=4, w=7)

# color by missing data:
missing.plot(hi.het1)


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# six pop:

pops <- data.frame(id = pops1$indiv,
                   pop = pops1$sixpop)
# identify ancestry informative markers (AIMS)
# Create a new vcfR object composed only of sites above the given allele frequency difference threshold
diff1 <- alleleFreqDiff(vcfR = vcf, 
                        pm = pops, 
                        p1 = "Coastal_Atlantic", 
                        p2 = "Offshore_Atlantic", 
                        difference = 0.7)
#[1] "696  sites passed allele frequency difference threshold"

diff2 <- alleleFreqDiff(vcfR = vcf, 
                        pm = pops, 
                        p1 = "Coastal_Gulf", 
                        p2 = "Offshore_Gulf", 
                        difference = 0.7)
#[1] "301 sites passed allele frequency difference threshold"

write.table(data.frame(snps=row.names(extract.gt(diff2))),
            quote = F,
            col.names = F,
            file="analysis/snp_set_gulfOnly.txt",
            sep="\t",
            row.names=F
)

diff2 <- alleleFreqDiff(vcfR = vcf, 
                        pm = pops, 
                        p1 = "Coastal_Gulf", 
                        p2 = "Offshore_Gulf", 
                        difference = 0.7)

# Calculate hybrid index and heterozygosity for each sample. 
# Values are returned in a data.frame
hi.het1_sixpop <- hybridIndex(vcfR = diff1, 
                       pm = pops, 
                       p1 = "Coastal_Atlantic", p2 = "Offshore_Atlantic")

hi.het2_sixpop <- hybridIndex(vcfR = diff2, 
                       pm = pops, 
                       p1 = "Coastal_Gulf", p2 = "Offshore_Gulf")

plot(hi.het1_sixpop$hybrid.index, hi.het2$hybrid.index)
plot(hi.het2_sixpop$heterozygosity, hi.het2$heterozygosity)
# they are largely consistent

# Generate colors (or leave blank to use default)

# View triangle plot
#triangle.plot(hi.het, colors = cols)
t1 <- triangle.plot(hi.het1_sixpop,colors = cols) + ggtitle("Coastal Atlantic vs. Offshore Atlantic")
t2 <- triangle.plot(hi.het2_sixpop,colors = cols) + ggtitle("Coastal Gulf  vs. Offshore Gulf")

ggsave(file="figures/triangle_plot_sixpop.png",
       ggpubr::ggarrange(t2, t1, common.legend = T),
       h=4, w=7)
ggsave(file="figures/triangle_plot_sixpop.pdf",
       ggpubr::ggarrange(t2, t1, common.legend = T),
       h=4, w=7)

# color by missing data:
missing.plot(hi.het1_sixpop)


#---------------------------------------------------------------------
#---------------------------------------------------------------------
# compare between the six pop and two pop

hi.het1_sixpop
hi.het2_sixpop

plot(hi.het1_sixpop$hybrid.index, hi.het1$hybrid.index)
plot(hi.het2_sixpop$hybrid.index, hi.het2$hybrid.index)

pops <- read.table("analysis/population_assignments_summary.txt", header=T)
hybs <- pops$indiv[pops$offshore_putative_hybrids == TRUE]

hi.het1_sixpop$pop[hi.het1_sixpop$id %in% hybs] <- "putative_hybrid"
hi.het2_sixpop$pop[hi.het2_sixpop$id %in% hybs] <- "putative_hybrid"

t1 <- triangle.plot(hi.het1_sixpop,colors = cols) + ggtitle("Coastal Atlantic vs. Offshore Atlantic")
t2 <- triangle.plot(hi.het2_sixpop,colors = cols) + ggtitle("Coastal Gulf  vs. Offshore Gulf")

ggsave(file="figures/triangle_plot_sixpop_hybrids.png",
       ggpubr::ggarrange(t2, t1, common.legend = T),
       h=4, w=7)
ggsave(file="figures/triangle_plot_hybrids.pdf",
       ggpubr::ggarrange(t2, t1, common.legend = T),
       h=4, w=7)

