### alignment stats
library(ggplot2)
dat <- read.csv("analysis/alignment_stats.csv")
dat$sample <- gsub("\\.merged\\.sorted", "", dat$sample)
dflib <- read.csv("merge_summary.csv", header=F)
colnames(dflib) <- c("sample", "lanes", "dir")
hist(dat$map_percent_q20)

all <- merge(dat, dflib, by="sample")

hist(dat$mapped_q20, breaks=30)

ggplot(all, aes(x = total_reads, fill = as.factor(lanes))) +
  geom_histogram(position = "identity", alpha = 0.9, bins = 40) +
  geom_vline(xintercept = quantile(all$total_mapped, c(0.025)), 
             linetype = "dashed", color = "black", size=2) +
  scale_fill_discrete(name = "Lanes") +
  theme_classic(base_size = 14) +
  labs(x = "Total Reads", 
       y = "Count",
       title = "Total # of Reads by Number of Lanes")

grep("b", all$sample)

# there are some samples originally included that are the wrong species. we knew this ahead of time. Drop them. 

## 41 are from brazil
all$sample[grep("^41",all$sample)]

all2 <- all[grep("^41",all$sample, invert=T),]
nrow(all2)
#385
# another species
all2$sample[grep("^157",all2$sample)]

all3 <- all2[grep("^157",all2$sample, invert=T),]
nrow(all3)
#380

# drop Tadu
all4 <- all3[grep("^Tadu",all3$sample, invert=T),]
nrow(all4)

all4$sample[grep("42193|78068",all3$sample, invert=F)]
#[1] "42193" "78068"
all5 <- all4[grep("42193|78068",all4$sample, invert=T),]

nrow(all5)
# 375 of the right species

all5[order(all5$total_reads),]

ggplot(all5, aes(x = total_reads, fill = as.factor(lanes))) +
  geom_histogram(position = "identity", alpha = 0.9, bins = 40) +
  geom_vline(xintercept = 1000000, 
             linetype = "dashed", color = "black", size=1.3) +
  scale_fill_discrete(name = "Lanes") +
  theme_classic(base_size = 12) +
  labs(x = "Total Reads", 
       y = "Count",
       title = "Total # of Reads by Number of Lanes")

table(all5$lanes)
ggsave("figures/bam_reads_hist.png", h=3,w=4)

length(all5$sample[which(all5$total_reads < 1000000)])
all5$sample[which(all5$total_reads < 1000000)]
length(all5$sample[which(all5$total_reads > 1000000)])
# 372

all5[which(all5$sample == "26Tt309"),]
all5[which(all5$sample == "26Tt313"),]
all5[which(all5$sample == "2Tt156"),]
all5[which(all5$sample == "7Tt278"),]

write.table(file="scripts/bam.list",
            paste0("/home/rbrennan/Tursiops-NC-PopulationAssignment-RAD/analysis/merged_bams/",
    all5$sample[which(all5$total_reads > 1000000)],".merged.sorted.bam"), row.names=F, quote=F, col.names =F)





















