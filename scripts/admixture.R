# admixture
library(tidyverse)
library(dplyr)
library(forcats)
library(patchwork)

### what is our most likely k?

# read in CV scores:
cvin <- read.csv("analysis/pop_structure/cv.txt", sep=":", header=F)
colnames(cvin) <- c("id", "cv")
# fix the formatting to get K into numeric format
cvin$K <- as.numeric(str_extract(cvin$id, "\\d+"))
cvin[which.min(cvin$cv),]

# plot the results
p <- ggplot(cvin,aes(x=K,y=cv)) +
  geom_point(size=3)  + geom_line(group=1)

p
ggsave("figures/cv.png", p, h=3, w=3)

# actual results:

samplelist <- read_delim("analysis/pop_structure/LDthin_numCorrect.fam",
                         col_names = c("individual", "id2", "a", "b", "c", "d"),
                         delim=" ")

# read in all date, in a loop
## first create an empty dataframe
all_data <- tibble(individual=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

# then add all results to this
for (k in 2:8){
  data <- read_delim(paste0("analysis/pop_structure/LDthin_numCorrect.",k,".Q"),
                     col_names = paste0("Q",seq(1:k)),
                     delim=" ")
  data$sample <- samplelist$individual
  data$k <- k
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k) -> data
  all_data <- rbind(all_data,data)
}

pops <- read.csv("Tursiops_RADseq_Metadata_new.csv")
pops$region <- ifelse(pops$Long > -81.7, "Atlantic",
                      ifelse(pops$Long > 80, "Gulf", NA))
pops$region <- ifelse(pops$Long < -81.7, "Gulf", pops$region)

depth <- read.csv("depths.csv", header=T)

df <- merge(pops, depth, by.x="Lab.ID", by.y="id")

#all_data$IDs <- gsub("b", "",all_data$sample)

result <- all_data %>%
  left_join(df %>% select( Lab.ID, region, corrected_depth), by = c("sample" = "Lab.ID"))

all_data <- result
#pop_sub <- pops[pops$Lab.ID.. %in% result$IDs,]

# order samples by population:
sampleorder <- df$Lab.ID[order(df$region, df$corrected_depth)]
all_data$IDs <- as.factor(all_data$sample)
all_data$IDs <- factor(all_data$IDs, levels=sampleorder)
all_data$sample <- all_data$IDs

all_data %>%
  filter(k == 2) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_rug(aes(x=sample, y=value, color=region)) +
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  scale_color_manual(values=c("#eac435","#557fc3", "#03cea4", "#fb4d3d"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",                    
                    labels=c("1","2"))

# within each population, order indivs by q val
new_dat<- all_data[all_data$k == 2 & all_data$Q == "Q1",]

new_dat2 <- new_dat %>%
  mutate(sample = fct_reorder2(sample, region, value)) %>%
  arrange(region, value)

all_data$sample <- factor(all_data$sample, levels=c(new_dat2$sample))
all_data$k <- as.numeric(all_data$k)

p2 <-  
  all_data %>%
  filter(k < 8) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  geom_rug(aes(x=sample, y=value, color=region),
           linewidth = 1,
           length= unit(0.06, "npc"),
           sides="b",
           outside = TRUE) +
  coord_cartesian(clip = "off")+
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_color_manual(values = c("black","grey55"),
                     name = "Population") +
  scale_fill_brewer(palette = "Set1", name = "K") +
  facet_wrap(~k,ncol=1)
p2


combined_plot <- wrap_plots(p, p2, heights = c(0.15, 1), ncol=1)
combined_plot

ggsave("figures/Admixture_plot.pdf", combined_plot, width = 7, height = 8, units="in")
ggsave("figures/Admixture_plot.png", combined_plot, width = 7, height = 8, units="in")

#-----------------------------------------------------
#-----------------------------------------------------
# how to best summarize these populations?
# most likely k is 4, lets try to assign indivs based on this.
k <- 4
data <- read_delim(paste0("analysis/pop_structure/LDthin_numCorrect.",k,".Q"),
                   col_names = paste0("Q",seq(1:k)),
                   delim=" ")
data$sample <- samplelist$individual
head(data)
idxq1 <- which(data$Q1> 0.7)
idxq2 <- which(data$Q2> 0.7)
idxq3 <- which(data$Q3> 0.5)
idxq4 <- which(data$Q4> 0.7)
length(idxq1)
length(idxq2)
length(idxq3)
length(idxq4)

sub_id <- data$sample[sample(idxq3, size=60, replace=F)]
rest_id <- data$sample[which(data$Q2 < 0.7)]

write.table(file= "subset.txt", data.frame(V1= c(sub_id, rest_id), V2= c(sub_id, rest_id)),
            col.names=F, row.names=F, quote=F, sep="\t")



#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
# analysis of subset data, to make sure consistent

# read in CV scores:
cvin <- read.csv("analysis/pop_structure/cv_subset.txt", sep=":", header=F)
colnames(cvin) <- c("id", "cv")
# fix the formatting to get K into numeric format
cvin$K <- as.numeric(str_extract(cvin$id, "\\d+"))
cvin[which.min(cvin$cv),]

# plot the results
p <- ggplot(cvin,aes(x=K,y=cv)) +
  geom_point(size=3)  + geom_line(group=1)

p
ggsave("figures/cv_subset.png", p, h=3, w=3)

# actual results:

samplelist <- read_delim("analysis/pop_structure/LDthin_numCorrect_subset.fam",
                         col_names = c("individual", "id2", "a", "b", "c", "d"),
                         delim=" ")

# read in all date, in a loop
## first create an empty dataframe
all_data <- tibble(individual=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

# then add all results to this
for (k in 2:8){
  data <- read_delim(paste0("analysis/pop_structure/LDthin_numCorrect_subset.",k,".Q"),
                     col_names = paste0("Q",seq(1:k)),
                     delim=" ")
  data$sample <- samplelist$individual
  data$k <- k
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k) -> data
  all_data <- rbind(all_data,data)
}

pops <- read.csv("Tursiops_RADseq_Metadata_new.csv")
pops$region <- ifelse(pops$Long > -81.7, "Atlantic",
                      ifelse(pops$Long > 80, "Gulf", NA))
pops$region <- ifelse(pops$Long < -81.7, "Gulf", pops$region)

depth <- read.csv("depths.csv", header=T)

df <- merge(pops, depth, by.x="Lab.ID", by.y="id")

#all_data$IDs <- gsub("b", "",all_data$sample)

result <- all_data %>%
  left_join(df %>% select( Lab.ID, region, corrected_depth), by = c("sample" = "Lab.ID"))

all_data <- result
#pop_sub <- pops[pops$Lab.ID.. %in% result$IDs,]

# order samples by population:
sampleorder <- df$Lab.ID[order(df$region, df$corrected_depth)]

df$corrected_depth[order(df$corrected_depth)]
all_data$IDs <- as.factor(all_data$sample)
all_data$IDs <- factor(all_data$IDs, levels=sampleorder)
all_data$sample <- all_data$IDs

all_data %>%
  filter(k == 2) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_rug(aes(x=sample, y=value, color=region)) +
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  scale_color_manual(values=c("#eac435","#557fc3", "#03cea4", "#fb4d3d"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",                    
                    labels=c("1","2"))

# within each population, order indivs by q val
#new_dat<- all_data[all_data$k == 2 & all_data$Q == "Q1",]

#new_dat2 <- new_dat %>%
 # mutate(sample = fct_reorder2(sample, region, value)) %>%
#  arrange(region, value)

#all_data$sample <- factor(all_data$sample, levels=c(new_dat2$sample))
#all_data$k <- as.numeric(all_data$k)

p2 <-  
  all_data %>%
  filter(k < 8) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  geom_rug(aes(x=sample, y=value, color=region),
           linewidth = 1,
           length= unit(0.06, "npc"),
           sides="b",
           outside = TRUE) +
  coord_cartesian(clip = "off")+
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_color_manual(values = c("black","grey55"),
                     name = "Population") +
  scale_fill_brewer(palette = "Set1", name = "K") +
  facet_wrap(~k,ncol=1)
p2


combined_plot <- wrap_plots(p, p2, heights = c(0.15, 1), ncol=1)
combined_plot

ggsave("figures/Admixture_plot_subset.pdf", combined_plot, width = 7, height = 8, units="in")
ggsave("figures/Admixture_plot_subset.png", combined_plot, width = 7, height = 8, units="in")


