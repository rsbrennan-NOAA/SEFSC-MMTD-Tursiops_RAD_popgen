# new hybrids

# https://github.com/bwringe/hybriddetective

# https://groups.google.com/g/dartr/c/dvpgv6C6Ras

library(tidyverse)
pofz_path <- "analysis/pop_structure/newhybrids/"
inds_file <- "analysis/variants/indivs_gulfVSoffshore.txt"

pofz_files <- list.files(path = pofz_path, pattern = "aa-PofZ_\\d+\\.txt", full.names = TRUE)
indiv_names <- read_delim(inds_file, delim = '\t', col_names = FALSE) %>% pull(X1)
length(indiv_names)

allres <- pofz_files %>%
  # getrun number from each file path
  set_names(str_extract(., "(?<=aa-PofZ_)\\d+(?=\\.txt)")) %>%
  # Read all files and combine 
  map_dfr(~read_delim(.x, 
                      delim = '\t', 
                      skip = 1,
                      col_names = c("indNR", "indiv", "Pure1", "Pure2", "F1", "F2", "Backcross1", "Backcross2")),
          .id = "run_number") %>%
  # Convert run_number to num
  mutate(run_number = as.numeric(run_number)) %>%
  # Group by individual number
  group_by(indNR) %>%
  # Assign the id to each
  mutate(indiv = indiv_names[as.numeric(indNR)]) %>%
  ungroup()


pops <- read_delim("analysis/population_assignments_summary.txt")

allres <- allres %>%
  left_join(pops, by = "indiv")


library(ggbeeswarm)
# For each individual in each run, find the highest probability category
allres_with_highest <- allres %>%
  mutate(
    highest_value = pmax(Pure1, Pure2, F1, F2, Backcross1, Backcross2),
    highest_category = case_when(
      highest_value == Pure1 ~ "Pure1",
      highest_value == Pure2 ~ "Pure2",
      highest_value == F1 ~ "F1",
      highest_value == F2 ~ "F2",
      highest_value == Backcross1 ~ "Backcross1",
      highest_value == Backcross2 ~ "Backcross2"
    )
  )

head(allres_with_highest)

# Analyze consistency across runs for each individual
consistency_analysis <- allres_with_highest %>%
  group_by(indiv, sixpop) %>%
  summarize(
    num_runs = length(highest_category),
    # Count occurrences of each category
    Pure1_count = sum(highest_category == "Pure1"),
    Pure2_count = sum(highest_category == "Pure2"),
    F1_count = sum(highest_category == "F1"),
    F2_count = sum(highest_category == "F2"),
    Backcross1_count = sum(highest_category == "Backcross1"),
    Backcross2_count = sum(highest_category == "Backcross2"),
    # Find most common category
    most_common_category = names(which.max(c(
      Pure1 = Pure1_count, 
      Pure2 = Pure2_count,
      F1 = F1_count,
      F2 = F2_count,
      Backcross1 = Backcross1_count,
      Backcross2 = Backcross2_count
    ))),
    # Calculate consistency percentage
    consistency_percent = 100 * max(Pure1_count, Pure2_count, F1_count, F2_count, Backcross1_count, Backcross2_count) / num_runs,
    # Calculate mean of highest values
    mean_highest_value = mean(highest_value),
    .groups = "drop"
  ) %>%
  mutate(is_consistent = consistency_percent == 100)
      
sum(consistency_analysis$is_consistent)


#summary stats
consistency_summary <- consistency_analysis %>%
  group_by(sixpop) %>%
  summarize(
    total_individuals = n(),
    fully_consistent_individuals = sum(is_consistent),
    percent_consistent = 100 * fully_consistent_individuals / total_individuals,
    .groups = "drop"
  )
      



# Calculate median assignment probabilities for each individual
median_assignments <- allres %>%
  group_by(indiv) %>%
  summarize(
    Pure1_median = median(Pure1),
    Pure2_median = median(Pure2),
    F1_median = median(F1),
    F2_median = median(F2),
    Backcross1_median = median(Backcross1),
    Backcross2_median = median(Backcross2),
    # Keep the dapc_population information (assuming it's the same for each indiv)
    sixpop = first(sixpop),
    .groups = "drop"
  ) %>%
  # Find the highest median category
  mutate(
    highest_value = pmax(Pure1_median, Pure2_median, F1_median, F2_median, Backcross1_median, Backcross2_median),
    highest_category = case_when(
      highest_value == Pure1_median ~ "Pure1",
      highest_value == Pure2_median ~ "Pure2",
      highest_value == F1_median ~ "F1",
      highest_value == F2_median ~ "F2",
      highest_value == Backcross1_median ~ "Backcross1",
      highest_value == Backcross2_median ~ "Backcross2"
    )
  )

median_assignments$sig <- "Low_Confidence"
median_assignments$sig[median_assignments$highest_value > 0.95] <- "High_Confidence"

category_counts <- median_assignments %>%
  group_by(sixpop, highest_category, sig) %>%
  summarize(count = n(), .groups = "drop")

p <- ggplot(median_assignments, aes(x = sixpop, y = highest_value, 
                               fill = highest_category,
                               color = highest_category)) +
  geom_boxplot(alpha = 0.7) +
  # Add count labels (just the number)
  geom_text(data = category_counts, 
            aes(label = count, y = 0.05, group = highest_category),
            position = position_dodge(width = 0.75),
            size = 3.5, fontface = "bold") +
  scale_fill_brewer(palette = "Set1", name = "Hybrid Category") +
  scale_color_brewer(palette = "Set1", guide = "none") +
  labs(
    title = "NewHybrids Posteriors",
    x = "Population",
    y = "Median posterior of 10 runs"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  facet_wrap(vars(sig))

ggsave(file="figures/newhybrids.png", p, 
       h=4, w=7)

# write output:
write.table(median_assignments, 
            file="analysis/pop_structure/newhybrids/newhybrids_posteriors.txt",
            quote = F,
            row.names = F, sep="\t")

filtered_dat <- median_assignments %>%
  filter(highest_category != "Pure1" & highest_category != "Pure2")



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
diff1 <- alleleFreqDiff(vcfR = vcf, 
                        pm = pops, 
                        p1 = "Coastal_Gulf", 
                        p2 = "Offshore", 
                        difference = 0.7)


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



pops <- read.table("analysis/population_assignments_summary.txt", header=T)
hybs <- pops$indiv[pops$offshore_putative_hybrids == TRUE]

hi.het1$pop[hi.het1$id %in% hybs] <- "putative_hybrid"
hi.het2$pop[hi.het2$id %in% hybs] <- "putative_hybrid"

t1 <- triangle.plot(hi.het1,colors = cols) + ggtitle("Coastal Atlantic vs. Offshore Atlantic")
t2 <- triangle.plot(hi.het2,colors = cols) + ggtitle("Coastal Gulf  vs. Offshore Gulf")


pops <- read.table("analysis/population_assignments_summary.txt", header=T)

pops2 <- read.table("analysis/pop_structure/newhybrids/newhybrids_posteriors.txt", header=T)

hybs <- pops$indiv[pops$offshore_putative_hybrids == TRUE]
newhybrids <- pops2[pops2$sig == "High_Confidence" & pops2$highest_category != "Pure1" & pops2$highest_category != "Pure2" ,]

F2 <- newhybrids$indiv[newhybrids$highest_category == "F2"]
Backcross2  <- newhybrids$indiv[newhybrids$highest_category == "Backcross2"]
Backcross1 <- newhybrids$indiv[newhybrids$highest_category == "Backcross1"]

hi.het1$pop[hi.het1$id %in% hybs] <- "putative_hybrid"
hi.het2$pop[hi.het2$id %in% hybs] <- "putative_hybrid"

hi.het1$pop[hi.het1$id %in% F2] <- "NewHybrids_F2"
hi.het2$pop[hi.het2$id %in% F2] <- "NewHybrids_F2"

hi.het1$pop[hi.het1$id %in% Backcross2] <- "NewHybrids_Backcross2"
hi.het2$pop[hi.het2$id %in% Backcross2] <- "NewHybrids_Backcross2"

hi.het1$pop[hi.het1$id %in% Backcross1] <- "NewHybrids_Backcross1"
hi.het2$pop[hi.het2$id %in% Backcross1] <- "NewHybrids_Backcross1"

library(paletteer)
cols <- c("#E41A1CFF", "#377EB8FF", "#4DAF4AFF", "#984EA3FF", 
          "#FF7F00FF", "#EEC229FF", "#A65628FF", "#F781BFFF", 
          "#999999FF", "black")


t1 <- triangle.plot(hi.het1,colors = cols) + ggtitle("Coastal Atlantic vs. Offshore Atlantic")
t2 <- triangle.plot(hi.het2,colors = cols) + ggtitle("Coastal Gulf  vs. Offshore Gulf")

ggsave(file="figures/triangle_plot_sixpop_Newhybrids.png",
       ggpubr::ggarrange(t2, t1, common.legend = T),
       h=4, w=7.5)
ggsave(file="figures/triangle_plot_Newhybrids.pdf",
       ggpubr::ggarrange(t2, t1, common.legend = T),
       h=4, w=7.5)




#--------------------------------------------------------------
# add the  hybrids to the map:
#--------------------------------------------------------------

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(ggspatial)
library(marmap)
library(RColorBrewer)



# Get the world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# read in locations
dat <- read.csv("Tursiops_RADseq_Metadata_new.csv")

head(dat)

usa <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

pops <- read.table("analysis/population_assignments_summary.txt", header=T)
pops2 <- read.table("analysis/pop_structure/newhybrids/newhybrids_posteriors.txt", header=T)

dat <- read.csv("Tursiops_RADseq_Metadata_new.csv")

out <- merge(dat, pops, by.x="Lab.ID", by.y="indiv")
nhybs <- merge(dat, pops2, by.x="Lab.ID", by.y="indiv")


F2 <- nhybs %>% filter(sig == "High_Confidence" & highest_category == "F2")
Backcross1 <- nhybs %>% filter(sig == "High_Confidence" & highest_category == "Backcross1")
Backcross2 <- nhybs %>% filter(sig == "High_Confidence" & highest_category == "Backcross2")

p <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_point(data = out, 
             aes(x = Long, y = Lat, fill=fourpop, shape=fourpop),
             size = 2,
             alpha=1,
             color="black") +
  #geom_point(data = out_hyb, 
  #           aes(x = Long, y = Lat),
  #           fill="red", 
  #           shape=21,
  #           size = 4,
  #           alpha=1,
  #           color="black") +
  geom_point(data = Backcross1, 
             aes(x = Long, y = Lat),
             fill="red", 
             shape=21,
             size = 4,
             alpha=1,
             color="black") +
  geom_point(data = Backcross2, 
             aes(x = Long, y = Lat),
             fill="red", 
             shape=21,
             size = 4,
             alpha=1,
             color="black") +
  geom_point(data = F2, 
             aes(x = Long, y = Lat),
             fill="red", 
             shape=21,
             size = 4,
             alpha=1,
             color="black") +
  geom_text(data = F2,
            aes(x = Long, y = Lat, label = "F2"),
            size = 2,
            color = "black") +
  geom_text(data = Backcross1,
            aes(x = Long, y = Lat, label = "B1"),
            size = 2,
            color = "black") +
  geom_text(data = Backcross2,
            aes(x = Long, y = Lat, label = "B2"),
            size = 2,
            color = "black") +
  coord_sf() +
  theme_bw() +
  theme(
    #panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    #axis.text = element_blank(),
    #axis.ticks = element_blank()
    legend.position = "top",
    legend.title=element_blank()) +
  xlab("Longitude")+
  ylab("Latitude") +
  scale_shape_manual(values=c(21,22,23, 24)) +
  coord_sf(xlim = c(-100, -60), ylim = c(23, 47), expand = FALSE)+
  annotation_scale()

p

ggsave("figures/map_allpops_hybrids.pdf", p, h=4, w=5)
ggsave("figures/map_allpops_hybrids.png", p, h=4, w=5)


# save output:

pops <- read.table("analysis/population_assignments_summary.txt", header=T)
pops2 <- read.table("analysis/pop_structure/newhybrids/newhybrids_posteriors.txt", header=T)

dat <- read.csv("Tursiops_RADseq_Metadata_new.csv")

out <- merge(dat, pops, by.x="Lab.ID", by.y="indiv")
nhybs <- merge(dat, pops2, by.x="Lab.ID", by.y="indiv")


F2 <- nhybs %>% filter(sig == "High_Confidence" & highest_category == "F2")
Backcross1 <- nhybs %>% filter(sig == "High_Confidence" & highest_category == "Backcross1")
Backcross2 <- nhybs %>% filter(sig == "High_Confidence" & highest_category == "Backcross2")


pops$newhybrids_category <- NA
pops$newhybrids_category[pops$indiv %in% F2$Lab.ID] <- "F2"
pops$newhybrids_category[pops$indiv %in% Backcross1$Lab.ID] <- "Backcross1"
pops$newhybrids_category[pops$indiv %in% Backcross2$Lab.ID] <- "Backcross2"
table(pops$newhybrids_category)


# add in the low confidence offshore indivs that all cluster together. remember to give the rationale for this in the MS

offhybs <- pops$indiv[which(pops$offshore_putative_hybrids == TRUE)]
length(offhybs)
# 16
#allhybs <- nhybs[nhybs$highest_category != "Pure1" & nhybs$highest_category != "Pure2",]
allhybs <- nhybs[nhybs$highest_category == "Backcross1" | nhybs$highest_category == "Backcross2" | nhybs$highest_category == "F2",]
nrow(allhybs)
#26
highconf <- allhybs[allhybs$sig == "High_Confidence",]
nrow(highconf)
# 20
allhybs %>% group_by(sixpop, highest_category, sig) %>% summarise(count = n())


#low_conf_offshore <- 

nrow(highconf[highconf$Lab.ID %in% offhybs,])
# 12 of the 16
low_conf_offshore <- allhybs[allhybs$Lab.ID %in% offhybs & allhybs$sig == "Low_Confidence",]
nrow(low_conf_offshore)
# the other 4 are hybs, but considered low confidence. add these in
allNewHybs <- (c(highconf$Lab.ID, low_conf_offshore$Lab.ID))
length(allNewHybs)
#24

for(i in 1:nrow(low_conf_offshore)){
  pops$newhybrids_category[pops$indiv == low_conf_offshore$Lab.ID[i]] <- low_conf_offshore$highest_category[i]
}

sum(table(pops$newhybrids_category))
#24

write.table(pops,"analysis/population_assignments_hybrids_summary.txt", sep="\t",
            row.names = FALSE, quote = FALSE)

# for the ones that look like hybrids, but new hybrids doesn't assign, what is going on?
# they're mostly just low confidence. 

####### write a clust file without the hybrids

onlyhybs <- pops[!is.na(pops$newhybrids_category),]

write.table(onlyhybs$indiv,
            file="analysis/pop_structure/newhybrids/hybrids.txt", 
            sep="\t", quote=F, col.names=F, row.names=F)




# write hybrids to drop from vcf:
write.table(filtered_dat$indiv, 
            file="analysis/pop_structure/newhybrids/hybrids.txt",
            quote = F,
            row.names = F, sep="\t", col.names=F)




       
