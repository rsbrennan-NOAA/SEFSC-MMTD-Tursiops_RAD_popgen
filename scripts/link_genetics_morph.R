#morphology overlap with genetics

library(ggplot2)
library(dplyr)
library(tidyverse)

# use list from the AllATL_matadata. these are the micros we need
# then make sure the morphology ones from ana are in there.
# use this table to get the base microsat list: https://docs.google.com/spreadsheets/d/1M-0P6VRTiMthxtVIjww3YzKK-o4VTK6nbPbpP3MiMg0/edit?gid=60770286#gid=60770286

# we think there are 47 radseq ones missing from the micro data. 
  # check that my number agrees

# do we really need to re-run what we've already done.

# DO the microsat assignments agree with the RAD asisngments with k=4
# add haplotype infor to the 

# end should be a big table with the microsat loci, if RAD, if morphology

skull_ids <-read.csv("analysis/skull_labID.csv", header=T)
nrow(skull_ids)
# 182
sum(!is.na(skull_ids$LabID))
# 105

morph <- read.csv("analysis/morphology_list.txt", sep=" ", header=T)
head(morph)

#4Tt328 is 4Tt184 so substitute them in the morph list:

morph$lab_id[morph$lab_id == "4Tt328"] <- "4Tt184"

# add the new samples
new_morph <- read.csv("analysis/ListofSamplesMeasured-Mote.csv")

all_morph_ids <- c(morph$lab_id, new_morph$Lab.ID)

dat <- read.csv("metadata_FINAL.csv")

dat$id <- gsub("b$", "", dat$id)

dat[(dat$id %in% all_morph_ids),]

# read in pca coords:
pccoord <- read.csv("analysis/TUR-ATL-PCAcoordinates.txt", sep="\t")
nrow(pccoord)
#180

skull_ids$MuseumID %in% pccoord$sample_name
pccoord2 <- merge(skull_ids, pccoord, by.x="MuseumID", by.y="sample_name")
pccoord <- pccoord2
nrow(pccoord)
#180
sum(!is.na(pccoord$LabID))
#105

dat[(dat$id %in% pccoord$LabID),]


# I need to check the allsamp, because sometimes there are different lab ids for genotyping than for morph
# but these are the same samples. 

# look through allsamp, to find the matching rows with pccoord$MuseumID:
# add a column for the lab id match. do they agree? 

# add lab ids to the pc coords:
allsamp <- read.csv("analysis/DatabaseDump_2025July30.txt", sep="\t")

library(tidyr)

allsamp$Lab_ID <- gsub("-","",allsamp$Lab_ID)
allsamp$Field_ID <- gsub("-","",allsamp$Field_ID)
allsamp$Lab_ID <- gsub("\n", "", allsamp$Lab_ID)

matching_rows_all <- list()
matching_rows_check <- c()
pc_match_id <- c()
pc_match_labid <- c()

pccoord <- pccoord %>% rename(sample_name = MuseumID)

# find if any pca indivs are in the database... most should be

#make new df to save output.
# new dataframe with sample_name and LabID from pccoord
new_pccoord <- data.frame(
  morph_Field_ID = pccoord$sample_name,
  morph_Lab_ID = pccoord$LabID
)

# empty dataframe to store matched allsamp rows
matched_allsamp <- allsamp[0, ]

for(pcID in 1:length(pccoord$sample_name)) {
  matching_rows_tmp <- apply(allsamp, 1, function(x) any(grepl(pccoord$sample_name[pcID], x,
                                                               ignore.case = TRUE))) 
  
  if(sum(matching_rows_tmp) > 1){
    print(cat(pccoord$sample_name[pcID],"more than 1 match"))
    if(pccoord$sample_name[pcID] == "SC0039"){
      print("dropping SC0039, keeping 4Tt095")
      matching_rows_tmp <- allsamp$Field_ID == "FB811"
    }
    if(pccoord$sample_name[pcID] == "SC9707"){
      print("dropping CMSC9707")
      matching_rows_tmp <- allsamp$Field_ID == "SC9707"
    }
    if(pccoord$sample_name[pcID] == "SC0704"){
      print("dropping SC0704F ")
      matching_rows_tmp <- allsamp$Field_ID == "SC0704"
    }
  }
  
  # add  matching row or NA to matched_allsamp
  if(sum(matching_rows_tmp) > 0) {
    matched_allsamp <- rbind(matched_allsamp, allsamp[which(matching_rows_tmp)[1], ])
  } else {
    # make NA row
    na_row <- allsamp[1, ]
    na_row[1, ] <- NA
    matched_allsamp <- rbind(matched_allsamp, na_row)
  }
}


matched_df <- cbind(new_pccoord, matched_allsamp)

nrow(matched_df)
# 180

sum(!is.na(matched_df$Lab_ID))
sum(!is.na(matched_df$morph_Lab_ID ))

matched_df[which(matched_df$morph_Lab_ID != matched_df$Lab_ID),]

na_in_one_not_both <- (is.na(matched_df$morph_Lab_ID) & !is.na(matched_df$Lab_ID)) | 
  (!is.na(matched_df$morph_Lab_ID) & is.na(matched_df$Lab_ID))
# some aren't being match correctly
matched_df[na_in_one_not_both, ]

sum(na_in_one_not_both)


# fix the ones where the morph field IDs are wrong so we're not getting lab ids.

for(i in which(na_in_one_not_both)) {
  morph_lab_id <- matched_df$morph_Lab_ID[i]
  
  # search allsamp for this morph_Lab_ID
  matching_rows_tmp <- apply(allsamp, 1, function(x) any(grepl(morph_lab_id, x, ignore.case = TRUE)))
  
  # if we find a match, replace the NA row with the matching allsamp row
  if(sum(matching_rows_tmp) > 0) {
    # Replace the allsamp columns (everything except morph_Field_ID and morph_Lab_ID)
    matched_df[i, 3:ncol(matched_df)] <- allsamp[which(matching_rows_tmp)[1], ]
  }
}

# Check the corrected rows
matched_df[na_in_one_not_both, ]

nrow(matched_df)
# 180 

sum(!(is.na(matched_df$Lab_ID)))
# 105

#--------------------------------------------------
# link to RAD id:
#--------------------------------------------------

# from matched_df, figure out which have RAD genotypes:
matched_df$RAD_ID <- NA
matched_df$RAD_fourpop <- NA
matched_df$newhybrids_category <- NA

# from out_match, figure out which have RAD genotypes:
dat <- read.csv("metadata_FINAL.csv")

dat$id <- gsub("b$", "", dat$id)

rad_match <- c()
for(i in dat$id){
  tmp_match <- grep(i, matched_df$Lab_ID)
  if(length(tmp_match) > 0){
    rad_match <- c(rad_match, i)
  }
}
rad_match
# only 5 match

# For each rad_match ID, find matching rows and add the data
for(i in rad_match) {
  tmp_match <- grep(i, matched_df$Lab_ID)
  if(length(tmp_match) > 0) {
    matched_df$RAD_ID[tmp_match] <- i
    # Find the corresponding row in dat
    dat_row <- which(dat$id == i)
    if(length(dat_row) > 0) {
      matched_df$RAD_fourpop[tmp_match] <- dat$fourpop[dat_row]
      matched_df$newhybrids_category[tmp_match] <- dat$newhybrids_category[dat_row]
    }
  }
}

matched_df[!is.na(matched_df$RAD_ID), ]
# 2 coastal gulf, 3 coastal atlantic


# add the micro assignments from nikki:
# basically, figure out generally if the RAD and micro assignments agree

micro_data <- read.csv("analysis/RADseqMicroData_n343_2025_NEW-RAD&MicroStructureData_n343.csv")
# 343 indivs

extract_ids <- function(lab_id_string) {
  # split by "="
  ids <- trimws(unlist(strsplit(as.character(lab_id_string), "=")))
  return(ids)
}

# make expanded dataframe with all possible ID matches
micro_expanded <- data.frame()
micro_data <- cbind(micro_data$Lab.ID, micro_data)
colnames(micro_data) <- c("ID.Original",colnames(micro_data)[2:ncol(micro_data)])
for(i in 1:nrow(micro_data)) {
  possible_ids <- extract_ids(micro_data$Lab.ID[i])
  
  # make a row for each possible ID
  for(id in possible_ids) {
    temp_row <- micro_data[i, ]
    temp_row$Lab.ID <- id
    micro_expanded <- rbind(micro_expanded, temp_row)
  }
}

micro_matched_RAD <- merge(micro_expanded, dat, by.x="Lab.ID", by.y="id")
nrow(micro_matched_RAD)
table(micro_matched_RAD$fourpop.y, 
      micro_matched_RAD$ClumppK4.0.50.cutoff, 
      useNA = "ifany")

#                     1   2   3   4 unassigned
#Coastal_Atlantic   0   0 164   1          1
#Coastal_Gulf       0  38   1   7          0
#Intermediate       0   0   0  59          1
#Offshore          65   0   0   0          6

# what is the accuracy:
343-2-8-1-6
326/346



micro_matched_RAD$combined_category <- paste(micro_matched_RAD$fourpop.y, 
                                             micro_matched_RAD$newhybrids_category, 
                                             sep = "_")
micro_matched_RAD$combined_category <- gsub("_NA", "",micro_matched_RAD$combined_category)

table(micro_matched_RAD$combined_category, 
                              micro_matched_RAD$ClumppK4.0.50.cutoff, 
                              useNA = "ifany")

#                            1   2   3   4 unassigned
# Coastal_Atlantic          0   0 164   1          1
# Coastal_Gulf              0  38   1   7          0
# Intermediate              0   0   0  51          1
# Intermediate_Backcross2   0   0   0   4          0
# Intermediate_F2           0   0   0   4          0
# Offshore                 55   0   0   0          0
# Offshore_Backcross1       3   0   0   0          3
# Offshore_F2               7   0   0   0          3



#---------------------------------------------------------
# now read in full microsat dataset:
#---------------------------------------------------------
allmicro <- read.csv("analysis/RADseqMicroData_n343_2025_NEW-MicroStructureData_n4189.csv")
nrow(allmicro)
# 4189

table(allmicro$Source)
sum(allmicro$Source == "stranding")

allmicro <- cbind(allmicro$Lab.ID, allmicro)
colnames(allmicro) <- c("ID.Original",colnames(allmicro)[2:ncol(allmicro)])

#  empty dataframe to store expanded rows
allmicro_expanded <- allmicro[0, ]

for(i in 1:nrow(allmicro)) {
  current_row <- allmicro[i, ]
  
  # Split Lab.ID column
  lab_id_value <- as.character(current_row$Lab.ID)
  
  # NA or empty values
  if(is.na(lab_id_value) || lab_id_value == "") {
    splits <- NA
  } else {
    splits <- trimws(unlist(strsplit(lab_id_value, "=")))
    # Handle case where split results in empty character vector
    if(length(splits) == 0) {
      splits <- NA
    }
  }
  
  # Create multiple rows for each Lab.ID split
  for(j in 1:length(splits)) {
    new_row <- current_row
    
    # Replace Lab.ID with the split value
    if(!is.na(splits[j])) {
      new_row$Lab.ID <- splits[j]
    } else {
      new_row$Lab.ID <- NA
    }
    
    allmicro_expanded <- rbind(allmicro_expanded, new_row)
  }
}
head(allmicro_expanded)
nrow(allmicro_expanded)
# 4530 indivs

#----------------------------------------------------
#----------------------------------------------------
# merge with morphology data
matched_df <-cbind(matched_df$Lab_ID, matched_df)
colnames(matched_df) <- c("ID.Original",colnames(matched_df)[2:ncol(matched_df)])

matched_df_expanded <- matched_df[0, ]
# loop through each row
for(i in 1:nrow(matched_df)) {
  current_row <- matched_df[i, ]
  
  # get Lab_ID
  lab_id_value <- as.character(current_row$Lab_ID)
  
  # NA or empty values
  if(is.na(lab_id_value) || lab_id_value == "") {
    matched_df_expanded <- rbind(matched_df_expanded, current_row)
    next
  }
  
  # check if Lab_ID contains parentheses and split if so
  if(grepl("\\(", lab_id_value)) {
    # get the part before parentheses and the part inside parentheses
    before_paren <- trimws(sub("\\(.*\\)", "", lab_id_value))
    inside_paren <- trimws(gsub(".*\\((.*)\\).*", "\\1", lab_id_value))
    
    # make first row with value before parentheses
    new_row1 <- current_row
    new_row1$Lab_ID <- before_paren
    matched_df_expanded <- rbind(matched_df_expanded, new_row1)
    
    # make second row with value inside parentheses
    new_row2 <- current_row
    new_row2$Lab_ID <- inside_paren
    matched_df_expanded <- rbind(matched_df_expanded, new_row2)
    
  } else if(grepl("=", lab_id_value)) {
    # has "=" but no parentheses
    splits <- trimws(unlist(strsplit(lab_id_value, "=")))
    
    for(split_val in splits) {
      new_row <- current_row
      new_row$Lab_ID <- split_val
      matched_df_expanded <- rbind(matched_df_expanded, new_row)
    }
    
  } else {
    # no parentheses or "=" just add the row as is
    matched_df_expanded <- rbind(matched_df_expanded, current_row)
  }
}


head(matched_df_expanded)
nrow(matched_df_expanded)
# 188

merged_micro_morph <- merge(allmicro_expanded, matched_df_expanded, 
                   by.x = "Lab.ID", by.y = "Lab_ID", 
                   all = TRUE)

merged_micro_morph[!is.na(merged_micro_morph$ClumppK4.0.50.cutoff),]

merged_micro_morph_PCA <- merge(merged_micro_morph, pccoord, by.x="morph_Field_ID", by.y="sample_name")

# remove the dups from the previous step:
merged_micro_morph_PCA_unique <- merged_micro_morph_PCA[!duplicated(merged_micro_morph_PCA$morph_Field_ID), ]

nrow(merged_micro_morph_PCA_unique)
# 180 indivs
# how many new and how many added?
tmp_overlap <- merged_micro_morph_PCA_unique$Lab.ID %in% new_morph$Lab.ID
merged_micro_morph_PCA_unique[tmp_overlap,]
merged_micro_morph_PCA_unique$RAD_fourpop[tmp_overlap]
merged_micro_morph_PCA_unique$ClumppK4.0.50.cutoff[tmp_overlap]

new_morph$Lab.ID %in% merged_micro_morph_PCA_unique$Lab.ID

write.csv(merged_micro_morph_PCA_unique, 
            file="analysis/micro_morph_PCA_plotting.csv")

none_data <- merged_micro_morph_PCA_unique[is.na(merged_micro_morph_PCA_unique$ClumppK4.0.50.cutoff), ]
other_data <- merged_micro_morph_PCA_unique[!(is.na(merged_micro_morph_PCA_unique$ClumppK4.0.50.cutoff)), ]
table(other_data$sixpop)

#pccoord_subset_CG <-pccoord_subset[pccoord_subset$fourpop == "Coastal_Gulf" ,]
#pccoord_subset_CA <-pccoord_subset[pccoord_subset$fourpop == "Coastal_Atlantic" ,]
# plot na points first, then others

fill_colors <- c(
  "Coastal\nAtlantic" = "#4782d4",
  "Coastal\nGulf" = "#e1526b",
  "Intermediate\nAtlantic" = "#B4ED50",
  "Intermediate\nGulf" = "#2E8B57", 
  "Offshore\nAtlantic" = "#FFDD33", 
  "Offshore\nGulf" = "#C49E45"
)

none_data$group <- "5"  

p <- ggplot() +
  geom_point(data = none_data, aes(x=PC1, y=PC2, fill=group, shape=group, size="microsatellite"), color="grey") +
  geom_point(data = other_data, aes(x=PC1, y=PC2, fill=ClumppK4.0.50.cutoff, shape=ClumppK4.0.50.cutoff, size="microsatellite")) +
  #geom_point(data=pccoord_subset_CG, aes(x=PC1, y=PC2, shape=RAD_fourpop, fill=RAD_fourpop, size="RADseq"),
  #           color="black", shape=21, fill="#e1526b") +
  #geom_point(data=pccoord_subset_CA, aes(x=PC1, y=PC2, shape=RAD_fourpop, fill=RAD_fourpop, size="RADseq"),
  #           color="black", shape=21, fill="#4782d4") +
  scale_fill_manual(values = c( "1" = "#FFDD33", "2" = "#e1526b", "3" = "#4782d4", "4" = "#61BA5C", "5" = "grey"),
                    labels = c("1" = "Offshore", "2" = "Coastal Gulf", "3" = "Coastal Atlantic", "4" = "Intermediate", "5" = "Not genotyped"),
                    breaks = c("3", "2", "4", "1", "5"),
                    name = "Genotype group") +
  scale_shape_manual(values=c("1"=24, "2"=21, "3"=21, "4"=22, "5"=16),
                     labels = c("1" = "Offshore", "2" = "Coastal Gulf", "3" = "Coastal Atlantic", "4" = "Intermediate", "5" = "Not genotyped"),
                     breaks = c("3", "2", "4", "1", "5"),
                     name = "Genotype group") +
  scale_size_manual(values = c("microsatellite" = 3, "RADseq" = 6),
                    name = "Genotype method") +
  theme_classic() +
  guides(fill = guide_legend(override.aes = list(shape = c(21,21,22,24,16))),
         shape = guide_legend(override.aes = list(fill = c("#4782d4", "#e1526b", "#61BA5C", "#FFDD33", "grey"),
                                                  size=4)),
         size = guide_legend(override.aes = list(fill = "black", color = "black", shape = 21)))

p

ggsave("figures/morphology_pca_micros.png", p, h=4, w=6)
ggsave("figures/morphology_pca_micros.pdf", p, h=2.75, w=4.5)


pccoord_subset


#---------------------------------------------------------
# make it the 6 pop, equivalent to the genetics pca.

dat <- read.csv("metadata_FINAL.csv")
library(sf)
library(rnaturalearth)
library(ggspatial)

# First need to make map to make sure I'm assigning things to the correct location:
# Get the world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# read in locations
usa <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))

#----------------------------------------------
# add pop colors

# coastal
#56B4E9
#004488

#Offshore:
#F0B800
#B65A00

#intermediate:
#1B9E77
#66A61E
other_data$region <- "Atlantic"
other_data$region[other_data$Long < -81] <- "Gulf"
other_data$fourpop <- NA
other_data$fourpop[other_data$ClumppK4.0.50.cutoff =="1"] <- "Offshore"
other_data$fourpop[other_data$ClumppK4.0.50.cutoff =="2"] <- "Coastal_Gulf"
other_data$fourpop[other_data$ClumppK4.0.50.cutoff =="3"] <- "Coastal_Atlantic"
other_data$fourpop[other_data$ClumppK4.0.50.cutoff =="4"] <- "Intermediate"
  
other_data$sixpop <- paste(other_data$fourpop, other_data$region, sep= "_")
table(other_data$sixpop)

other_data$sixpop <- gsub("Coastal_Atlantic_Atlantic", "Coastal_Atlantic",other_data$sixpop)
other_data$sixpop <- gsub("Coastal_Gulf_Gulf", "Coastal_Gulf",other_data$sixpop)

ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_point(data = other_data, 
             aes(x = Long, y = Lat, color=sixpop, shape=sixpop),
             #aes(x = Long, y = Lat),
             size = 2.5,
             alpha=1) +
  coord_sf() +
  theme_bw() +
  theme(
    #panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    #axis.text = element_blank(),
    #axis.ticks = element_blank()
    legend.position = "left",
    legend.title=element_blank()) +
  xlab("Longitude")+
  ylab("Latitude") +
  #scale_shape_manual(values=c(21,22,23,24)) +
  #scale_fill_manual(values=c("#56B4E9","#004488","#66A61E","#F0B800")) +
  #scale_shape_manual(values=c(21,21,22,22,24,24))+
  #scale_fill_manual(values=c("#4782d4", "#e1526b","#B4ED50","#2E8B57", "#FFDD33", "#C49E45")) +
  coord_sf(xlim = c(-100, -70), ylim = c(22, 43), expand = FALSE)+
  annotation_scale()


# ggsave("figures/map_allpops.pdf", p1, h=4, w=5)
# ggsave("figures/map_allpops.png", p1, h=4, w=5)

# flip the x axis so it corresponds to genetics pca

none_data$PC1 <- none_data$PC1*-1
other_data$PC1 <- other_data$PC1*-1

p_m_sixpop <- ggplot() +
  geom_point(data = none_data, aes(x=PC1, y=PC2, fill=group, shape=group), color="grey",size = 2.5) +
  geom_point(data = other_data, aes(x=PC1, y=PC2, fill=sixpop, shape=sixpop),color="black",size = 2.5) +
  #geom_point(data=pccoord_subset_CG, aes(x=PC1, y=PC2, shape=RAD_fourpop, fill=RAD_fourpop, size="RADseq"),
  #           color="black", shape=21, fill="#e1526b") +
  #geom_point(data=pccoord_subset_CA, aes(x=PC1, y=PC2, shape=RAD_fourpop, fill=RAD_fourpop, size="RADseq"),
  #           color="black", shape=21, fill="#4782d4") +
  scale_fill_manual(values = c( "Coastal_Atlantic" = "#4782d4", "Coastal_Gulf" = "#e1526b", 
                                "Intermediate_Atlantic" = "#B4ED50", "Intermediate_Gulf" = "#2E8B57", 
                                "Offshore_Atlantic" = "#FFDD33",
                                "5" = "grey"),
                    labels = c("Coastal_Atlantic" = "Coastal Atlantic", "Coastal_Gulf" = "Coastal Gulf", 
                    "Intermediate_Atlantic" = "Intermediate Atlantic", "Intermediate_Gulf" = "Intermediate Gulf", 
                    "Offshore_Atlantic" = "Offshore Atlantic",
                    "5" = "Not genotyped"),
                    breaks = c("Coastal_Atlantic", "Coastal_Gulf", "Intermediate_Atlantic", "Intermediate_Gulf",
                               "Offshore_Atlantic", "5"),
                    name = "Genotype group") +

  scale_shape_manual(values=c("Coastal_Atlantic"=21, "Coastal_Gulf"=21, 
                              "Intermediate_Atlantic" = 22, "Intermediate_Gulf" = 22, 
                              "Offshore_Atlantic" = 24,
                              "5"=16),
                     labels = c("Coastal_Atlantic" = "Coastal Atlantic", "Coastal_Gulf" = "Coastal Gulf", 
                                "Intermediate_Atlantic" = "Intermediate Atlantic", "Intermediate_Gulf" = "Intermediate Gulf", 
                                "Offshore_Atlantic" = "Offshore Atlantic",
                                "5" = "Not genotyped"),                    
                     breaks = c("Coastal_Atlantic", "Coastal_Gulf", "Intermediate_Atlantic", "Intermediate_Gulf",
                                "Offshore_Atlantic", "5"),
                     name = "Genotype group") +
  #scale_size_manual(values = c("microsatellite" = 3, "RADseq" = 6),
  #                  name = "Genotype method") +
  theme_classic() +
  guides(fill = guide_legend(override.aes = list(shape = c(21,21,22,22,24,16))),
         shape = guide_legend(override.aes = list(fill = c("#4782d4", "#e1526b","#B4ED50", "#2E8B57", "#FFDD33", "grey"),
                                                  size=4)),
         size = guide_legend(override.aes = list(fill = "black", color = "black", shape = 21)))

p_m_sixpop

#ggsave("figures/morphology_pca_micros.png", p, h=4, w=6)
#ggsave("figures/morphology_pca_micros.pdf", p, h=2.75, w=4.5)
