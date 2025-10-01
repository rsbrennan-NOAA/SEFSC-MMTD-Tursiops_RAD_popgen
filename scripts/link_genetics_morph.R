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


#--------------------------------------------------
# link to RAD id:
#--------------------------------------------------

# from out_match, figure out which have RAD genotypes:
dat <- read.csv("metadata_FINAL.csv")

dat$id <- gsub("b$", "", dat$id)
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

# add the micro assignments from nikki:

micro_data <- read.csv("analysis/RADseqMicroData_n343_2025_NEW-RAD&MicroStructureData_n343.csv")

extract_ids <- function(lab_id_string) {
  # Split by "="
  ids <- trimws(unlist(strsplit(as.character(lab_id_string), "=")))
  return(ids)
}

# Create expanded dataframe with all possible ID matches
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


micro_matched_RAD$combined_category <- paste(micro_matched_RAD$fourpop.y, 
                                             micro_matched_RAD$newhybrids_category, 
                                             sep = "_")
micro_matched_RAD$combined_category <- gsub("_NA", "",micro_matched_RAD$combined_category)

table(micro_matched_RAD$combined_category, 
                              micro_matched_RAD$ClumppK4.0.50.cutoff, 
                              useNA = "ifany")



#---------------------------------------------------------
# now read in entire microsat dataset:
#---------------------------------------------------------
allmicro <- read.csv("analysis/RADseqMicroData_n343_2025_NEW-MicroStructureData_n4189.csv")

allmicro <- cbind(allmicro$Lab.ID, allmicro)
colnames(allmicro) <- c("ID.Original",colnames(allmicro)[2:ncol(allmicro)])

# Create empty dataframe to store expanded rows
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

#----------------------------------------------------
# merge with morphology data
matched_df <-cbind(matched_df$Lab_ID, matched_df)
colnames(matched_df) <- c("ID.Original",colnames(matched_df)[2:ncol(matched_df)])

matched_df_expanded <- matched_df[0, ]
# Loop through each row
for(i in 1:nrow(matched_df)) {
  current_row <- matched_df[i, ]
  
  # Get Lab_ID value
  lab_id_value <- as.character(current_row$Lab_ID)
  
  # Handle NA or empty values
  if(is.na(lab_id_value) || lab_id_value == "") {
    matched_df_expanded <- rbind(matched_df_expanded, current_row)
    next
  }
  
  # Check if Lab_ID contains parentheses
  if(grepl("\\(", lab_id_value)) {
    # Extract the part before parentheses and the part inside parentheses
    before_paren <- trimws(sub("\\(.*\\)", "", lab_id_value))
    inside_paren <- trimws(gsub(".*\\((.*)\\).*", "\\1", lab_id_value))
    
    # Create first row with value before parentheses
    new_row1 <- current_row
    new_row1$Lab_ID <- before_paren
    matched_df_expanded <- rbind(matched_df_expanded, new_row1)
    
    # Create second row with value inside parentheses
    new_row2 <- current_row
    new_row2$Lab_ID <- inside_paren
    matched_df_expanded <- rbind(matched_df_expanded, new_row2)
    
  } else if(grepl("=", lab_id_value)) {
    # Contains "=" but no parentheses
    splits <- trimws(unlist(strsplit(lab_id_value, "=")))
    
    for(split_val in splits) {
      new_row <- current_row
      new_row$Lab_ID <- split_val
      matched_df_expanded <- rbind(matched_df_expanded, new_row)
    }
    
  } else {
    # No parentheses or "=", just add the row as is
    matched_df_expanded <- rbind(matched_df_expanded, current_row)
  }
}


head(matched_df_expanded)
nrow(matched_df_expanded)

merged_micro_morph <- merge(allmicro_expanded, matched_df_expanded, 
                   by.x = "Lab.ID", by.y = "Lab_ID", 
                   all = TRUE)

merged_micro_morph[!is.na(merged_micro_morph$ClumppK4.0.50.cutoff),]

merged_micro_morph_PCA <- merge(merged_micro_morph, pccoord, by.x="morph_Field_ID", by.y="sample_name")

# remove the dups from the previous step:
merged_micro_morph_PCA_unique <- merged_micro_morph_PCA[!duplicated(merged_micro_morph_PCA$morph_Field_ID), ]

nrow(merged_micro_morph_PCA_unique)
# 180

none_data <- merged_micro_morph_PCA_unique[is.na(merged_micro_morph_PCA_unique$ClumppK4.0.50.cutoff), ]
other_data <- merged_micro_morph_PCA_unique[!(is.na(merged_micro_morph_PCA_unique$ClumppK4.0.50.cutoff)), ]

pccoord_subset_CG <-pccoord_subset[pccoord_subset$fourpop == "Coastal_Gulf" ,]
pccoord_subset_CA <-pccoord_subset[pccoord_subset$fourpop == "Coastal_Atlantic" ,]
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
  scale_fill_manual(values = c( "1" = "#FFDD33", "2" = "#e1526b", "3" = "#4782d4", "4" = "#B4ED50", "5" = "grey"),
                    labels = c("1" = "Offshore", "2" = "Coastal Gulf", "3" = "Coastal Atlantic", "4" = "Intermediate", "5" = "Not genotyped"),
                    breaks = c("3", "2", "4", "1", "5"),
                    name = "Genotype group") +
  scale_shape_manual(values=c("1"=24, "2"=21, "3"=21, "4"=22, "5"=16),
                     labels = c("1" = "Offshore", "2" = "Coastal Gulf", "3" = "Coastal Atlantic", "4" = "Intermediate", "5" = "Not genotyped"),
                     breaks = c("3", "2", "4", "1", "5"),
                     name = "Genotype group") +
  scale_size_manual(values = c("microsatellite" = 3, "RADseq" = 6),
                    name = "Genotype method") +
  theme_classic(base_size=14) +
  guides(fill = guide_legend(override.aes = list(shape = c(21,21,22,24,16))),
         shape = guide_legend(override.aes = list(fill = c("#4782d4", "#e1526b", "#B4ED50", "#FFDD33", "grey"),
                                                  size=4)),
         size = guide_legend(override.aes = list(fill = "black", color = "black", shape = 21)))

p

ggsave("figures/morphology_pca_micros.png", p, h=4, w=6)


pccoord_subset

#------------------------------------------------------------------------#------------------------------------------------------------------------#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
# get missing microsats
#------------------------------------------------------------------------

# read in all microsat data:

library(googledrive)
library(readxl)
library(httr)

drive_auth()

folder_id <- "1-kVGWCDQo9jazFTabcRz-7zA33YX5qId"

# all excel files
xlsx_files <- drive_ls(path = as_id(folder_id), pattern = "\\.xlsx$")

print(paste("Found", nrow(xlsx_files), "Excel files:"))
print(xlsx_files$name)

excel_data_list <- list()

temp_dir <- tempdir()

# Loop through each file
for(i in 1:nrow(xlsx_files)) {
  file_name <- xlsx_files$name[i]
  file_id <- xlsx_files$id[i]
  
  cat("Processing file", i, "of", nrow(xlsx_files), ":", file_name, "\n")
  
  # make temp file path
  temp_path <- file.path(temp_dir, file_name)
  
  drive_download(as_id(file_id), path = temp_path, overwrite = TRUE)
  
  data <- read_excel(temp_path, sheet = 1)
  
  # add to list
  clean_name <- gsub("\\.xlsx$", "", file_name)
  excel_data_list[[clean_name]] <- data
  
  file.remove(temp_path)
}

cat(length(excel_data_list), "excel files read into list")


# add the DWH files:

folder_id <- "1vHWeOdqblht3HTQKWB1YKZLEkI3TnqHb"

# all excel files
xlsx_files <- drive_ls(path = as_id(folder_id), pattern = "\\.xlsx$")
print(paste("Found", nrow(xlsx_files), "Excel files:"))
print(xlsx_files$name)

# Loop through each file
for(i in 1:nrow(xlsx_files)) {
  file_name <- xlsx_files$name[i]
  file_id <- xlsx_files$id[i]
  
  cat("Processing file", i, "of", nrow(xlsx_files), ":", file_name, "\n")
  
  # Create temporary file path
  temp_path <- file.path(temp_dir, file_name)
  
  # Download file
  drive_download(as_id(file_id), path = temp_path, overwrite = TRUE)
  
  # Read first sheet only
  data <- read_excel(temp_path, sheet = 1)
  
  # Store in list
  clean_name <- gsub("\\.xlsx$", "", file_name)
  excel_data_list[[clean_name]] <- data
  
  # Clean up temp file
  file.remove(temp_path)
}

cat(length(excel_data_list), "excel files read into list")

# need to expand the matches:

excel_data_expanded <- list()

for(file_name in names(excel_data_list)) {
  df <- excel_data_list[[file_name]]
  lab_col <- grep("lab", colnames(df), ignore.case = TRUE)
  df <- cbind(df[lab_col], df)
  colnames(df) <- c("id_original", colnames(df)[2:ncol(df)])
  # Create empty dataframe to store expanded rows
  df_expanded <- df[0, ]
  
  # Loop through each row and expand lab column
  for(j in 1:nrow(df)) {
    current_row <- df[j, ]
    lab_value <- as.character(current_row[, lab_col])
    # Split lab_value by either parentheses or "="
    if(grepl("\\(", lab_value)) {
      # Extract values before and inside parentheses
      before_paren <- trimws(sub("\\(.*\\)", "", lab_value))
      inside_paren <- trimws(gsub(".*\\((.*)\\).*", "\\1", lab_value))
      
      # Create two rows
      new_row1 <- current_row
      new_row1[, lab_col+1] <- before_paren
      df_expanded <- rbind(df_expanded, new_row1)
      
      new_row2 <- current_row
      new_row2[, lab_col+1] <- inside_paren
      df_expanded <- rbind(df_expanded, new_row2)
      
    } else if(grepl("=", lab_value)) {
      # Split by "="
      splits <- trimws(unlist(strsplit(lab_value, "=")))
      for(split_val in splits){
        new_row <- current_row
        new_row[, lab_col+1] <- split_val
        df_expanded <- rbind(df_expanded, new_row)
        }
      } else {
      # No special characters, use as is
      df_expanded <- rbind(df_expanded, current_row)
      }
  }
  
  excel_data_expanded[[file_name]] <- df_expanded
  print(paste("done with ",file_name))
}


nrow(excel_data_expanded[[1]])
nrow(excel_data_list[[1]])
excel_data_expanded[[1]][1:5,1:5]

# use matched_df_expanded
matching_rows <- list()
for(i in 1:nrow(matched_df_expanded)){
  for(file_name in names(excel_data_expanded)) {
    df <- excel_data_expanded[[file_name]]
    lab_col <- grep("lab", colnames(df), ignore.case = TRUE)
    
    tmp_id_search <- matched_df_expanded$morph_Lab_ID[i]
    matching_rows_tmp <- which(df[, lab_col] == tmp_id_search)
    
    #if(length(matching_rows_tmp) == 1){
    #  matching_rows[[file_name]] <- rbind(matching_rows[[file_name]], df[matching_rows_tmp,])
    #}
    if(length(matching_rows_tmp) > 0){
      matching_rows[[file_name]] <- rbind(matching_rows[[file_name]], df[matching_rows_tmp,])
      }
  }
  print(i)
}


# summary
if(length(matching_rows) > 0) {
  cat("\nSummary by file:\n")
  for(file_name in names(matching_rows)) {
    cat(file_name, ":", nrow(matching_rows[[file_name]]), "match\n")
  }
}
cat("Total across all files:", 
    sum(sapply(matching_rows, nrow)))
# Total across all files: 122

# Show the actual matching rows for each file
for(file_name in names(matching_rows)) {
  cat("- ", file_name, " -\n")
  print(matching_rows[[file_name]])
  cat("\n")
}

find_loci_column <- function(df) {
  loci_col <- grep("Loci", colnames(df), ignore.case = TRUE)
  if(length(loci_col) > 0) {
    return(loci_col[1])  # Return first match if multiple
  } else {
    return(NULL)
  }
}

# show results for each file
standardized_matches <- list()

for(file_name in names(matching_rows)) {
  df <- matching_rows[[file_name]]
  
  # Get first 8 columns
  first_8 <- df[, 1:3]

  # Find Loci column
  loci_col_num <- find_loci_column(df)
  result <- cbind(first_8, df[, loci_col_num, drop = FALSE])
  colnames(result)<- c("id_origina","lab_id", "field_id", "# of Loci")
  print(result)
  
  # Store for combining
  standardized_matches[[file_name]] <- result
}

# Now rbind
combined_matches <- do.call(rbind, standardized_matches)

nrow(combined_matches)
# 123
    
head(combined_matches, n=20)
sum(combined_matches$`# of Loci` > 0)
# 103

sum(combined_matches$`# of Loci` >18)
#98
sum(combined_matches$`# of Loci` ==19)
#87
sum(combined_matches$`# of Loci` >19)
#11


# need to pull out microsat loci for analysis. write to new csv


# show results for each file
standardized_matches <- list()

for(file_name in names(matching_rows)) {
  df <- matching_rows[[file_name]]
  
  # Get first 8 columns
  first_2 <- df[, 1:3]
  
  # lab id is always first col
  # haplotype col:
  hap_col <- grep("haplotype", colnames(df), ignore.case = TRUE)
  source_col <- grep("Source", colnames(df), ignore.case = TRUE)
  col_date_col <- grep("Collection", colnames(df), ignore.case = TRUE)
  microstart_col <- grep("Ttr", colnames(df), ignore.case = TRUE)[1]
  
  
  # Find Loci column
  loci_col_num <- find_loci_column(df)
  result <- cbind(first_2, df[, c(col_date_col,source_col,loci_col_num,hap_col,c(microstart_col:ncol(df))), drop = FALSE])
  colnames(result)<- c("original_id", "lab_id", "field_id", "Collection_date", "Source", "# of Loci",
                       "Haplotype", colnames(df)[microstart_col:ncol(df)])
  #print(result)
  result$'# of Loci' <- as.character(result$'# of Loci')
  result$Collection_date <- as.character(result$Collection_date)
  microstart_result_col <- which(colnames(result) == colnames(df)[microstart_col])
  result[, microstart_result_col:ncol(result)] <- lapply(result[, microstart_result_col:ncol(result)], as.character)
  
  
  # Store for combining
  standardized_matches[[file_name]] <- result
}


micro_out_morph <- list_rbind(standardized_matches)
nrow(micro_out_morph)
sum(micro_out_morph$`# of Loci` =="19")
micro_out_morph[micro_out_morph$`# of Loci` =="19",]
head(micro_out_morph)
micro_out_morph[1:50,1:5]
# the problem here, is that I have multiple matches because some indivs were sampled multiple times. 
# how to parse it down? 
# I can have either duplicate lab ids, or duplicate field ids. 
# because lab1=lab2 OR lab2=lab1. So these go into two rows. 
# So I want the one that is the lab id from Ana AND was the stranding. 

# first easy step, just find the ones where lab ids match:
nrow(micro_out_morph)
#122
micro_out_morph_f1 <- micro_out_morph[(micro_out_morph$lab_id %in% pccoord$LabID |
                                         micro_out_morph$field_id %in% pccoord$sample_name),]
nrow(micro_out_morph_f1)
#108

micro_out_morph_f1[1:15,1:3]

micro_out_morph2 <- micro_out_morph %>%
  group_by(lab_id, field_id, Collection_date, Source) %>%
  slice(1) %>%  # Keep only first occurrence of each unique combination
  ungroup()


write.table(micro_out_morph,
            file="analysis/morphology_microsats.txt", 
            sep="\t", quote = FALSE, row.names=FALSE)

micro_out_morph[micro_out_morph$`# of Loci` == 14,]



29Tt101
26Tt338

allsamp[grep("29Tt100",allsamp$Lab_ID),]
allsamp[grep("4Tt037",allsamp$Lab_ID),]

length(unique(micro_out_morph$lab_id))
# 112
length(unique(micro_out_morph$field_id))
#109
df_step1 <- micro_out_morph[!duplicated(micro_out_morph$lab_id), ]
df_final <- df_step1[!duplicated(df_step1$field_id), ]
nrow(df_final)
#105

tmp1 <- micro_out_morph[unique(micro_out_morph$lab_id) %in% micro_out_morph$lab_id,]
tmp2 <- tmp1[unique(tmp1$field_id) %in% tmp1$field_id,]
nrow(tmp2)

micro_out_morph$lab_id == "29Tt101"
micro_out_morph$Source

micro_out_morph_f1 <- micro_out_morph[(micro_out_morph$lab_id %in% pccoord$LabID),]

micro_out_morph_f1[1:10,1:5]
sum(pccoord$LabID %in% micro_out_morph$lab_id)

nrow(micro_out_morph_f1)
length(unique(micro_out_morph_f1$lab_id))
length(unique(micro_out_morph_f1$field_id))

todrop <- which(micro_out_morph_f1$Source[duplicated(micro_out_morph_f1$lab_id)] != "stranding")

micro_out_morph_f2 <- micro_out_morph_f1[-todrop,]
nrow(micro_out_morph_f2)
# need to get rid of the duplicates. 
## this happens because there are duplicates in the micros. one indiv sampled mult times
## I want to keep the stranding one, 
## the () ones, I don't know what is going on.
micro_out_morph

length(unique(micro_out_morph_f2$lab_id))
todrop <- which(duplicated(micro_out_morph_f1$lab_id) & micro_out_morph_f1$Source != "stranding")
micro_out_morph_f2 <- micro_out_morph_f1[-todrop,]
nrow(micro_out_morph_f2)


write.table(micro_out_morph_f2,
            file="analysis/morphology_microsats_dedup.txt", 
            sep="\t", quote = FALSE, row.names=FALSE)

# figure out which duplicates to drop.


#--------------------------------------------
#--------------------------------------------
#--------------------------------------------
#--------------------------------------------
#--------------------------------------------
# from out_match, figure out which have RAD genotypes:
dat <- read.csv("metadata_FINAL.csv")

dat$id <- gsub("b$", "", dat$id)

rad_match <- c()
for(i in dat$id){
  tmp_match <- grep(i, out_match$Lab_ID)
  if(length(tmp_match) > 0){
    rad_match <- c(rad_match, i)
  }
}
rad_match
[1] "2Tt156" "2Tt700" "4Tt184"

# here were the original ones I found. why don't they overlap?
#[1] "USNM571388" "USNM572828" "LA813"  "LA823"   
# USNM571388  2Tt700
# USNM572828  2Tt156
# LA813 26Tt309
# LA823 26Tt313

# 26Tt313  and 26Tt309

allsamp[grep("26Tt313",allsamp$Lab_ID),]
# JSH20130429LA001
allsamp[grep("LA813",allsamp$Field_ID),]
allsamp[grep("26Tt309",allsamp$Lab_ID),]
allsamp[grep("4Tt184",allsamp$Lab_ID),]


pccoord$sample_name == "LA813"

dat[dat$id == "4Tt184",]

out_match$Catalog_num == "USNM572828"
out_match[grep("26Tt309", out_match$Field_ID),]
tmp_check <- apply(out_match, 1, function(x) any(grepl("LA813", x,
                                                       ignore.case = TRUE))) 
tmp_check <- apply(allsamp, 1, function(x) any(grepl("26Tt309", x,
                                                     ignore.case = TRUE))) 
allsamp[tmp_check,]

tmp_check <- apply(allsamp, 1, function(x) any(grepl("26Tt313", x,
                                                     ignore.case = TRUE))) 
allsamp[tmp_check,]

LA823 26Tt313





#-----
# try with Ana's ids directly:

pccoord$sample_name %in% morph$lab_id
morph$Field_id[sum(morph$Field_id %in% pccoord$sample_name)]

pccoord$lab_id <- NA
for(i in 1:nrow(pccoord)){
  tmp_mtch <- which(morph$Field_id == pccoord$sample_name[i])
  if(length(tmp_mtch) > 0){
    pccoord$lab_id[i] <-morph$lab_id[tmp_mtch] 
  }
}

head(pccoord$lab_id)
sum(!is.na(pccoord$lab_id))
# 88

rad_match <- c()
for(i in dat$id){
  tmp_match <- grep(i, out_match$Lab_ID)
  if(length(tmp_match) > 0){
    rad_match <- c(rad_match, i)
  }
}
out_match

sum(pccoord$lab_id %in% dat$id)

# next microsats:
# 88 should have matches


rad_match <- c()
for(i in dat$id){
  tmp_match <- grep(i, out_match$Lab_ID)
  if(length(tmp_match) > 0){
    rad_match <- c(rad_match, i)
  }
}
out_match


pccoord$sample_name[!pccoord$sample_name %in% combined_matches$field_id]

length(pccoord$sample_name[!pccoord$sample_name %in% combined_matches$field_id])
# 126

pccoord$lab_id_1


sum(as.numeric(morph$Loci) > 0 & morph$TL > 0 & !(is.na(morph$TL)))
                                                  
sum(as.numeric(morph$Loci) > 0, na.rm=T)

pccoord$lab_id_1


# in ana's table: USNM504121 2Tt697
allsamp[grep("7Tt497",allsamp$Lab_ID),]
# 7Tt497  USNM504656
#!!!!!!!!! they don't match
allsamp[grep("USNM504121",allsamp$Field_ID),]

# in ana's table: USNM504482 9Tt168
allsamp[grep("9Tt168",allsamp$Lab_ID),]
# 7Tt497  USNM504656
#!!!!!!!!! they don't match
allsamp[grep("USNM504482",allsamp$Field_ID),]


# Catalog_num is sometimes used instead of field number in ana's table. 
# so then I can't directly compare even to our RAD data. 


























