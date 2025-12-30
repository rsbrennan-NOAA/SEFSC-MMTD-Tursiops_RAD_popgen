library(dplyr)
# merge data for supplemental summary table:

df_dist <- read.csv("analysis/distance.csv")

df_depth <- read.csv("depths.csv")

df_depthDist <- merge(df_dist, df_depth, by="id")
#id

df_meta <- read.csv("Tursiops_RADseq_Metadata_new_forMerge.csv")
# id

# pop assignments
df_hybs <- read.table("analysis/population_assignments_hybrids_summary.txt", header=T)
# indiv

final_df <- df_meta %>%
  left_join(df_depth, by = "id") %>%
  left_join(df_dist, by = "id") %>%
  left_join(df_hybs, by = c("id" = "indiv"))

head(final_df)
nrow(final_df)
nrow(df_meta)


########################## ---------------------------------------
# add ncbi location/id to this table

# ncbi accessions

new <- read.csv("indiv_fastq_summary.txt", sep="\t")# this is all files that were sequenced.
new$IndividualID
nc_acc <- read.csv("NC_accessions.csv") # these are the NC samples already uploaded. 

#new$IndividualID <- gsub("b$", "", new$IndividualID)  # Remove "b" at the end
#new$IndividualID <- gsub("-rep", "", new$IndividualID)  # Remove "-rep"

# of the new files, which are in the old NC samples. 
sum(new$IndividualID %in% nc_acc$sample_name)
#142
nrow(nc_acc)

#142
nc_acc$sample_name[! nc_acc$sample_name %in% new$IndividualID]

#should be 142, but only getting 137. because I dropped some. 
gsub("b$", "",final_df$id) %in% nc_acc$sample_name


ids <- gsub("b$","", new$IndividualID)
nc_acc$sample_name[! nc_acc$sample_name %in% ids]

final_df$NCBI.Accession[gsub("b$", "",final_df$id) %in% nc_acc$sample_name] <- "PRJNA1227152"
final_df$NCBI.Accession[!gsub("b$", "",final_df$id) %in% nc_acc$sample_name] <- "PRJNA1289198"


# non NC: # add bioproject to final_df PRJNA1289198
# NC PRJNA1227152

#208 are

final_df$species <- ifelse(final_df$fourpop %in% c("Coastal_Atlantic", "Coastal_Gulf", "Intermediate"), 
                              "Tursiops erebennus", 
                              "Tursiops truncatus")

head(final_df)


write.csv(final_df, file="metadata_FINAL.csv", row.names=F, quote = F)


# add species:
noNC_merged <- merge(final_df, noNC, by.x="id", by.y = "IndividualID")
nrow(noNC_merged)
nrow(noNC)
noNC_merged$species <- ifelse(noNC_merged$fourpop %in% c("Coastal_Atlantic", "Coastal_Gulf", "Intermediate"), 
                              "Tursiops erebennus", 
                              "Tursiops truncatus")

write.csv(noNC_merged, file="ncbi_toadd.csv", row.names=F)

accessions <- read.delim("ncbi_accessions.tsv")

out<- merge(noNC_merged, accessions, by.x="id", by.y="sample_name")
write.csv(out, file="ncbi_toadd_accessions.csv", row.names=F)



!!!!!!!!!!!!!!!!
# add actual accessions and bioproject ids to the metadata
!!!!!!!!!!!!!!!!

# new ones in
  ## bioproject PRJNA1289198
accessions

# old ones:
nc_acc





# make sure that all indivs in the metadata are in the existing NCBI accessions and to the to add:
noNC$IndividualID[!noNC$IndividualID %in% final_df$id]

noNC$IndividualID[noNC$IndividualID %in% nc_acc$sample_name]
nc_acc$sample_name[nc_acc$sample_name %in% noNC$IndividualID]

sum(noNC$IndividualID %in% final_df$id)
#208
sum(nc_acc$sample_name %in% final_df$id)
#137






