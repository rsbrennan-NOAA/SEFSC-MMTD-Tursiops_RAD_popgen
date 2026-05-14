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

final_df <- df_meta |> 
  left_join(df_depth, by = "id") |> 
  left_join(df_dist, by = "id") |> 
  left_join(df_hybs, by = c("id" = "indiv"))

head(final_df)
nrow(final_df)
nrow(df_meta)


########################## ---------------------------------------
# add ncbi location/id to this table

# ncbi accessions

all <- read.csv("indiv_fastq_summary.txt", sep="\t")# this is all files that were sequenced.

all2 <- all[grep("^41",all$IndividualID, invert=T),]
nrow(all2)
#385
all3 <- all2[grep("^157",all2$IndividualID, invert=T),]
nrow(all3)
#380

# drop Tadu
all4 <- all3[grep("^Tadu",all3$IndividualID, invert=T),]
nrow(all4)

#[1] "42193" "78068"
all5 <- all4[grep("42193|78068",all4$IndividualID, invert=T),]
nrow(all5)
#375

nc_acc <- read.csv("NC_accessions.csv") # these are the NC samples already uploaded. 

new$IndividualID <- gsub("b$", "", new$IndividualID)  # Remove "b" at the end
#new$IndividualID <- gsub("-rep", "", new$IndividualID)  # Remove "-rep"
nc_acc$sample_name <- gsub("b$", "", nc_acc$sample_name)
# of the new files, which are in the old NC samples. 
sum(new$IndividualID %in% nc_acc$sample_name)
#142
nrow(nc_acc)

new$IndividualID[!new$IndividualID %in% final_df$id]
#142
nc_acc$sample_name[! nc_acc$sample_name %in% new$IndividualID]

#should be 142, but only getting 137. because I dropped some. 
# make sure:
dat <- read.csv("analysis/out.imiss", 
                header=T, sep="\t")
tout <- table(gsub("b","",dat$INDV))
tout[order(tout)]
dropped_sm <- dat[which(dat$F_MISS > 0.75),]

gsub("b$", "",final_df$id) %in% nc_acc$sample_name

sum(new$IndividualID %in% nc_acc$sample_name)

#indivs that were in original nikki data
all_sm <- gsub("b$", "",new$IndividualID[new$IndividualID %in% nc_acc$sample_name])

#indivs in current filtered data
included_sm <- gsub("b$", "",final_df$id[gsub("b$", "",final_df$id) %in% nc_acc$sample_name])

miss_sm <- all_sm[!all_sm %in% included_sm]

# now check that the missing ones were ones I dropped:
dropped_sm[dropped_sm$INDV%in% miss_sm,]
# yes, all 5 were dropped bc high missing data

final_df$NCBI.Bioproject[gsub("b$", "",final_df$id) %in% nc_acc$sample_name] <- "PRJNA1227152"
final_df$NCBI.Bioproject[!gsub("b$", "",final_df$id) %in% nc_acc$sample_name] <- "PRJNA1289198"

# non NC: # add bioproject to final_df PRJNA1289198
# NC PRJNA1227152

#208 are

final_df$species <- ifelse(final_df$fourpop %in% c("Coastal_Atlantic", "Coastal_Gulf", "Intermediate"), 
                              "Tursiops erebennus", 
                              "Tursiops truncatus")

head(final_df)

# add mt haplotype:

grouping <- read_csv("analysis/trees/RADseq_mtDNA_labels.csv")

final_df2 <- final_df %>%
  left_join(
    grouping %>% select(`indiv`, Haplotype),
    by = c("id" = "indiv")
  ) %>%
  relocate(Haplotype, .after = Collection.Date)

write.csv(file="meta_table_forNCBI.csv", final_df2, quote = F, row.names=F)


# write fasta for submission to ncbi:

genbank <- read.csv("AllATL_mtDNA_SeqInfo_metadata.csv")

genbank[51,]
# need to get species for each

genbank_meta <- final_df2 %>%
  right_join(
    genbank,
    by = c("id" = "indiv")
  ) 

#drop the non-sequenced one:
genbank_meta <- genbank_meta[genbank_meta$Haplotype != "NA-Never Sequenced",]

nrow(genbank_meta)
#203 to add

#>SeqID [organism=Tursiops sp] [HAPLOTYPE] mitochondrial control region

for(i in 1:nrow(genbank_meta)){
  tmpin <- genbank_meta[i,]
  line1 <- paste0(">", tmpin$id, " [organism=", tmpin$species, "] [haplotype=", 
                  tmpin$Haplotype, 
                  "] ",tmpin$species, 
                  " tRNA-Pro gene, partial sequence; control region, partial sequence; mitochondrial.")
  line2 <- tmpin$Sequence
  cat(line1, line2, file = "tursiops_submit.fasta", sep = "\n", append = TRUE)
}

write.csv(file="genbank_metadata_toSubmit.csv", genbank_meta, quote = F, row.names=F)


# write feature table:


# control region starts at: "GAAAAA"
# CAACCC is the end of the haplotype. Include this in the length. 
# this needs to be in the source file, not feature. 
# just check for now that they're all the same

for(i in 1:nrow(genbank_meta)){
  tmpin <- genbank_meta[i,]
  seqlen <- nchar(tmpin$Sequence)
  title_ln <- paste(">Feature", tmpin$id, sep=" ")
  tRNA_end <- regexpr("GAAAAA", tmpin$Sequence, ignore.case = T)[1]-1
  total_len <- paste0(">",seqlen)
  feat_name <- "D-loop"
  qual_name <- "note"
  qual_val <- "control region"
  line1_1 <- paste(tRNA_end, ">1", "tRNA", "", "", sep="\t")
  line2_1 <- paste("", "", "", "product", "tRNA-Pro", sep="\t")
  line1 <- paste(tRNA_end+1, total_len, feat_name, "", "", sep="\t")
  line2 <- paste("", "", "", qual_name, qual_val, sep="\t")
  cat(title_ln,
      line1_1, line2_1,
      line1, line2, file = "tursiops_submit.feature_table.txt", sep = "\n", append = TRUE)
  cat(print(regexpr("CAACCC", tmpin$Sequence, ignore.case = T)[1] + 6 - as.numeric(tRNA_end +1)),
      file = "tursiops_submit.CR_Hap_len.txt", sep = "\n", append = TRUE)
}







# add mtdna genbank info:

mtbank <- read.csv("Genbank_mt_seqIDs.txt", sep="\t", header=F)
colnames(mtbank) <- c("indiv", "GenBank")

dim(mtbank)
head(mtbank)
# need to get the haplotype info and add it to mtbank

genbank <- read.csv("AllATL_mtDNA_SeqInfo.csv")
dim(genbank)
head(genbank)

genbank$GenBank <- ifelse(
  genbank$GenBank =="" & genbank$indiv %in% mtbank$indiv,
  mtbank$GenBank[match(genbank$indiv, mtbank$indiv)],
  genbank$GenBank
)

genbank[which(genbank$indiv == "FB196"),]
sum((genbank$GenBank == ""))


# add the genbank to the df:

final_df3 <- final_df2 %>%
  left_join(
    genbank %>% select(indiv,GenBank),
    by = c("id" = "indiv")
  ) %>%
  relocate(GenBank, .after = Haplotype)

colnames(final_df3)[colnames(final_df3) == "Haplotype"] <- "mt_haplotype"
colnames(final_df3)[colnames(final_df3) == "GenBank"] <- "mt_GenBank"

final_df3[final_df3$Haplotype != final_df3$Haplotype2,]
nrow(final_df3)

final_df3[51,]

write.csv(final_df3, file="metadata_FINAL.csv", row.names=F, quote = F)

#final_df4 <- final_df3 %>% distinct(Haplotype, .keep_all = TRUE)
#final_df4[,1:5]

#write.csv(final_df4, file="metadata_mtGenbank.csv", row.names=F, quote = F)

