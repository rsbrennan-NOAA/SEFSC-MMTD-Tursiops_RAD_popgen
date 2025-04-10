#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
# write out files for treemix:
# need to have relatively equal sample sizes

# 1: 4 pops from admixture
# 2: separate gulf and atlantic for intermediate and offshore.


summary_dat <- read.table(file="analysis/population_assignments_summary.txt", header=T)

#------------------------------------------------------------
# 1: 4 pop
table(summary_dat$fourpop)
# Coastal_Atlantic     Coastal_Gulf     Intermediate         Offshore 
# 168               46               60               71 

write.table(file="analysis/pop_structure/fourpop_all.clust", 
            data.frame(samp1 = summary_dat$indiv,
                       group = summary_dat$fourpop),
            sep="\t", quote=F, col.names=F, row.names=F)


# We have huge overrepresentation from NC. so we want to target those indivs.
# we have 168 coastal ATL. But 118 are from NC
# Coastal_Atlantic     Coastal_Gulf     Intermediate         Offshore 
# 168               46               60               71 
datMeta <- read.csv("Tursiops_RADseq_Metadata_new.csv")

ncsamp <- datMeta$Lab.ID[datMeta$Lat > 33.32 & datMeta$Lat < 36.16]
length(ncsamp)
#136

# get just the ones that are Coastal_Atlantic

ncsamp2 <- ncsamp[ncsamp %in% summary_dat$indiv[summary_dat$dapc_population == "Coastal_Atlantic"]]
ncsamp <- ncsamp2
length(ncsamp)
# 118

# there are 50 non NC samples. 
# subsamp down to 20 NC samples. 

coastal_atl_subset_idx <- sample(ncsamp, 20)
coastal_atl_subset <-  summary_dat[summary_dat$indiv %in% coastal_atl_subset_idx,]
nrow(coastal_atl_subset)

other_dat <- summary_dat[!summary_dat$indiv %in% ncsamp,]
nrow(other_dat)

combodat <- rbind(coastal_atl_subset, other_dat)

# Verify new counts
table(combodat$fourpop)
#Coastal_Atl Coastal_Gulf Intermediate     Offshore 
#70           47           54           74 
sum(table(combodat$fourpop))
# 247

write.table(file="analysis/pop_structure/fourpop.clust", 
            data.frame(samp1 = combodat$indiv,
                       group = combodat$fourpop),
            sep="\t", quote=F, col.names=F, row.names=F)


###----------------------------------------------------------------------------
# 2: separate gulf and atlantic for intermeidate and offshore.

# add region to the pops

table(summary_dat$sixpop)
# Coastal_Atlantic          Coastal_Gulf Intermediate_Atlantic     Intermediate_Gulf 
# 168                    46                    44                    16 
# Offshore_Atlantic         Offshore_Gulf 
# 33                    38 

# write the output, 

write.table(file="analysis/pop_structure/sixpop_all.clust", 
            data.frame(samp1 = summary_dat$indiv,
                       group = summary_dat$sixpop),
            sep="\t", quote=F, col.names=F, row.names=F)



# Subset Coastal_Atl IDs, keep all others

#coastal_atl_subset_idx <- sample(ncsamp, 20)
#coastal_atl_subset2 <-  summary_dat[summary_dat$indiv %in% coastal_atl_subset_idx,]
nrow(coastal_atl_subset)

# then another 26 from the other coastal_atl that aren't the NC pop
ncsamp
summary_dat

cst_atl <- summary_dat$indiv[summary_dat$fourpop == "Coastal_Atlantic"]
length(cst_atl)
#168
cst_atl <- cst_atl[!cst_atl %in% ncsamp]
length(cst_atl)
cst_sub <- sample(cst_atl, 26)
length(cst_sub)
cst_sub_df <-  summary_dat[summary_dat$indiv %in% cst_sub,]
cst_merge <- rbind(cst_sub_df,coastal_atl_subset )
nrow(cst_merge)
#46

other_dat <- summary_dat[summary_dat$sixpop != "Coastal_Atlantic",]
nrow(other_dat)
#177
combodat <- rbind(cst_merge, other_dat)
nrow(combodat)
#223
# Verify new counts
table(combodat$sixpop)

#Coastal_Atlantic          Coastal_Gulf Intermediate_Atlantic     Intermediate_Gulf 
#46                    46                    44                    16 
#Offshore_Atlantic         Offshore_Gulf 
#33                    38 
sum(table(combodat$sixpop))
# 223

write.table(file="analysis/pop_structure/sixpop.clust", 
            data.frame(samp1 = combodat$indiv,
                       group = combodat$sixpop),
            sep="\t", quote=F, col.names=F, row.names=F)

