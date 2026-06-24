library(strataG)
# original script from Ana Costa.

# Coastal-ATL vs Coastal-Gulf vs Offshore vs Intermediate (without hybrids)    

seqs <- read.fasta("./analysis/mt_pd_neid/AllAtl_CR_4pops_NoHybrids_n307.fasta")
df <- readGenData("./analysis/mt_pd_neid/AllAtl_CR_4pop_NoHybrids_n307.csv")

df.2 <- df[, c("id", "type", "id")]
colnames(df.2)[3] <- "dLoop"
g <- df2gtypes(df.2, ploidy = 1, sequences = seqs)
g

g.haps <- labelHaplotypes(g)

nuc.div <- nucleotideDivergence(g.haps, model = "TN93")
nuc.div

write.csv(nuc.div$between[,c(2,3,4)], file="analysis/mt_pd_neid/nei_da_output.csv", row.names=F, quote=F)
