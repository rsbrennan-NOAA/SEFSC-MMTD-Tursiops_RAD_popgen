library(rfPermute)

# original script from Ana Costa.

createSiteFactors <- function(df) {
  as.data.frame(
    sapply(colnames(df), function(x) {
      if(x == "type") {
        factor(df[[x]])
      } else {
        factor(df[[x]], levels = c(LETTERS, "-"))
      }
    }, simplify = FALSE)
  )
}

T.hyb <- read.table(
  "./analysis/mt_pd_neid/AllAtl_CR_3pops_NoHybrids_n259RF.txt", 
  header = TRUE, 
  colClasses = "character",
  stringsAsFactors = FALSE
)


T.hyb <- createSiteFactors(T.hyb)
set.seed(20)

RF_Tt <- T.hyb[, -1]

ntree <- 10000
num.cores <- 8
nrep <- 1000

sampsize.Tt <- balancedSampsize(RF_Tt$type)

rf.Tt <- rfPermute(
  type ~ .,
  data = RF_Tt,
  sampsize = sampsize.Tt,
  ntree = ntree,
  replace = FALSE,
  proximity = TRUE,
  nrep = nrep,
  num.cores = num.cores
)

rf.Tt.3 <- rf.Tt


# write summary table:
cf_out <- as.data.frame(confusionMatrix(rf.Tt.3))
pc_out <- pctCorrect(rf.Tt, pct = c(0.8, 0.95))
cf_out$pct.correct_0.8 <- pc_out$pct.correct_0.8
cf_out$pct.correct_0.95 <- pc_out$pct.correct_0.95

write.csv(cf_out, file="analysis/mt_pd_neid/PD_noInt.csv", row.names=F, quote=F)

#votes

rf.Tt$rf$votes
summary(rf.Tt)
# write.table(rf.Tt$rf$votes, file = "RF_votes.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Plotting

plotTrace(rf.Tt)
plotInbag(rf.Tt, replace = FALSE, sampsize = sampsize.Tt)
plotVotes(rf.Tt)
plotProximity(rf.Tt)
plotImportance(rf.Tt)




#----------------------------------------------------------
## full data:

T.hyb <- read.table(
  "./analysis/mt_pd_neid/AllAtl_CR_3pops_NoHybrids_n307RF.txt", 
  header = TRUE, 
  colClasses = "character",
  stringsAsFactors = FALSE
)

table(T.hyb$type)

pops <- read.csv("metadata_FINAL.csv", header=T)
sum(table(pops$dapc_population))
table(pops$dapc_population)

T.hyb <- createSiteFactors(T.hyb)
set.seed(20)

RF_Tt <- T.hyb[, -1]

ntree <- 10000
num.cores <- 8
nrep <- 1000

sampsize.Tt <- balancedSampsize(RF_Tt$type)

rf.Tt <- rfPermute(
  type ~ .,
  data = RF_Tt,
  sampsize = sampsize.Tt,
  ntree = ntree,
  replace = FALSE,
  proximity = TRUE,
  nrep = nrep,
  num.cores = num.cores
)

rf.Tt.4 <- rf.Tt
rf.Tt.4


#votes

rf.Tt$rf$votes
summary(rf.Tt.4)

plotConfMat(rf.Tt.4)

# write summary table:
cf_out <- as.data.frame(confusionMatrix(rf.Tt.4))
pc_out <- pctCorrect(rf.Tt.4, pct = c(0.8, 0.95))
cf_out$pct.correct_0.8 <- pc_out$pct.correct_0.8
cf_out$pct.correct_0.95 <- pc_out$pct.correct_0.95

write.csv(cf_out, file="analysis/mt_pd_neid/PD_all.csv", row.names=F, quote=F)

# write.table(rf.Tt$rf$votes, file = "RF_votes.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Plotting
plotInbag(rf.Tt, replace = FALSE, sampsize = sampsize.Tt)

plotTrace(rf.Tt)
plotInbag(rf.Tt, replace = FALSE, sampsize = sampsize.Tt)
plotVotes(rf.Tt)
plotProximity(rf.Tt)
plotImportance(rf.Tt)
