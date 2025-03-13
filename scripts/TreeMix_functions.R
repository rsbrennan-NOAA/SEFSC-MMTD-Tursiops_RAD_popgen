######################### R functions for the Treemix Pipeline ############################
###########################################################################################
# 
# adapted from Carolin Dahms from https://github.com/carolindahms/TreeMix/blob/main/TreeMix_functions.R
# CD originally adapted Zecca et al. (2019)


pckgs <- c("ape","phytools");
lapply(pckgs, library, character.only = TRUE)

maxLL<-function (input_stem, nt, uel=FALSE){
  lliknames<-c(paste(input_stem,1:nt,'.llik',sep='')) 
  d<-unname(sapply(lliknames, function(x) rbind(read.table(x,header = F, sep = ":")[2,2])))
  d<-as.data.frame(cbind(1:nt,d)); names(d)[1:2]<- c('treeN','LL')
  d<-d[order(-d$LL,d$treeN),]
  bestLL<-d[,"LL"] %in% d[1,"LL"]
  bestN<- d$treeN[bestLL]	
  bestlikelihood<-unique(max(d$LL))
  bestnames<-c(paste0(input_stem,bestN,'.treeout.gz'))
  BKLTs<-list()
  class(BKLTs)<-"multiPhylo"
  for (t in 1:length(bestnames)){
    tree<-read.tree(file = bestnames[t],keep.multi = TRUE)				
    BKLTs[[t]]<-tree[[1]]												
    BKLTs[[t]][5]<- bestN[t]}
  uniqueT<- unique.multiPhylo(BKLTs, use.edge.length = uel, use.tip.label = TRUE)               
  if(length(BKLTs)==1){
    keep<-sapply(uniqueT,'[[',5)
    UniNames<- c(paste0(input_stem,keep,'.treeout.gz'))
    write.table(d,file='TreeLLs.txt',sep=' ',dec = ".",eol = "\n", row.names = F, col.names =T ,quote = F, na = "NA")
    cat('\nNumber of input trees: ', nt,
        '\nNumber of best-known ML trees: ' ,length(BKLTs),
        '\nML trees: ', bestnames,
        '\nBest likelihood of ML tree(s): ',bestlikelihood,
        '\nNumber of unique ML tree(s): ',length(keep),
        '\nRun of ML tree(s) with unique topology: ',UniNames)
  }else{
    keep<-sapply(uniqueT,'[[',5)
    UniNames<- c(paste0(input_stem,keep,'.treeout.gz'))
    write.table(d,file='TreeLLs.txt',sep=' ',dec = ".",eol = "\n", row.names = F, col.names =T ,quote = F, na = "NA")
    cat('\nNumber of input trees: ', nt,
        '\nNumber of best-known ML trees: ' ,length(BKLTs),
        '\nML trees: ', bestnames,
        '\nBest likelihood of ML tree(s): ',bestlikelihood,
        '\nNumber of unique ML tree(s): ',length(uniqueT),
        '\nRun of ML tree(s) with unique topology: ',UniNames)
  }
  
  # Load all trees to compare topologies
  all_trees <- list()
  class(all_trees) <- "multiPhylo"
  for (t in 1:nt){
    tree_file <- paste0(input_stem, t, '.treeout.gz')
    if(file.exists(tree_file)) {
      tree <- read.tree(file = tree_file, keep.multi = TRUE)
      all_trees[[t]] <- tree[[1]]
      all_trees[[t]][5] <- t  # Store the original tree number
    }
  }
  
  # Find unique topologies among all trees
  uniqueT_all <- unique.multiPhylo(all_trees, use.edge.length = uel, use.tip.label = TRUE)
  
  # Create a dataframe to store results
  result_df <- data.frame(
    UniqueTopology = integer(length(uniqueT_all)),
    Count = integer(length(uniqueT_all)),
    FirstTree = integer(length(uniqueT_all)),
    AllTrees = character(length(uniqueT_all)),
    stringsAsFactors = FALSE
  )
  
  # For each unique topology, find all trees with that topology
  for (i in 1:length(uniqueT_all)) {
    trees_in_group <- c()
    
    for (j in 1:length(all_trees)) {
      if (!is.null(all_trees[[j]]) && 
          identical(all.equal(uniqueT_all[[i]], all_trees[[j]], use.edge.length = uel, use.tip.label = TRUE), TRUE)) {
        trees_in_group <- c(trees_in_group, j)
      }
    }
    
    result_df$UniqueTopology[i] <- i
    result_df$Count[i] <- length(trees_in_group)
    result_df$FirstTree[i] <- min(trees_in_group)
    result_df$AllTrees[i] <- paste(trees_in_group, collapse = ",")
  }
  
  # Identify which unique topology corresponds to the best likelihood trees
  for (i in 1:length(uniqueT_all)) {
    for (j in 1:length(bestN)) {
      if (identical(all.equal(uniqueT_all[[i]], all_trees[[bestN[j]]], use.edge.length = uel, use.tip.label = TRUE), TRUE)) {
        result_df$IsBestML[i] <- TRUE
        break
      } else {
        result_df$IsBestML[i] <- FALSE
      }
    }
  }
  
  # Sort by count in descending order
  result_df <- result_df[order(-result_df$Count), ]
  
  # Find which unique tree corresponds to the best likelihood tree
  best_tree_group <- NA
  for (i in 1:length(uniqueT_all)) {
    for (j in 1:length(bestN)) {
      if (identical(all.equal(uniqueT_all[[i]], all_trees[[bestN[j]]], use.edge.length = uel, use.tip.label = TRUE), TRUE)) {
        best_tree_group <- i
        break
      }
    }
    if (!is.na(best_tree_group)) break
  }
  # Get the top 10 (or fewer if nt < 10) trees by likelihood
  top_trees <- d[1:min(10, nrow(d)), ]
  
  # Create the final result list containing both dataframes
  final_result <- list(
    unique_topologies = result_df,
    top_likelihood_trees = top_trees
  )
  # Print summary information
  cat('\nNumber of input trees: ', nt)
  cat('\nNumber of unique tree topologies: ', length(uniqueT_all))
  cat('\nUnique tree topology count summary:\n')
  #print(result_df)
  cat('\n Best', bestN[1], 'tree belongs to group: ', best_tree_group, '\n')
  
  # Return the dataframe
  return(final_result)
}

cfTrees<-function (input_stem,nt,p=1, uel=FALSE, m='PH85'){
  lliknames<-c(paste(input_stem,1:nt,'.llik',sep='')) 
  d<-unname(sapply(lliknames, function(x) rbind(read.table(x,header = F, sep = ":")[2,2])))
  d<-as.data.frame(cbind(1:nt,d)); names(d)[1:2]<- c('treeN','LL')
  d<-d[order(-d$LL,d$treeN),]
  bestLL<-d[,"LL"] %in% d[1,"LL"]
  bestN<- d$treeN[bestLL]	
  bestnames<-c(paste0(input_stem,bestN,'.treeout.gz'))
  BKLTs<-list()
  class(BKLTs)<-"multiPhylo"
  for (t in 1:length(bestnames)){
    tree<-read.tree(file = bestnames[t],keep.multi = TRUE)				
    BKLTs[[t]]<-tree[[1]]												
    BKLTs[[t]][5]<- bestN[t]}
  if(length(BKLTs)==1){
    uniqueT<-BKLTs
  }else{
    uniqueT<- unique(BKLTs, use.edge.length = uel ,use.tip.label = TRUE)}	
  keep<-sapply(uniqueT,'[[',5)
  jump<-bestN[!bestN%in%keep]	
  UniNames<- c(paste0(input_stem,keep,'.treeout.gz'))
  if(length(uniqueT)==1){
    tab<- matrix(0); 
    colnames(tab)<-c(paste0('tree',keep)); rownames(tab)<-c(paste0('tree',keep))
  }else{
    tab<-dist.topo(unroot(uniqueT),method = m);
    attributes(tab)$Diag<-TRUE;attributes(tab)$Labels<-paste0('tree',as.character(keep))}
  cons<-consensus(uniqueT, p=p, check.labels=T)
  write.table(d,file='TreeLLs.txt',sep=' ',dec = ".",eol = "\n", row.names = F, col.names =T ,quote = F, na = "NA")
  write.tree(BKLTs,file = "AllBestTrees.newick",tree.names =FALSE,append =FALSE)
  write.tree(uniqueT,file='UniqueTop.newick')
  write(UniNames,sep='',ncolumns=1, file='UniqueNames.txt');	
  write.table(as.matrix(tab),file='PairwDist.txt',sep='\t',row.names = TRUE, quote = FALSE, col.names=NA)
  write.tree(cons,file='Consensus.newick')
  plot(cons,label.offset=0.2)	
  cat('\nInput trees: ',nt,
      '\nBest-known ML trees: ' ,length(BKLTs),
      '\n\nNumber of duplicates in the Best-known ML trees: ',length(jump),
      '\nDuplicate trees #:\n ')
  cat(jump, sep = " ", fill =110);
  cat('\nNumber of unique ML tree(s): ',length(uniqueT),
      '\nUnique ML tree(s) #: ',UniNames,
      '\n\nPairwise topological distance ( method =',m,') between unique trees:\n');
  print(tab)
  cat('\n\nJob done.\n\n')
}

ImportData <- function(skipL){
  filevector <- list.files(path = ".", pattern = "\\.treeout.gz$",full.names = TRUE)  
  mydata <- do.call('rbind',lapply(filevector, function(x) read.table(x,skip =skipL, stringsAsFactors=FALSE))) 
  mydata <- lapply(mydata[5:6], function (x) paste0('(',x,');'))
  From <- paste0("(",sapply(read.newick(text=mydata$V5),function(x) paste(sort(x$tip.label),sep='', collapse=',')),");"); 
  To <- paste0("(",sapply(read.newick(text=mydata$V6), function(x) paste(sort(x$tip.label),sep='', collapse=',')),");"); 
  FromTo <- list(From=From, To=To,filevector=filevector, mydata=mydata)
}


f.MSets <- function(m_from, m_to,inputfile){
  if(m_from!="all" && m_to!="all"){
    MSets <- vector(mode = 'list', 2);names(MSets) <- c("From", "To")
    MSets$From <- m_from 
    MSets$To <- m_to   
  } else{
    MSets <- scan(file= inputfile, skip=1, what =list('','')); names(MSets)<-c('From','To')}
  cat("\n\n")	
  MSets
}


f.fileout <- function(x,fileout){
  write.table(x,file=fileout,quote = FALSE, sep = "   ",row.names = FALSE,col.names = TRUE)
}


GetPairsOfSets <- function(outputfile='PairsOfSets.txt',skipL=1){
  FromTo <- ImportData (skipL)
  Migrations <- unique(data.frame(FromTo[c("From","To")], stringsAsFactors=FALSE)) 
  f.fileout(Migrations,outputfile)
}


GetMigrSupp<-function(m_from="all", m_to="all",inputfile='PairsOfSets.txt', outputfile='MigrSupp.txt', skipL=1){
  FromTo <- ImportData(skipL)
  MSets <- f.MSets(m_from, m_to,inputfile)
  boots<-vector(mode = 'list',1); names(boots)<- 'MS'
  for (i in seq_along(MSets$From)){
    FT <- sapply(as.list(FromTo$From), function(x) identical(MSets$From[[i]],x))
    TT <- sapply(as.list(FromTo$To), function(x) identical(MSets$To[[i]],x))
    boots$MS[[i]]<-round((sum(FT & TT, na.rm=TRUE)/length(FromTo$filevector))*100)} 
  f.fileout(c(MSets,boots),outputfile)
  print(boots)
}

GetExtMigrSupp_ed<-function(skipL, min_n, nmigr, m_from, m_to, fixed, inputfile='PairsOfSets.txt'){
  FromTo <- ImportData (skipL)
  MSets <- f.MSets(m_from, m_to,inputfile)
  MSets2 <- lapply(MSets, function(x) strsplit(gsub("[();]", "", x),","))
  From <- lapply(read.newick(text=FromTo$mydata$V5),function(x) sort(x$tip.label))
  To <- lapply(read.newick(text=FromTo$mydata$V6),function(x) sort(x$tip.label))
  boots<-vector(mode = 'list', 1); names(boots)<- 'ExtMS'
  counts<-vector(mode = 'list', 2); names(counts)<- c('TotalCount','CorrectedCount')
  FT <- vector(mode = 'list', length(MSets$From));TT <- vector(mode = 'list', length(MSets$From)); 
  tot<-vector(mode = 'integer', length(FromTo$filevector))
  for (i in seq_along(MSets$From)){
    if(fixed[i] == "To"){
      FT[[i]] <- sapply(From, function(x) ifelse(length(intersect(MSets2$From[[i]],x))>= min_n[i] && length(x)==sum(x %in% MSets2$From[[i]]), TRUE,FALSE))  # Only  proper subsets allowed
      TT[[i]] <- sapply(To, function(x) identical(MSets2$To[[i]],x))
    }else if (fixed[i] == "From"){
      TT[[i]] <- sapply(To, function(x) ifelse(length(intersect(MSets2$To[[i]],x))>= min_n[i] && length(x)==sum(x %in% MSets2$To[[i]]), TRUE,FALSE))  # Only  proper subsets allowed
      FT[[i]] <- sapply(From, function(x) identical(MSets2$From[[i]],x))
    }else{
      stop('Please, provide a valid value for the argument <fixed>. Valid options are: "From", "To"  or a character vector that contains any sensible combination of "From" and "To".', call. = FALSE)}
    fctr<-rep(1:length(FromTo$filevector), each=nmigr)			
    splFT<-split(FT[[i]], fctr)									
    splTT<-split(TT[[i]], fctr)									
    for(k in seq_along(FromTo$filevector)){ 					# Correct for multiple matches
      tot[k]<-sum(splFT[[k]] & splTT[[k]], na.rm=TRUE)}	
    counts$TotalCount[[i]]<-sum(tot)							
    tot[tot>1]<-1												
    counts$CorrectedCount[[i]]<-sum(tot)
    boots$ExtMS[[i]] <-round((counts$CorrectedCount[[i]]/length(FromTo$filevector))*100)}
  n_migr<-paste0('migr#',as.character(seq_along(MSets$From)))
  d<-as.character(boots$ExtMS)
}

GetMS_MSe<-function(skipL=1, min_n, nmigr, m_from=all, m_to=all, fixed, inputfile='MigrSupp.txt', outputfile='MS_MSE.txt'){
  MS_df<-read.table(inputfile,header=T) #read in data table created with GetMigrSupp
  MSval<-as.numeric(MS_df[2,3:ncol(MS_df)])
  MS_df<-MS_df[,-3:-ncol(MS_df)]
  MS_df$MSval<-MSval
  MSE<-list()
  for (i in 1:nrow(MS_df)){
    MS_df$MSE[[i]]<-unlist(GetExtMigrSupp_ed(skipL, min_n, nmigr, m_from=MS_df$From[i], m_to=MS_df$To[i], fixed, inputfile='PairsOfSets.txt'))
  }
  MS_df$MSE <- vapply(MS_df$MSE, paste, collapse = ", ", character(1L))
  f.fileout(MS_df,outputfile) 
}

GetMigrStats<-function(input_stem, nt, skipL=1, outputfile='MigrationStats.txt'){
  treenames<-c(paste(input_stem,1:nt,'.treeout.gz',sep='')) 
  myfiles = lapply(treenames, read.table, header=F, sep="",fill=TRUE, skip=skipL)
  d<-rbindlist(myfiles);names(d)[1:6]<- c('W','Wj','SEj','pval','From','To')
  d$From<-gsub("[^ a-zA-Z',]", "", d$From)
  d$To<-gsub("[^ a-zA-Z',]", "", d$To)
  d$pval<-gsub("<", "", d$pval)
  d$pval<-as.numeric(d$pval)
  
  d$From<- sapply(d$From, function(x) paste(sort(unlist(strsplit(x, ","))),
                                            collapse = ","))
  d$To<- sapply(d$To, function(x) paste(sort(unlist(strsplit(x, ","))),
                                        collapse = ","))
  d$pair<- with(d, paste(From,To,sep="->"))
  setDT(d)[, meanW := mean(W), by = pair]
  setDT(d)[, meanWj := mean(Wj), by = pair]
  setDT(d)[, meanSEj := mean(SEj), by = pair]
  setDT(d)[, maxpval := max(pval), by = pair]
  d$N<-nt
  Results<-d[!duplicated(d[,c('pair')]),];Results<-Results[,-1:-6] 
  f.fileout(Results, outputfile)
  
}