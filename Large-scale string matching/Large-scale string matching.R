s <- Sys.time()

library(stringr)
library(doSNOW)
library(progress)
library(tcltk)

PMC <- read.csv(file = "Abstract-WOS.csv")
PMC <- PMC[, c(2,13)]
colnames(PMC) <- c("PMID", "FullText")
PMC$PMID <- c(1:nrow(PMC))
PMC <- PMC[which(!is.na(PMC$FullText)), ]
# PMC <- PMC[1:100,]
PMC$FullText <- tolower(PMC$FullText)

chem_Synonym <- read.table("CID-Synonym-filtered.gz", header=F, stringsAsFactors = F,sep="\t",quote = "")
chem_Synonym$V2 <- tolower(chem_Synonym$V2)

remIndex <- grep("[0-9]{5}",chem_Synonym$V2)
chem_Synonym <- chem_Synonym[-remIndex,]

synonym_exclusion_list <- readLines("Compound synonym exclusion list.txt")
chem_Synonym <- chem_Synonym[which(!chem_Synonym$V2 %in% synonym_exclusion_list),]

colnames(chem_Synonym) <- c("CID", "Synonym")
chem_Synonym <- chem_Synonym[which(nchar(chem_Synonym$Synonym)>3),]

qfunc <- function(Synonym){
  chem.name <- gsub("\\\\", "", Synonym)
  chem.name <- gsub("\\[", "\\\\[", chem.name)
  chem.name <- gsub("\\]", "\\\\]", chem.name)
  chem.name <- gsub("\\{", "\\\\{", chem.name)
  chem.name <- gsub("\\}", "\\\\}", chem.name)
  chem.name <- gsub("\\(", "\\\\(", chem.name)
  chem.name <- gsub("\\)", "\\\\)", chem.name)
  chem.name <- gsub("\\++", "\\\\++", chem.name)
  chem.name <- gsub("\\*", "\\\\*", chem.name)
  chem.name <- gsub("\\?", "\\\\?", chem.name)
  chem.name <- gsub("\\|", "\\\\|", chem.name)
  return(chem.name)
}

chem_Synonym <- data.frame(CID = chem_Synonym$CID, Synonym = chem_Synonym$Synonym, Synonym_Q = qfunc(chem_Synonym$Synonym))

chem_Synonym <- data.frame(CID = chem_Synonym$CID, Synonym = chem_Synonym$Synonym, Synonym_Q = chem_Synonym$Synonym_Q, MergedID = paste(chem_Synonym$CID,chem_Synonym$Synonym,chem_Synonym$Synonym_Q, sep = ".....Delimiter....."))


################################################
cl <- makeSOCKcluster(4)
registerDoSNOW(cl)

# progress bar ------------------------------------------------------------
iterations <- nrow(PMC)
pb <- progress_bar$new(
  format = ":letter [:bar] :elapsed | Remaining time: :eta <br>",
  total = iterations,
  width = 120)
# allowing progress bar to be used in foreach -----------------------------
progress <- function(n){
  pb$tick(tokens = list(letter = "Progress of CID Match"))
}
opts <- list(progress = progress)

matchfunc <- function(i){
  posi <- stringr::str_locate_all(PMC$FullText[i],chem_Synonym$Synonym_Q)
  names(posi) <- chem_Synonym$MergedID
  posi <- posi[lapply(posi, length)>0]
  if (length(posi) == 0){
    matched_res <- matrix(data = NA,nrow = 0, ncol = 7)
  }
  else {
    matched_res <- matrix(data = NA,nrow = 0, ncol = 7)
    for (j in c(1:length(posi))){
      names.j <- unlist(strsplit(names(posi)[j],split = ".....Delimiter....."))
      names.j <- as.data.frame(matrix(data = rep(names.j,each = nrow(posi[[j]])),nrow = nrow(posi[[j]]),ncol = 3))
      matched_res <- rbind(matched_res, as.matrix(cbind(PMC$PMID[i], "...", names.j, posi[[j]])))
    }
    matched_res[,2] <- PMC$FullText[i]
  }
  return(matched_res)
}

abstractRes <- foreach(i=1:nrow(PMC), .options.snow=opts, .combine='rbind') %dopar% matchfunc(i)
abstractRes <- as.data.frame(abstractRes)
colnames(abstractRes) <- c("PMID", "PMC", "CID", "Synonym", "Synonym_Q", "Start", "End")
PMCfulltext.index <- which(abstractRes$PMC != "...")
PMCfulltext.index <- c(PMCfulltext.index,nrow(abstractRes)+1)
abstractRes$PMCIndex <- NA
pb <- tkProgressBar("Synonym","rate of progress %", 0, 100)
for (i in c(1:nrow(abstractRes))){
  info<- sprintf("rate of progress %d%%", round(i*100/nrow(abstractRes)))
  setTkProgressBar(pb, i*100/nrow(abstractRes), sprintf("Synonym", info),info)
  abstractRes$PMCIndex[i] <- PMCfulltext.index[which((PMCfulltext.index-i) > 0)[1]-1]
}
close(pb)

closeAllConnections()

###############################################################
checkFunc <- function(Start, End, PMC){
  Nearby <- substr(PMC,as.numeric(Start)-9,as.numeric(End)+9)
  # print(paste(check.start.value, "&", check.end.value))
  return(Nearby)
}

leftValueFunc <- function(Start, PMC){
  check.start.index <- as.numeric(Start)-1
  check.start.value <- substr(PMC,check.start.index,check.start.index)
  return(check.start.value)
}

rightValueFunc <- function(End, PMC){
  check.end.index <- as.numeric(End)+1
  check.end.value <- substr(PMC,check.end.index,check.end.index)
  return(check.end.value)
}

PMC_2 <- abstractRes$PMC[abstractRes$PMCIndex]
abstractRes <- data.frame(abstractRes, 
                          Nearby = checkFunc(abstractRes$Start, abstractRes$End, PMC_2),
                          leftValue = leftValueFunc(abstractRes$Start, PMC_2),
                          rightValue = rightValueFunc(abstractRes$End, PMC_2))

abstractRes <- abstractRes[which((abstractRes$leftValue == " " | abstractRes$leftValue == "," | abstractRes$leftValue == ".") &
                                   (abstractRes$rightValue == " " | abstractRes$rightValue == "," | abstractRes$rightValue == ".")),]

e <- Sys.time()
print(e-s)

write.csv(abstractRes, file = "Abstract-Large-scale match.csv")
