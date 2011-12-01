BTables2Sif <-
function(BTable, writeSif=TRUE){
  #BTables are boolean tables inferred from data
  
  # i signals
  # j stimuli
  # k inhibitors
  sifFile <- NULL
  for (i in 1:length(BTable$namesSignals)){
    for (j in 1:dim(BTable$tables[[i]])[2]){
      if (sum(BTable$tables[[i]][,j]) == dim(BTable$tables[[i]])[1]){
        link <- cbind(colnames(BTable$tables[[i]])[j], 1, BTable$namesSignals[[i]])
        sifFile = rbind(sifFile, link)
      }else{
        for (k in 1:dim(BTable$tables[[i]])[1]){
          if (BTable$tables[[i]][k,j] == 2){
            link1 <- cbind(colnames(BTable$tables[[i]])[j], 1, rownames(BTable$tables[[i]])[k])
            sifFile = rbind(sifFile, link1)
            link2 <- cbind(rownames(BTable$tables[[i]])[k], 1, BTable$namesSignals[[i]])
            sifFile = rbind(sifFile, link2)
          }
        }
      }
      
    }
  }
  
  sifFile.new <- NULL
  sifFile.new <- rbind(sifFile.new, sifFile[1,])
  for (i in 2:dim(sifFile)[1]){
    ck <- 0
    for (j in 1:i-1){
      if (sum(sifFile[i,] == sifFile[j,]) == 3){
        ck <- 1
      }
    }
    if (ck == 0){
      sifFile.new <- rbind(sifFile.new, sifFile[i,])
    }
  }
  sifFile <- sifFile.new
  
  if (writeSif == TRUE){
	write.table(sifFile, file="CauseEffectNet.sif", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
	
  return(sifFile)
}