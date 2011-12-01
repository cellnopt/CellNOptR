MapBTables2Model <-
function(BTable,Model,optimRes=NA,allInter=TRUE){
  # BTable are the Bolean Tables inferred from data using the function makeBTables
  # Model is the model optimised using CNO and than compressed
  
  ####for now links are added to the scaffol..THINK if to add them to the CNO optimized model
  #raw<-model2sif(Model=Model,optimRes=T1opt,writeSif=FALSE)
  #g<-sif2graph(raw)
  
  g<-sif2graph(model2sif(Model=Model, optimRes=optimRes))
  indexIntegr<-length(Model$reacID)
  
  ck=0
  # i signals
  # j stimuli
  # k inhibitors
  for (i in 1:length(BTable$namesSignals)){
    for (j in 1:dim(BTable$tables[[i]])[2]){
      
      print(paste("---------------------------------------------"))
      print(paste('stimulus: ', colnames(BTable$tables[[i]])[j], '     signal: ', BTable$namesSignals[i]))
      print(paste("Added links:"))
      
            
      # noNodes are the inhibitors that in the Btable of signal i, have a 1 in the column
      # corresponding to simulus j, meaning that the inhibitor has no effect
      noNodes <- rownames(as.matrix(which(BTable$tables[[i]][,j]==1)))
      noNodes <- setdiff(noNodes,BTable$namesSignals[i])
      
      ix <- as.vector(which(BTable$tables[[i]][,j]==2))
      if (length(ix) > 0){
        for (k in 1:length(ix)){
          # connection between the stimulus and the inhibitor
          cue <- colnames(BTable$tables[[i]])[j]
          signal <- rownames(BTable$tables[[i]])[ix[k]]
          
          ck <- SearchLinkGraph(node1 = cue,node2 = signal, graph=g, noNodes=noNodes)
          #ck <- SearchLink(node1 = cue,node2 = signal, noNodes=noNodes, Model = Model)
          if (ck == 0){
            # if there is no connection between the cue and the signal in the Model
            # all possible connections are added between nodes downstream the cue
            # (but upstream the following cue) and those upstream the signal
            # (but downstream the previous signal)
            if (allInter == TRUE){
              # 1. create a vector containing all nodes downstream the cue
              # (until reaching a cue)
              AllCues <- unique(c(colnames(BTable$tables[[i]]),rownames(BTable$tables[[i]])))
              CueDown <- downCueGraph(cue=cue, graph=g, stopNodes=AllCues)
              
              # 2. create a vector containing all nodes upstream the signal
              # (until reaching a signal)
              # AllSignals <- BTable$namesSignals
              AllSignals <-unique(c(BTable$namesSignals, noNodes))
              stopNodes<-unique(c(AllSignals,colnames(BTable$tables[[i]])))
              SigUp <- upSignalGraph(signal=signal, graph=g, stopNodes=stopNodes)

            }else{
              CueDown <- as.vector(cue)
              SigUp <- as.vector(signal)    
            }
            for (ix_tmp1 in 1:length(CueDown)){
              for (jx_pmp1 in 1:length(SigUp)){
                print(paste(CueDown[ix_tmp1], "->", SigUp[jx_tmp1],sep=""))
                Model <- AddLink(CueDown[ix_tmp1], SigUp[jx_tmp1], Model)
                #g<-sif2graph(model2sif(Model=Model))
              }
            }
          }
    
      
          # connection between the inhibitor and the signal
          cue <- rownames(BTable$tables[[i]])[ix[k]]
          signal <- BTable$namesSignals[i]
          #ck <- SearchLink(node1 = cue,node2 = signal, noNodes=noNodes, Model = Model)
          ck <- SearchLinkGraph(node1 = cue,node2 = signal, graph=g, noNodes=noNodes)
          if (ck == 0){
            # if there is no connection between the cue and the signal in the Model
            # all possible connections are added between nodes downstream the cue
            # (but upstream the following cue) and those upstream the signal
            # (but downstream the previous signal)
            if (allInter == TRUE){
              # 1. create a vector containing all nodes downstream the cue
              # (until reaching a cue)
              AllCues <- unique(c(colnames(BTable$tables[[i]]),rownames(BTable$tables[[i]])))
              CueDown <- downCueGraph(cue=cue, graph=g, stopNodes=AllCues)
              
              # 2. create a vector containing all nodes upstream the signal
              # (until reaching a signal)
              # AllSignals <- BTable$namesSignals
              AllSignals <-unique(c(BTable$namesSignals, noNodes))
              stopNodes<-unique(c(AllSignals,colnames(BTable$tables[[i]])))
              SigUp <- upSignalGraph(signal=signal, graph=g, stopNodes=stopNodes)
            }else{
              CueDown <- as.vector(cue)
              SigUp <- as.vector(signal)    
            }
            for (ix_tmp1 in 1:length(CueDown)){
              for (jx_tmp1 in 1:length(SigUp)){
                print(paste(CueDown[ix_tmp1], "->", SigUp[jx_tmp1],sep=""))
                Model <- AddLink(CueDown[ix_tmp1], SigUp[jx_tmp1], Model)
                #g<-sif2graph(model2sif(Model=Model))
              }
            }      
          }
        }
      }
              
      # additional connection between the stimulus and the signal
	    if (sum(BTable$tables[[i]][,j]) == dim(BTable$tables[[i]])[2]){
	      cue <- colnames(BTable$tables[[i]])[j]
        signal <- BTable$namesSignals[i]
        #ck <- SearchLink(node1 = cue,node2 = signal, noNodes=noNodes, Model = Model)
        ck <- SearchLinkGraph(node1 = cue,node2 = signal, graph=g, noNodes=noNodes)
        if (ck == 0){
          # if there is no connection between the cue and the signal in the Model
          # all possible connections are added between nodes downstream the cue
          # (but upstream the following cue) and those upstream the signal
          # (but downstream the previous signal)
          if (allInter == TRUE){
            # 1. create a vector containing all nodes downstream the cue
            # (until reaching a cue)
            AllCues <- unique(c(colnames(BTable$tables[[i]]),rownames(BTable$tables[[i]])))
            CueDown <- downCueGraph(cue=cue, graph=g, stopNodes=AllCues)
            
            # 2. create a vector containing all nodes upstream the signal
            # (until reaching a signal)
            # AllSignals <- BTable$namesSignals
            AllSignals <-unique(c(BTable$namesSignals, noNodes))
            stopNodes<-unique(c(AllSignals,colnames(BTable$tables[[i]])))
            SigUp <- upSignalGraph(signal=signal, graph=g, stopNodes=stopNodes)            
          }else{
            CueDown <- as.vector(cue)
            SigUp <- as.vector(signal)    
          }
	        for (ix_tmp1 in 1:length(CueDown)){
            for (jx_tmp1 in 1:length(SigUp)){
              print(paste(CueDown[ix_tmp1], "->", SigUp[jx_tmp1],sep=""))
              Model <- AddLink(CueDown[ix_tmp1], SigUp[jx_tmp1], Model)
              #g<-sif2graph(model2sif(Model=Model))
            }
          }
        }
	    }
    print(paste("---------------------------------------------"))
    print(paste(""))
    }
  }
  Model$indexIntegr<-seq(from=indexIntegr+1, to=length(Model$reacID))
  return(Model)
}