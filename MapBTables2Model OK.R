MapBTables2Model <-
function(BTable,Model,allInter=TRUE){
  # BTable are the Bolean Tables inferred from data using the function makeBTables
  # Model is the model optimised using CNO and than compressed

  # i signals
  # j stimuli
  # k inhibitors
  for (i in 1:length(BTable$namesSignals)){
    for (j in 1:dim(BTable$tables[[i]])[2]){

      # noNodes are the inhibitors that in the Btable of signal i, have a 1 in the column
      # corresponding to simulus j, meaning that the inhibitor has no effect
      noNodes <- rownames(as.matrix(which(BTable$tables[[i]][,j]==1)))
      
      ix <- as.vector(which(BTable$tables[[i]][,j]==2))
      if (length(ix) > 0){
        for (k in 1:length(ix)){
          # connection between the stimulus and the inhibitor
          cue <- colnames(BTable$tables[[i]])[j]
          signal <- rownames(BTable$tables[[i]])[ix[k]]
          ck <- SearchLink(node1 = cue,node2 = signal, noNodes=noNodes, Model = Model)
          if (ck == 0){
            # if there is no connection between the cue and the signal in the Model
            # all possible connections are added between nodes downstream the cue
            # (but upstream the following cue) and those upstream the signal
            # (but downstream the previous signal)
            if (allInter == TRUE){
              # 1. create a vector containing all nodes downstream the cue
              # (until reaching a cue)
              AllCues <- unique(c(colnames(BTable$tables[[i]]),rownames(BTable$tables[[i]])))
              CueDown <- downCue(AllCues = AllCues, cue = cue, Model = Model)
              # 2. create a vector containing all nodes upstream the signal
              # (until reaching a signal)
              # AllSignals <- BTable$namesSignals
              AllSignals <-unique(c(BTable$namesSignals, noNodes))
              SigUp <- upSignal(AllSignals = AllSignals, AllStimuli = colnames(BTable$tables[[i]]), signal = signal, Model = Model)
            }else{
              CueDown <- as.vector(cue)
              SigUp <- as.vector(signal)    
            }
            Model <- AddLinks(CueDown, SigUp, Model)
          }
    
      
          # connection between the inhibitor and the signal
          cue <- rownames(BTable$tables[[i]])[ix[k]]
          signal <- BTable$namesSignals[i]
          ck <- SearchLink(node1 = cue,node2 = signal, noNodes=noNodes, Model = Model)
          if (ck == 0){
            # if there is no connection between the cue and the signal in the Model
            # all possible connections are added between nodes downstream the cue
            # (but upstream the following cue) and those upstream the signal
            # (but downstream the previous signal)
            if (allInter == TRUE){
              # 1. create a vector containing all nodes downstream the cue
              # (until reaching a cue)
              AllCues <- unique(c(colnames(BTable$tables[[i]]),rownames(BTable$tables[[i]])))
              CueDown <- downCue(AllCues = AllCues, cue = cue, Model = Model)
              # 2. create a vector containing all nodes upstream the signal
              # (until reaching a signal)
              # AllSignals <- BTable$namesSignals
              AllSignals <-unique(c(BTable$namesSignals, noNodes))
              SigUp <- upSignal(AllSignals = AllSignals, AllStimuli = colnames(BTable$tables[[i]]), signal = signal, Model = Model)
              Model <- AddLinks(CueDown, SigUp, Model)
            }else{
              CueDown <- as.vector(cue)
              SigUp <- as.vector(signal)    
            }
            Model <- AddLinks(CueDown, SigUp, Model)       
          }
        }
      }
              
      # additional connection between the stimulus and the signal
	  if (sum(BTable$tables[[i]][,j])>0){
	    cue <- colnames(BTable$tables[[i]])[j]
        signal <- BTable$namesSignals[i]
        ck <- SearchLink(node1 = cue,node2 = signal, noNodes=noNodes, Model = Model)
        if (ck == 0){
          # if there is no connection between the cue and the signal in the Model
          # all possible connections are added between nodes downstream the cue
          # (but upstream the following cue) and those upstream the signal
          # (but downstream the previous signal)
          if (allInter == TRUE){
            # 1. create a vector containing all nodes downstream the cue
            # (until reaching a cue)
            AllCues <- unique(c(colnames(BTable$tables[[i]]),rownames(BTable$tables[[i]])))
            CueDown <- downCue(AllCues = AllCues, cue = cue, Model = Model)
            # 2. create a vector containing all nodes upstream the signal
            # (until reaching a signal)
            # AllSignals <- BTable$namesSignals
            AllSignals <-unique(c(BTable$namesSignals, noNodes))
            SigUp <- upSignal(AllSignals = AllSignals, AllStimuli = colnames(BTable$tables[[i]]), signal = signal, Model = Model)
          }else{
            CueDown <- as.vector(cue)
            SigUp <- as.vector(signal)    
          }
		  Model <- AddLinks(CueDown, SigUp, Model)
        }
	  }

    }
  }
  return(Model)
}