upSignal <-
function(AllSignals,AllStimuli,signal,Model,SigUp=vector()){
  # creates an array of nodes upstream the given signal (but downstream the previous signal)
  
  # 1. find reactions in which the investigated signal is the target node
  ix <- which(Model$interMat[rownames(Model$interMat)==signal] == 1)
  # 2. check if the investigated node is already in the array of upstream nodes
  #    and if it is a signal, if it is none of these, it adds it to the array
  if ((sum(SigUp == signal) == 0) && (sum(AllStimuli == signal) == 0)){
    SigUp[length(SigUp)+1] <- signal
  }
  if (length(ix) > 0) {
    # 3. finds which species are source nodes in those reactions (see point 1.)
    index <- unique(as.matrix(which(Model$interMat[,ix]==-1, arr.ind=TRUE))[,1])
    # 4. and keeps seraching/adding upstream nodes... 
    for (i in 1:length(index)){
      #...until the previous signal is reached
      if (sum(rownames(Model$interMat)[index[i]] == AllSignals) == 0){
        signal <- rownames(Model$interMat)[index[i]]
        SigUp <- upSignal(AllSignals = AllSignals, AllStimuli = AllStimuli, signal = signal, Model = Model, SigUp = SigUp)
      }
    }
  }
  return(SigUp)
}