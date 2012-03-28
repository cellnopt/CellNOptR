MapCEN2Model <-
function(CEN,Model,CNOlist,allInter=TRUE){
  g<-sif2graph(model2sif(Model=Model))
  indexIntegr<-length(Model$reacID)
  
  namesStimuli<-CNOlist$namesStimuli
  namesInhibitors<-CNOlist$namesInhibitors
  namesSignals<-CNOlist$namesSignals
  #Others are the white nodes in the network, those that are not stimulated, nor inhibited, nor measures
  namesOthers<-setdiff(setdiff(setdiff(Model$namesSpecies, namesSignals), namesInhibitors), namesStimuli)
  
  for (i in 1:dim(CEN)[1]){
    noNodes <- namesOthers
    node1<-CEN[i,1]
    node2<-CEN[i,3]
    
    print(paste("---------------------------------------------"))
    print(paste('stimulus: ', node1, '     signal: ', node2))
    print(paste("Added links:"))
    
    
    ck <- SearchLinkGraph(node1 = node1,node2 = node2, graph=g, noNodes=noNodes)
    if (ck==0){
      if (allInter==TRUE){
        # 1. create a vector containing all nodes downstream the cue
        # (until reaching a cue or a noNode)
        AllCues <- union(namesStimuli,namesInhibitors)
        CueDown <- downCueGraph(cue=node1, graph=g, stopNodes=union(AllCues,noNodes))
        # 2. create a vector containing all nodes upstream the signal
        # (until reaching a signal, a Stimulus, a noNode)
        # AllSignals <- BTable$namesSignals
        stopNodes <-unique(union(namesSignals, noNodes))
        stopNodes<-unique(union(stopNodes,namesStimuli))
        SigUp <- upSignalGraph(signal=node2, graph=g, stopNodes=stopNodes)
      }else{
        CueDown <- as.vector(node1)
        SigUp <- as.vector(node2) 
      }
      for (ix_tmp1 in 1:length(CueDown)){
        for (jx_tmp1 in 1:length(SigUp)){
          print(paste(CueDown[ix_tmp1], "->", SigUp[jx_tmp1],sep=""))
          
          # the sign is the sign of the effect of prodict of the stimulis and the inhibitor to the protein
          Sign<-CEN[i,2]
          
          Model <- AddLink(CueDown[ix_tmp1], SigUp[jx_tmp1], Model, Sign=Sign)
        }
      }
    }
  }
  Model$indexIntegr<-seq(from=indexIntegr+1, to=length(Model$reacID))
  return(Model)
}