#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id$
mapBTables2model <-
function(BTable,model=NULL,optimRes=NA,allInter=TRUE, compressed=TRUE){
  # BTable are the Bolean Tables inferred from data using the function makeBTables
  # Model is the model optimised using CNO and than compressed
  
  if (!is.null(model)){
    g<-sif2graph(model2sif(model=model, optimRes=optimRes))
    #the vectr indexIntegr will contain the indexes of the added links
    indexIntegr<-length(model$reacID)
  }else{
    g<-new("graphNEL", edgemode="directed")
    indexIntegr<-0
    model<-list()
    model$reacID<-vector()
    model$namesSpecies<-unique(c(colnames(BTable$tables[[1]]), rownames(BTable$tables[[1]]), BTable$namesSignals))
    model$interMat<-matrix(data=0,nrow=length(model$namesSpecies), ncol=0)
    rownames(model$interMat)<-model$namesSpecies
	  model$notMat<-matrix(data=0,nrow=length(model$namesSpecies), ncol=0)
	  rownames(model$notMat)<-model$namesSpecies
  }
  
  namesStimuli<-colnames(BTable$tables[[1]])
  namesInhibitors<-rownames(BTable$tables[[1]])
  namesSignals<-BTable$namesSignals
  #Others are the white nodes in the network, those that are not stimulated, nor inhibited, nor measures
  namesOthers<-setdiff(setdiff(setdiff(model$namesSpecies, namesSignals), namesInhibitors), namesStimuli)
  
    
  ck=0
  # i signals
  # j stimuli
  # k inhibitors
  for (i in 1:length(namesSignals)){   
    for (j in 1:length(namesStimuli)){
      
      print(paste("---------------------------------------------"))
      print(paste('stimulus: ', namesStimuli[j], '     signal: ', namesSignals[i]))
      print(paste("Added links:"))
      
            
      # noNodes are the inhibitors that in the Btable of signal i, have a 1 in the column
      # corresponding to simulus j, meaning that the inhibitor has no effect
      noNodes <- rownames(as.matrix(which(BTable$tables[[i]][,j]==1)))
      if (compressed==TRUE){
        noNodes <- union(noNodes, namesOthers)
      }
      noNodes <- setdiff(noNodes,BTable$namesSignals[i])
      
      ix <- as.vector(which(BTable$tables[[i]][,j]==2))
      if (length(ix) > 0){
        for (k in 1:length(ix)){
          # SignStim is the sign of the effect of the Stimulus on the Signal
          SignStim<-1
          if (BTable$NotMatStim[[i]][ix[k],j]==1) {SignStim<-(-1)}
                
          # SignInhib is the sign of the effect of the Inhibitor on the Signal
          SignInhib<-1
          if (BTable$NotMatInhib[[i]][ix[k],j]==1) {SignInhib<-(-1)}
          
          # connection between the stimulus and the inhibitor
          cue <- colnames(BTable$tables[[i]])[j]
          signal <- rownames(BTable$tables[[i]])[ix[k]]
          
          ck <- searchLinkGraph(node1 = cue,node2 = signal, graph=g, noNodes=noNodes)
          #ck <- SearchLink(node1 = cue,node2 = signal, noNodes=noNodes, model = model)
          if (ck == 0){
            # if there is no connection between the cue and the signal in the model
            # all possible connections are added between nodes downstream the cue
            # (but upstream the following cue) and those upstream the signal
            # (but downstream the previous signal)
            if (allInter == TRUE && length(graph::nodes(g))>0){
              # 1. create a vector containing all nodes downstream the cue
              # (until reaching a cue or a noNode)
              AllCues <- union(namesStimuli,namesInhibitors)
              CueDown <- downCueGraph(cue=cue, graph=g, stopNodes=union(AllCues,noNodes))
            
              # 2. create a vector containing all nodes upstream the signal
              # (until reaching a signal, a Stimulus, a noNode)
              # AllSignals <- BTable$namesSignals
              stopNodes <-unique(union(namesSignals, noNodes))
              stopNodes<-unique(union(stopNodes,namesStimuli))
              SigUp <- upSignalGraph(signal=signal, graph=g, stopNodes=stopNodes)

            }else{
              CueDown <- as.vector(cue)
              SigUp <- as.vector(signal)    
            }
            for (ix_tmp1 in 1:length(CueDown)){
              for (jx_tmp1 in 1:length(SigUp)){
                print(paste(CueDown[ix_tmp1], "->", SigUp[jx_tmp1],sep=""))

                # the sign is the sign of the effect of prodict of the stimulis and the inhibitor to the protein
                Sign<-SignStim*SignInhib
                
                model <- addLink(CueDown[ix_tmp1], SigUp[jx_tmp1], model, Sign=Sign)
                #g<-sif2graph(model2sif(model=model))
              }
            }
          }
    
      
          # connection between the inhibitor and the signal
          cue <- rownames(BTable$tables[[i]])[ix[k]]
          signal <- BTable$namesSignals[i]
          #ck <- searchLink(node1 = cue,node2 = signal, noNodes=noNodes, model = model)
          ck <- searchLinkGraph(node1 = cue,node2 = signal, graph=g, noNodes=noNodes)
          if (ck == 0){
            # if there is no connection between the cue and the signal in the model
            # all possible connections are added between nodes downstream the cue
            # (but upstream the following cue) and those upstream the signal
            # (but downstream the previous signal)
            if (allInter == TRUE && length(graph::nodes(g))>0){
              # 1. create a vector containing all nodes downstream the cue
              # (until reaching a cue or a noNode)
              AllCues <- union(namesStimuli,namesInhibitors)
              CueDown <- downCueGraph(cue=cue, graph=g, stopNodes=union(AllCues,noNodes))
            
              # 2. create a vector containing all nodes upstream the signal
              # (until reaching a signal, a Stimulus, a noNode)
              # AllSignals <- BTable$namesSignals
              stopNodes <-unique(union(namesSignals, noNodes))
              stopNodes<-unique(union(stopNodes,namesStimuli))
              SigUp <- upSignalGraph(signal=signal, graph=g, stopNodes=stopNodes)
            }else{
              CueDown <- as.vector(cue)
              SigUp <- as.vector(signal)    
            }
            for (ix_tmp1 in 1:length(CueDown)){
              for (jx_tmp1 in 1:length(SigUp)){
                print(paste(CueDown[ix_tmp1], "->", SigUp[jx_tmp1],sep=""))
                
                # the sign is the sign of the effect of the inhibitor to the protein
                Sign<-SignInhib
                
                model <- addLink(CueDown[ix_tmp1], SigUp[jx_tmp1], model, Sign=Sign)
                #g<-sif2graph(model2sif(model=model))
              }
            }      
          }
          
          ###!!!!!!!!!!!!
          ## ADD HERE GENERALIZED AND GATES
          # possible negative and gates between the inhibited prot and one signal
          # check if the inhibited prot is also measured
          ixInh <- which(namesSignals==rownames(BTable$tables[[i]])[ix[k]])
          if (length(ixInh)>0){
            # if yes, look at which stimuli are affecting that signal (though that inhibitor or not
		    #ixSti1 <- which(BTable$tables[[i]][ix[k],]==2)
            ixSti1 <- which(BTable$tables[[i]][ix[k],]>0)
			# and which stimuli are affecting the inhibited (and measured) protein
            ixSti2 <- which(BTable$tables[[ixInh]][ix[k],]>0)
            ixStiDiff <- setdiff(ixSti2, ixSti1)
            # if they do not coincide, add an AND gate between the stimulus that is not affecting the protien
            # (but only the inhibitor) and the inhibitor - with negative sign for the stimulus
            if (length(ixStiDiff)>0){
              for (k2 in 1:length(ixStiDiff)){
                ck <- searchLinkGraph(node1 = cue,node2 = signal, graph=g, noNodes=noNodes)
				print(paste("!",namesStimuli[ixStiDiff[k2]], "+", namesInhibitors[ix[k]] ,"->", namesSignals[i],sep=""))
                model<-addLinkAND(namesStimuli[ixStiDiff[k2]], namesInhibitors[ix[k]], namesSignals[i], model)
              }
            }
          }
          ###!!!!!!!!!!!!
        }
      }
              
      # additional connection between the stimulus and the signal
	    if (sum(BTable$tables[[i]][,j]) == dim(BTable$tables[[i]])[1]){
	      cue <- colnames(BTable$tables[[i]])[j]
        signal <- BTable$namesSignals[i]
        
        SignStim<-1
        if (BTable$NotMatStim[[i]][1,j]==1) {SignStim<-(-1)}
        
        #ck <- SearchLink(node1 = cue,node2 = signal, noNodes=noNodes, model = model)
        ck <- searchLinkGraph(node1 = cue,node2 = signal, graph=g, noNodes=noNodes)
        if (ck == 0){
          # if there is no connection between the cue and the signal in the model
          # all possible connections are added between nodes downstream the cue
          # (but upstream the following cue) and those upstream the signal
          # (but downstream the previous signal)
          if (allInter == TRUE && length(graph::nodes(g))>0){
            # 1. create a vector containing all nodes downstream the cue
            # (until reaching a cue or a noNode)
            AllCues <- union(namesStimuli,namesInhibitors)
            CueDown <- downCueGraph(cue=cue, graph=g, stopNodes=union(AllCues,noNodes))
            
            # 2. create a vector containing all nodes upstream the signal
            # (until reaching a signal, a Stimulus, a noNode)
            # AllSignals <- BTable$namesSignals
            stopNodes <-unique(union(namesSignals, noNodes))
            stopNodes<-unique(union(stopNodes,namesStimuli))
            SigUp <- upSignalGraph(signal=signal, graph=g, stopNodes=stopNodes)            
          }else{
            CueDown <- as.vector(cue)
            SigUp <- as.vector(signal)    
          }
	        for (ix_tmp1 in 1:length(CueDown)){
            for (jx_tmp1 in 1:length(SigUp)){
              print(paste(CueDown[ix_tmp1], "->", SigUp[jx_tmp1],sep=""))
              
              model <- addLink(CueDown[ix_tmp1], SigUp[jx_tmp1], model, Sign=SignStim)
              #g<-sif2graph(model2sif(model=model))
            }
          }
        }
	    }
    print(paste("---------------------------------------------"))
    print(paste(""))
    }
  }
  model$indexIntegr<-seq(from=indexIntegr+1, to=length(model$reacID))
  return(model)
}
