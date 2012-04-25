#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv2 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-2.0.html
#
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/software.html
#
##############################################################################
# $Id: PPIweight.R 853 2012-03-28 15:09:06Z cokelaer $
PPIweight <-
  function(modelIntegr,PKNmodel,CNOlist,UniprotID,PPINigraph){
    
	count<-0
	saveShortestPath<-list()
	  
    #vector with the scores (score =1 means that links will have the same weight as the other links, final score will be 1+scorePPI)
    scoresVec<-rep(1,length(modelIntegr$indexIntegr))
    
    # transform the PKN CNOmodel in a graph
    PKNgraph<-sif2graph(model2sif(Model=PKNmodel, optimRes=NA))
    # I want to add the weight only to the integrated links (can be easily changed to all links)
    links2weight<-modelIntegr$reacID[modelIntegr$indexIntegr]
    
    namesCues<-CNOlist$namesCues
    namesStimuli<-CNOlist$namesStimuli
    namesSignals<-CNOlist$namesSignals
  
    # for each link
	ListPaths<-list()
    for (i in 1:length(links2weight)){
	  ListPaths[[i]]<-list()
	  saveShortestPath[[i]]<-vector()
      x<-links2weight[i]
      x<-gsub("!", "", x)
      x<-unlist(strsplit(x, "=", fixed=TRUE))
      signal<-x[2]
      # can contain more than 1 element in case of AND nodes
      cues<-unlist(strsplit(x[1], "+", fixed=TRUE))
      CueDown<-list()
      for (j in 1:length(cues)){
        # 1. create a vector for each cue containing all nodes downstream the cue
        # (until reaching a cue or a noNode)
        cue<-cues[j]
		CueDown[[j]] <- downCueGraph(cue=cue, graph=PKNgraph, stopNodes=namesCues)
      }
      # 2. create a vector containing all nodes upstream the signal
      # (until reaching a signal, a Stimulus)
      stopNodes <-unique(union(namesSignals,namesStimuli))
      SigUp <- upSignalGraph(signal=signal, graph=PKNgraph, stopNodes=stopNodes)
      node1<-as.list(SigUp)
      
      # ATTENTION: valid only for AND nodes with 2 input nodes (but for now I cannot infer AND with more than 2 inputs)
      # if there are AND nodes I have to create a vector of node2 with all possible couples of nodes
      if (length(CueDown)>1){
        node2<-list()
        count<-1
        for (j1 in 1:length(CueDown[[1]])){
          for (j2 in 1:length(CueDown[[2]])){
            node2[[count]]<-c(CueDown[[1]][j1],CueDown[[2]][j2])
            count<-count+1
          }
        }
      }else{
		#node2<-as.list(CueDown)    NOTE: THIS ONE FOR OLD WEIGHTS
        node2<-as.list(CueDown[[1]])
      }
      
      score<-0
	    cSP<-1
      for (j1 in 1:length(node1)){
        for (j2 in 1:length(node2)){
          res<-searchPPIigraph(node1=node1[[j1]], node2=node2[[j2]], UniprotID=UniprotID, igraph=PPINigraph, noPKNnodes=TRUE)
		      PPINigraph<-res$igraph

          #(res$ckmin+1) is the number of edges
		      scoreAdd<-res$ckmin+1
			    if (scoreAdd==0){scoreAdd<-1}
			    score<-score+(1/(scoreAdd))

			    ListPaths[[i]][[cSP]]<-res$path
			    saveShortestPath[[i]][cSP]<-scoreAdd
			    cSP<-cSP+1
        }
      }
      
      scoresVec[i]<-1+1/score
    }
    
    modelIntegr$IntegrPPIscores<-scoresVec
    return(list(modelIntegr=modelIntegr,PPINigraph=PPINigraph,saveShortestPath=saveShortestPath, ListPaths=ListPaths))
  }