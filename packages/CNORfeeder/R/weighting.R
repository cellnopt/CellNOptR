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
# $Id: weighting.R 853 2012-03-28 15:09:06Z eduati $
weighting <-
  function(modelIntegr,PKNmodel,CNOlist,integrFac=10,UniprotID=NULL,PPI=NULL){
    
    
	  
	if ((class(CNOlist)=="CNOlist")==FALSE){
		CNOlist = CellNOptR::CNOlist(CNOlist)
	}
	  
	# create a field with liks weight/reliability and allows different (higher) weigth for integrated links according to integrFac
	modelIntegr$linksWeights <- rep(1,length(modelIntegr$reacID))
	modelIntegr$linksWeights[modelIntegr$indexIntegr]<-integrFac * modelIntegr$linksWeights[modelIntegr$indexIntegr]
	
	# if PPI is set to FALSE we only add the field linksWeights (with additional penalty for integrated links) and stop here
	# if PPI is set to TRUE we do the weighting based on information from protein-protein interaction network
	if (!is.null(PPI)){
	  
	  requireNamespace("igraph")
		
		PPINigraph <- PPI
		if(igraph::is.igraph(PPI)==FALSE){
			stop("The provider PPI is not an igraph")
        }
		PPINigraph <- PPI
				
		count<-0
		#vector with the scores only for integrated links (score =1 means that links will have the same weight as the other links, final score will be 1+scorePPI)
		scoresVec<-rep(1,length(modelIntegr$indexIntegr))
		
		# transform the PKN CNOmodel in a graph
		PKNgraph<-sif2graph(model2sif(model=PKNmodel, optimRes=NA))
		# I want to add the weight only to the integrated links (can be easily changed to all links)
		links2weight<-modelIntegr$reacID[modelIntegr$indexIntegr]
		
		namesCues<-colnames(CNOlist@cues)
		namesStimuli<-colnames(CNOlist@stimuli)
		namesSignals<-colnames(CNOlist@signals[[1]])
		
		# for each link
		for (i in 1:length(links2weight)){
			
			x<-links2weight[i]
			print(x)
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
				whiteNodes <- setdiff(modelIntegr$namesSpecies, union(namesSignals,namesCues))
				stopNodes<-unique(union(namesCues, whiteNodes))
				CueDown[[j]] <- downCueGraph(cue=cue, graph=PKNgraph, stopNodes=stopNodes)
				#CueDown[[j]] <- downCueGraph(cue=cue, graph=PKNgraph, stopNodes=namesCues)  #before
			}
			# 2. create a vector containing all nodes upstream the signal
			# (until reaching a signal, a Stimulus) 
			#stopNodes <-unique(union(namesSignals,namesStimuli))  #before
			stopNodes <-unique(union(namesSignals,namesCues))
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
			
			for (j1 in 1:length(node1)){
				for (j2 in 1:length(node2)){
				
					
					ckmin<-searchPPIigraph(node1=node1[[j1]], node2=node2[[j2]], UniprotID=UniprotID, PPINigraph=PPINigraph, noPKNnodes=TRUE)

					#(ckmin+1) is the number of edges
					scoreAdd<-ckmin+1
					if (scoreAdd==0){scoreAdd<-1}
					score<-score+(1/(scoreAdd))
				}
			}
			
			scoresVec[i]<-1+1/score
		}
		
		modelIntegr$linksWeights[modelIntegr$indexIntegr]<-modelIntegr$linksWeights[modelIntegr$indexIntegr] * scoresVec
		
	}
	
    return(modelIntegr)
    
  }
