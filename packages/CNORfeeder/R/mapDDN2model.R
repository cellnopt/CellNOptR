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
# $Id: mapDDN2model.R 853 2012-03-28 15:09:06Z eduati $
mapDDN2model <-
function(DDN,model,CNOlist,allInter=TRUE){
	
	if ((class(CNOlist)=="CNOlist")==FALSE){
		CNOlist = CellNOptR::CNOlist(CNOlist)
	}
	
	g<-sif2graph(model2sif(model=model))
	indexIntegr<-length(model$reacID)
	
	namesStimuli<-colnames(CNOlist@stimuli)
	namesInhibitors<-colnames(CNOlist@inhibitors)
	namesSignals<-colnames(CNOlist@signals[[2]])
#Others are the white nodes in the network, those that are not stimulated, nor inhibited, nor measures
	namesOthers<-setdiff(setdiff(setdiff(model$namesSpecies, namesSignals), namesInhibitors), namesStimuli)
	
	for (i in 1:dim(DDN)[1]){
		noNodes <- namesOthers
		node1<-DDN[i,1]
		node2<-DDN[i,3]
		
		print(paste("---------------------------------------------"))
		print(paste('stimulus: ', node1, '     signal: ', node2))
		print(paste("Added links:"))
		
		
		ck <- searchLinkGraph(node1 = node1,node2 = node2, graph=g, noNodes=noNodes)
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
					Sign<-DDN[i,2]
					
					model <- addLink(CueDown[ix_tmp1], SigUp[jx_tmp1], model, Sign=Sign)
				}
			}
		}
	}
	model$indexIntegr<-seq(from=indexIntegr+1, to=length(model$reacID))
	return(model)
}
