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
addLinkAND <-
function(Stimulus, Inhibitor, Signal, model){
  
  # add a negative and gate to the model between the Sitimulus(-) and the Inhibitor(+)

  reaction<-paste("!",Stimulus,sep="")
  reaction<-paste(reaction,"+",Inhibitor,sep="")
  reaction<-paste(reaction,"=",Signal,sep="")
  if (!any(model$reacID==reaction)){
    model$reacID[length(model$reacID)+1] <- reaction
    tmp<-matrix(data=0,nrow=dim(model$interMat)[1], ncol=1)
    
	tmpNeg<-tmp
	souNeg<-match(Stimulus,model$namesSpecies)
	tmpNeg[souNeg,1]<-1
	model$notMat <- cbind(model$notMat,tmpNeg)
	sou1<-match(Stimulus,model$namesSpecies)
	sou2<-match(Inhibitor,model$namesSpecies)
	tar<-match(Signal,model$namesSpecies)
	tmp[sou1,1]<-(-1)
	tmp[sou2,1]<-(-1)
	tmp[tar,1]<-1
	model$interMat <- cbind(model$interMat,tmp)
	}
      
  colnames(model$interMat) <- model$reacID
  colnames(model$notMat) <- model$reacID
      
  return(model)
  
}
