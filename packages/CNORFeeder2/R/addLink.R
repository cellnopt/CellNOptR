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
addLink <-
function(CueDown, SigUp, model, Sign){
  
  # Update the vector of the name of reactions

  reaction<-paste(CueDown,"=",SigUp,sep="")
  if (Sign==(-1)){reaction<-paste("!",reaction,sep="")}
    
  if (!any(model$reacID==reaction)){
    model$reacID[length(model$reacID)+1] <- reaction
    
	  tmp<-matrix(data=0,nrow=dim(model$interMat)[1], ncol=1)
    sou<-match(CueDown,model$namesSpecies)
	  if (Sign==(-1)){tmp[sou,1]<-1}
	  model$notMat <- cbind(model$notMat,tmp)
    
	  tmp<-matrix(data=0,nrow=dim(model$interMat)[1], ncol=1)
	  sou<-match(CueDown,model$namesSpecies)
    tar<-match(SigUp,model$namesSpecies)
    tmp[sou,1]<-(-1)
	  tmp[tar,1]<-1
    model$interMat <- cbind(model$interMat,tmp)
  }
      
  colnames(model$interMat) <- model$reacID
  colnames(model$notMat) <- model$reacID
      
  return(model)
  
}
