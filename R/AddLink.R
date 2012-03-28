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
# $Id$
AddLink <-
function(CueDown, SigUp, Model, Sign){
  
  # Update the vector of the name of reactions

  reaction<-paste(CueDown,"=",SigUp,sep="")
  if (Sign==(-1)){reaction<-paste("!",reaction,sep="")}
    
  if (!any(Model$reacID==reaction)){
    Model$reacID[length(Model$reacID)+1] <- reaction
    
	  tmp<-matrix(data=0,nrow=dim(Model$interMat)[1], ncol=1)
    sou<-match(CueDown,Model$namesSpecies)
	  if (Sign==(-1)){tmp[sou,1]<-1}
	  Model$notMat <- cbind(Model$notMat,tmp)
    
	  tmp<-matrix(data=0,nrow=dim(Model$interMat)[1], ncol=1)
	  sou<-match(CueDown,Model$namesSpecies)
    tar<-match(SigUp,Model$namesSpecies)
    tmp[sou,1]<-(-1)
	  tmp[tar,1]<-1
    Model$interMat <- cbind(Model$interMat,tmp)
  }
      
  colnames(Model$interMat) <- Model$reacID
  colnames(Model$notMat) <- Model$reacID
      
  return(Model)
  
}
