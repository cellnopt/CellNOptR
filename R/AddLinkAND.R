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
AddLinkAND <-
function(Stimulus, Inhibitor, Signal, Model){
  
  # add a negative and gate to the model between the Sitimulus(-) and the Inhibitor(+)

  reaction<-paste("!",Stimulus,sep="")
  reaction<-paste(reaction,"+",Inhibitor,sep="")
  reaction<-paste(reaction,"=",Signal,sep="")
  if (!any(Model$reacID==reaction)){
    Model$reacID[length(Model$reacID)+1] <- reaction
    tmp<-matrix(data=0,nrow=dim(Model$interMat)[1], ncol=1)
    
	tmpNeg<-tmp
	souNeg<-match(Stimulus,Model$namesSpecies)
	tmpNeg[souNeg,1]<-1
	Model$notMat <- cbind(Model$notMat,tmpNeg)
	sou1<-match(Stimulus,Model$namesSpecies)
	sou2<-match(Inhibitor,Model$namesSpecies)
	tar<-match(Signal,Model$namesSpecies)
	tmp[sou1,1]<-(-1)
	tmp[sou2,1]<-(-1)
	tmp[tar,1]<-1
	Model$interMat <- cbind(Model$interMat,tmp)
	}
      
  colnames(Model$interMat) <- Model$reacID
  colnames(Model$notMat) <- Model$reacID
      
  return(Model)
  
}
