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
simulateTN<-function(CNOlist,model,simResT1,bStringT2,bStringTimes,simList,indexList,timeIndex){

#Made to look more like new SimulateT1 (6/26/12 edits) but old version is left in comments below
  
  modelCut<-cutModel(model, bStringT2)
  
  
  simListCut<-simList
  simListCut$finalCube<-simListCut$finalCube[as.logical(bStringT2),]
  simListCut$ixNeg<-simListCut$ixNeg[as.logical(bStringT2),]
  simListCut$ignoreCube<-simListCut$ignoreCube[as.logical(bStringT2),]
  simListCut$maxIx<-simListCut$maxIx[as.logical(bStringT2)]
  
  if(is.null(dim(simListCut$finalCube))){
    simListCut$finalCube<-matrix(simListcut$finalCube,ncol=1)
    simListCut$ixNeg<-matrix(simListcut$ixNeg,ncol=1)
    simListCut$ignoreCube<-matrix(simListcut$ignoreCube,ncol=1)
  }
  
  
#   bitString2<-bStringT2
#   modelCut <- model
#   modelCut$interMat <- modelCut$interMat[,as.logical(bitString2)]
#   modelCut$notMat <- modelCut$notMat[,as.logical(bitString2)]
#   modelCut$reacID <- modelCut$reacID[as.logical(bitString2)]
  modelCut$times <- bStringTimes[which(bStringTimes != 0)]
  
#   simListCut <- cutSimList(simList, bitString2)
  
  # simulate
  SimT2 <- simulatorTN(
    simResultsT1=simResT1,
    CNOlist=CNOlist,
    model=modelCut,
    simList=simListCut,
    indexList=indexList,
    timeIndex=timeIndex)
  
  return(SimT2)
}
