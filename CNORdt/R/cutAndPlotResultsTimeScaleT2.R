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
# $Id: $

cutAndPlotResultsTimeScaleT2 <- function (model, bStringT1, bStringT2, simList,
CNOlist, indexList, boolUpdates, divTime, lowerB=lowerB, upperB=upperB) {
    
    library(abind)
   	
   	modelCut <- model
    modelCut$interMat <- modelCut$interMat[, as.logical(bStringT1)]
    modelCut$notMat <- modelCut$notMat[, as.logical(bStringT1)]
    modelCut$reacID <- modelCut$reacID[as.logical(bStringT1)]
    
    simListCut <- simList
    simListCut$finalCube <- simListCut$finalCube[as.logical(bStringT1),]
    simListCut$ixNeg <- simListCut$ixNeg[as.logical(bStringT1),]
    simListCut$ignoreCube <- simListCut$ignoreCube[as.logical(bStringT1),]
    simListCut$maxIx <- simListCut$maxIx[as.logical(bStringT1)]
    
    simT1 <- simulatorTimeScale(CNOlist=CNOlist, model=modelCut, 
	simList=simListCut, indexList=indexList, boolUpdates=boolUpdates[1])
    simResT1 <- simT1[,indexList$signals,]
    
    getFitDataT1 <- getFitTimeScale(simList=simListCut, CNOlist=CNOlist, 
	model=modelCut, indexList=indexList, boolUpdates=boolUpdates, 
	divTime=divTime, lowerB=lowerB, upperB=upperB)
    
    xCoords1 <- getFitDataT1$xCoords
    bitString2 <- bStringT1
    bitString2[which(bStringT1 == 0)] <- bStringT2
    BStimes <- bStringT1
    BStimes[which(bStringT1 == 0)] <- bStringT2 * 2
    
    modelCut <- Model
    modelCut$interMat <- modelCut$interMat[, as.logical(bitString2)]
    modelCut$notMat <- modelCut$notMat[, as.logical(bitString2)]
    modelCut$reacID <- modelCut$reacID[as.logical(bitString2)]
    modelCut$times <- BStimes[which(BStimes != 0)]
    
    simListCut <- SimList
    simListCut$finalCube <- simListCut$finalCube[as.logical(bitString2),]
    simListCut$ixNeg <- simListCut$ixNeg[as.logical(bitString2),]
    simListCut$ignoreCube <- simListCut$ignoreCube[as.logical(bitString2),]
    simListCut$maxIx <- simListCut$maxIx[as.logical(bitString2)]
    
    simT2 <- simulatorTimeScaleT2(simResultsT1=SimT1[,,dim(simT1)[3]], 
	CNOlist=CNOlist, Model=modelCut, simList=simListCut, 
	indexList=indexList, boolUpdates=boolUpdates[2])
    simResT2 <- simT2[,indexList$signals,]
    
    getFitDataT2 <- getFitTimeScale(simList=simListCut, CNOlist=CNOlist, 
	model=modelCut, indexList=indexList, boolUpdates=boolUpdates, 
	divTime=divTime, timeSplit="late", simResultsT1=simT1, 
	lowerB=lowerB, upperB=upperB)
    
    xCoords2 <- CNOlist$timeSignals[which(CNOlist$timeSignals == 
	divTime) + 1] + getFitDataT2$xCoords
    simResAll = abind(simResT1, simResT2, along=3)

    yInterAll = abind(getFitDataT1$yInter, getFitDataT2$yInter, 
	along=3)
    
    plotOptimResultsTimeScale(simResults=simResAll, yInterpol=yInterAll, 
	xCoords=c(xCoords1, xCoords2), CNOlist=CNOlist)
}

