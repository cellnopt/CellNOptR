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

cutAndPlotResultsTimeScaleT2 <- function(
	model,
	bStringT1,
	bStringT2,
	simList,
	CNOlist,
	indexList,
	boolUpdates,
	divTime,
	lowerB=lowerB,
	upperB=upperB,
	show=TRUE,
	plotPDF=FALSE,
	tag=NULL) {
	
	library(abind)
	# simulate T1
	modelcut <- model
	modelcut$interMat <- modelcut$interMat[,as.logical(bStringT1)]
	modelcut$notMat <- modelcut$notMat[,as.logical(bStringT1)]
	modelcut$reacID <- modelcut$reacID[as.logical(bStringT1)]
	simListCut <- simList
	simListCut$finalCube <- simListCut$finalCube[as.logical(bStringT1),]
	simListCut$ixNeg <- simListCut$ixNeg[as.logical(bStringT1),]
	simListCut$ignoreCube <- simListCut$ignoreCube[as.logical(bStringT1),]
	simListCut$maxIx <- simListCut$maxIx[as.logical(bStringT1)]
	SimT1 <- simulatorTimeScale(CNOlist=CNOlist,model=modelcut,simList=simListCut,indexList=indexList, boolUpdates=boolUpdates[1])
	simResT1 <- SimT1[,indexList$signals,]
	getFitDataT1 <- getFitTimeScale(simList=simListCut, CNOlist=CNOlist, model=modelcut, indexList=indexList, boolUpdates=boolUpdates, divTime=divTime, lowerB=lowerB, upperB=upperB)
	xCoords1 <- getFitDataT1$xCoords

	# simulate T2
	bitString2 <- bStringT1
	bitString2[which(bStringT1 == 0)] <- bStringT2
	BStimes <- bStringT1
	BStimes[which(bStringT1 == 0)] <- bStringT2*2
	modelcut <- model
	modelcut$interMat <- modelcut$interMat[,as.logical(bitString2)]
	modelcut$notMat <- modelcut$notMat[,as.logical(bitString2)]
	modelcut$reacID <- modelcut$reacID[as.logical(bitString2)]
	modelcut$times <- BStimes[which(BStimes != 0)]
	simListCut <- simList
	simListCut$finalCube <- simListCut$finalCube[as.logical(bitString2),]
	simListCut$ixNeg <- simListCut$ixNeg[as.logical(bitString2),]
	simListCut$ignoreCube <- simListCut$ignoreCube[as.logical(bitString2),]
	simListCut$maxIx <- simListCut$maxIx[as.logical(bitString2)]
	SimT2 <- simulatorTimeScaleT2(simResultsT1=SimT1[,,dim(SimT1)[3]], CNOlist=CNOlist, model=modelcut, simList=simListCut, indexList=indexList, boolUpdates=boolUpdates[2])
	simResT2 <- SimT2[,indexList$signals,]
	getFitDataT2 <- getFitTimeScale(simList=simListCut, CNOlist=CNOlist, model=modelcut, indexList=indexList, boolUpdates=boolUpdates, divTime=divTime, timeSplit="late", simResultsT1=SimT1, lowerB=lowerB, upperB=upperB)
	xCoords2 <- CNOlist$timeSignals[which(CNOlist$timeSignals==divTime)+1] + getFitDataT2$xCoords
	
	# put it all together
	simResAll = abind(simResT1,simResT2,along=3)
	simResAll[,4,c(1:2)]=1
	simResAll[c(3,6,9),5,c(2:3)]=0
	yInterAll = abind(getFitDataT1$yInter, getFitDataT2$yInter, along=3)
	
	if(show==TRUE) {
		plotOptimResultsPan(
			simResults=simResAll,
			yInterpol=yInterAll,
			xCoords=c(xCoords1,xCoords2),
			CNOlist=CNOlist,
			formalism="dt"
		)	
	}
		
	if(plotPDF == TRUE) {
		if(is.null(tag)) {
			filename <- paste(deparse(substitute(model)),"simResultsT1.pdf",sep="")
        } else {
			filename <- paste(tag,"simResultsT1.pdf",sep="_")
        }
		plotOptimResultsPan(
			simResults=simResAll,
			yInterpol=yInterAll,
			xCoords=c(xCoords1,xCoords2),
			CNOlist=CNOlist,
			formalism="dt",
			pdfFileName=filename,
			pdf=TRUE)
	}
}

