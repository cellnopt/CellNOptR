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

cutAndPlotResultsTimeScale <- function(
	Model,
	bString,
	SimList,
	CNOlist,
	indexList,
	boolUpdates=boolUpdates,
	divTime=NULL,
	lowerB=lowerB,
	upperB=upperB, 
	show=TRUE,
	plotPDF=FALSE,
	tag=NULL) {
	
	Modelcut <- Model
	Modelcut$interMat <- Modelcut$interMat[,as.logical(bString)]
	Modelcut$notMat <- Modelcut$notMat[,as.logical(bString)]
	Modelcut$reacID <- Modelcut$reacID[as.logical(bString)]
	SimListcut <- SimList
	SimListcut$finalCube <- SimListcut$finalCube[as.logical(bString),]
	SimListcut$ixNeg <- SimListcut$ixNeg[as.logical(bString),]
	SimListcut$ignoreCube <- SimListcut$ignoreCube[as.logical(bString),]
	SimListcut$maxIx <- SimListcut$maxIx[as.logical(bString)]
	
	boolUpdates = boolUpdates[1]
	SimRes <- simulatorTimeScale(CNOlist=CNOlist, Model=Modelcut, SimList=SimListcut, indexList=indexList, boolUpdates=boolUpdates)
	SimRes = SimRes[,indexList$signals,]
	getFitData <- getFitTimeScale(SimList=SimListcut, CNOlist=CNOlist, Model=Modelcut, indexList=indexList, boolUpdates=boolUpdates, divTime=divTime, lowerB=lowerB, upperB=upperB)
	
	if(show==TRUE) {
		plotOptimResultsPan(SimResults=SimRes,
		yInterpol=getFitData$yInter,
		xCoords=getFitData$xCoords,
		CNOlist=CNOlist,
		formalism="dt"
		)
	}
	if(plotPDF == TRUE) {
		if(is.null(tag)) {
			filename <- paste(deparse(substitute(Model)),"SimResultsT1.pdf",sep="")
        } else {
			filename <- paste(tag,"SimResultsT1.pdf",sep="_")
        }
		plotOptimResultsPan(SimResults=SimRes,
			yInterpol=getFitData$yInter,
			xCoords=getFitData$xCoords,
			CNOlist=CNOlist,
			formalism="dt",
			pdfFileName=filename,
			pdf=TRUE
		)
	}
}

