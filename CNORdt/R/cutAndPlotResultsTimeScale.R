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

cutAndPlotResultsTimeScale <- function (model, bString, simList, CNOlist, indexList, boolUpdates=boolUpdates, 
divTime=NULL, lowerB=lowerB, upperB=upperB) {

	modelCut <- model
	modelCut$interMat <- modelCut$interMat[, as.logical(bString)]
	modelCut$notMat <- modelCut$notMat[, as.logical(bString)]
	modelCut$reacID <- modelCut$reacID[as.logical(bString)]
	
	simListCut <- simList
	simListCut$finalCube <- simListCut$finalCube[as.logical(bString),]
	simListCut$ixNeg <- simListCut$ixNeg[as.logical(bString),]
	simListCut$ignoreCube <- simListCut$ignoreCube[as.logical(bString),]

	simListCut$maxIx <- simListCut$maxIx[as.logical(bString)]
	boolUpdates = boolUpdates[1]
	simRes <- simulatorTimeScale(CNOlist = CNOlist, model = modelCut, 
	simList = simListCut, indexList = indexList, boolUpdates = boolUpdates)
	simRes = simRes[, indexList$signals, ]
	getFitData <- getFitTimeScale(simList = simListCut, CNOlist = CNOlist, 
	model = modelCut, indexList = indexList, boolUpdates = boolUpdates, 
	divTime = divTime, lowerB = lowerB, upperB = upperB)
    
    plotOptimResultsTimeScale(simResults = simRes, yInterpol = getFitData$yInter, 
	xCoords = getFitData$xCoords, CNOlist = CNOlist)
}

