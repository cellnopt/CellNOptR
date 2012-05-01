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

cutAndPlotResultsTimeScaleT2 <- function (Model, bStringT1, bStringT2, SimList, CNOlist, indexList, 
    boolUpdates, divTime, lowerB = lowerB, upperB = upperB) 
{
    library(abind)
    Modelcut <- Model
    Modelcut$interMat <- Modelcut$interMat[, as.logical(bStringT1)]
    Modelcut$notMat <- Modelcut$notMat[, as.logical(bStringT1)]
    Modelcut$reacID <- Modelcut$reacID[as.logical(bStringT1)]
    SimListcut <- SimList
    SimListcut$finalCube <- SimListcut$finalCube[as.logical(bStringT1), 
        ]
    SimListcut$ixNeg <- SimListcut$ixNeg[as.logical(bStringT1), 
        ]
    SimListcut$ignoreCube <- SimListcut$ignoreCube[as.logical(bStringT1), 
        ]
    SimListcut$maxIx <- SimListcut$maxIx[as.logical(bStringT1)]
    SimT1 <- simulatorTimeScale(CNOlist = CNOlist, Model = Modelcut, 
        SimList = SimListcut, indexList = indexList, boolUpdates = boolUpdates[1])
    SimResT1 <- SimT1[, indexList$signals, ]
    getFitDataT1 <- getFitTimeScale(SimList = SimListcut, CNOlist = CNOlist, 
        Model = Modelcut, indexList = indexList, boolUpdates = boolUpdates, 
        divTime = divTime, lowerB = lowerB, upperB = upperB)
    xCoords1 <- getFitDataT1$xCoords
    bitString2 <- bStringT1
    bitString2[which(bStringT1 == 0)] <- bStringT2
    BStimes <- bStringT1
    BStimes[which(bStringT1 == 0)] <- bStringT2 * 2
    Modelcut <- Model
    Modelcut$interMat <- Modelcut$interMat[, as.logical(bitString2)]
    Modelcut$notMat <- Modelcut$notMat[, as.logical(bitString2)]
    Modelcut$reacID <- Modelcut$reacID[as.logical(bitString2)]
    Modelcut$times <- BStimes[which(BStimes != 0)]
    SimListcut <- SimList
    SimListcut$finalCube <- SimListcut$finalCube[as.logical(bitString2), 
        ]
    SimListcut$ixNeg <- SimListcut$ixNeg[as.logical(bitString2), 
        ]
    SimListcut$ignoreCube <- SimListcut$ignoreCube[as.logical(bitString2), 
        ]
    SimListcut$maxIx <- SimListcut$maxIx[as.logical(bitString2)]
    SimT2 <- simulatorTimeScaleT2(SimResultsT1 = SimT1[, , dim(SimT1)[3]], 
        CNOlist = CNOlist, Model = Modelcut, SimList = SimListcut, 
        indexList = indexList, boolUpdates = boolUpdates[2])
    SimResT2 <- SimT2[, indexList$signals, ]
    getFitDataT2 <- getFitTimeScale(SimList = SimListcut, CNOlist = CNOlist, 
        Model = Modelcut, indexList = indexList, boolUpdates = boolUpdates, 
        divTime = divTime, timeSplit = "late", SimResultsT1 = SimT1, 
        lowerB = lowerB, upperB = upperB)
    xCoords2 <- CNOlist$timeSignals[which(CNOlist$timeSignals == 
        divTime) + 1] + getFitDataT2$xCoords
    SimResALL = abind(SimResT1, SimResT2, along = 3)
    SimResALL[, 4, c(1:2)] = 1
    SimResALL[c(3, 6, 9), 5, c(2:3)] = 0
    yInterALL = abind(getFitDataT1$yInter, getFitDataT2$yInter, 
        along = 3)
    expResults <- CNOlist$valueSignals
    plotOptimResultsTimeScale(SimResults = SimResALL, yInterpol = yInterALL, 
        xCoords = c(xCoords1, xCoords2), CNOlist = CNOlist)
}

