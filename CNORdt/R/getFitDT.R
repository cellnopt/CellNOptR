# This file is part of the CNO software
# 
# Copyright (c) 2011-2012 - EBI
# 
# File author(s): CNO developers (cno-dev@ebi.ac.uk)
# 
# Distributed under the GPLv2 License.  See accompanying file LICENSE.txt or copy at
# http://www.gnu.org/licenses/gpl-2.0.html
# 
# CNO website: http://www.ebi.ac.uk/saezrodriguez/software.html
# 
# $Id: $

getFitDT <- function(simResults, CNOlist, model, indexList, sizeFac = 1e-04, NAFac = 1, nInTot, boolUpdates, divTime = NULL, 
    lowerB, upperB) {
    
    # ignore for the moment and return to for TN if(is.null(divTime)) {divTime =
    # CNOlist$timeSignals[length(CNOlist$timeSignals)]} if(timeSplit=='early') { timeExper =
    # CNOlist$timeSignals[CNOlist$timeSignals <= divTime] timeIndex = which(CNOlist$timeSignals <= divTime) boolUpdates =
    # boolUpdates[1] } else if(timeSplit=='late') { timeExper = CNOlist$timeSignals[CNOlist$timeSignals > divTime] timeIndex
    # = which(CNOlist$timeSignals > divTime) boolUpdates = boolUpdates[2] } if ((class(CNOlist)=='CNOlist')==FALSE){ CNOlist
    # = CellNOptR::CNOlist(CNOlist) } make sure timeExper[1] = 0 for fitting
    timeExper = CNOlist$timeSignals - CNOlist$timeSignals[1]
    
    # dimensions, time points
    times = length(timeExper)
    sigs = dim(CNOlist$valueSignals[[1]])
    
    # simulator will need to add simulatorDT-TN here, ignore for now if(timeSplit=='early') { yBool =
    # simulatorTimeScale(CNOlist, model, simList, indexList, boolUpdates) } else if(timeSplit=='late') { yBool =
    # simulatorTimeScaleT2(simResultsT1[,,dim(simResultsT1)[3]], CNOlist, model, simList, indexList, boolUpdates) }
    
    simResults <- simResults[, indexList$signals, ]
    
    
    splineStore = list()
    splineAdd = 1
    
    for (nExper in 1:dim(CNOlist$valueSignals[[1]])[1]) {
        for (nSig in 1:dim(CNOlist$valueSignals[[1]])[2]) {
            yTest = c()
            for (a in 1:times) {
                yTest = c(yTest, CNOlist$valueSignals[[a]][nExper, nSig])
            }
            
            if (!is.na(yTest[1])) {
                cS = splinefun(timeExper, yTest)
                splineStore[splineAdd] = list(cS)
            } else {
                cS = splinefun(timeExper, rep(0, times))
                splineStore[splineAdd] = list(cS)
            }
            splineAdd = splineAdd + 1
        }
    }
    
    # what to optimize
    taufinder <- function(whatScale) {
        
        ySilico = array(dim = dim(simResults))
        numberPoints = dim(simResults)[3]
        xCoords = seq(0, by = whatScale, length.out = numberPoints)
        count = 1
        
        for (nExper in 1:dim(simResults)[1]) {
            for (nSig in 1:dim(simResults)[2]) {
                
                yOut = splineStore[[count]](xCoords)
                ySilico[nExper, nSig, ] = yOut
                count = count + 1
            }
        }
        
        ErrorVector = as.vector(ySilico) - as.vector(simResults)
        sse = sum(ErrorVector^2)
        return(sse)
    }
    
    seed = 0.99
    myEstimate = optim(seed, taufinder, method = "L-BFGS-B", lower = lowerB, upper = upperB)
    yFinal = array(dim = dim(simResults))
    xCoords = seq(0, by = myEstimate$par, length.out = boolUpdates)
    
    count = 1
    for (nExper in 1:dim(simResults)[1]) {
        for (nSig in 1:dim(simResults)[2]) {
            
            yOut = splineStore[[count]](xCoords)
            yFinal[nExper, nSig, ] = yOut
            count = count + 1
        }
    }
    
    diff <- (simResults - yFinal)  # * my.weights
    r <- diff^2
    deviationPen <- sum(r[!is.na(r)])/times
    NAPen <- NAFac * length(which(is.na(simResults)))
    
    # CHANGE HERE
    nDataPts <- dim(CNOlist$valueSignals[[1]])[1] * dim(CNOlist$valueSignals[[1]])[2]
    
    nInputs <- length(which(model$interMat == -1))
    
    # nInTot: number of inputs of expanded model nInputs: number of inputs of cut model
    sizePen <- (nDataPts * sizeFac * nInputs)/nInTot
    
    score <- deviationPen + NAPen + sizePen
    return(list(score = score, estimate = myEstimate$par, xCoords = xCoords, yInter = yFinal, simResults = simResults))
    
}

convert2array <- function(x, nRow, nCol, nBool) {
    v1 = c(x)
    count = 1
    out1 = array(NA, dim = c(nRow, nCol, nBool))
    for (d in 1:nBool) {
        for (a in 1:nRow) {
            for (b in 1:nCol) {
                out1[a, b, d] = v1[count]
                count = count + 1
            }
        }
    }
    
    return(out1)
}