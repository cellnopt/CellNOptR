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

getFitDT <- function(simResults, CNOlist, model, indexList, sizeFac = 1e-04, NAFac = 1, nInTot, boolUpdates,  
    lowerB, upperB) {
    	    
    if ((class(CNOlist) == "CNOlist") == FALSE) {
        CNOlist = CellNOptR::CNOlist(CNOlist)
    }
      
    # dimensions, time points
    nTimes = length(CNOlist@timepoints)
    sigs = dim(CNOlist@signals[[1]])
    
    # cut simResults to view signals only
    simResults <- simResults[, indexList$signals, ]
    
    # interpolate experimental data so it can be compared to boolean simulation  
    splineStore = list()
    splineAdd = 1
    
    for (nExper in 1:dim(CNOlist@signals[[1]])[1]) {
        for (nSig in 1:dim(CNOlist@signals[[1]])[2]) {
            yTest = c()
            for (a in 1:nTimes) {
                yTest = c(yTest, CNOlist@signals[[a]][nExper, nSig])
            }
            
            if (!is.na(yTest[1])) {
                cS = splinefun(CNOlist@timepoints, yTest)
                splineStore[splineAdd] = list(cS)
            } else {
                cS = splinefun(CNOlist@timepoints, rep(0, nTimes))
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
    yInter = array(dim = dim(simResults))
    xCoords = seq(0, by = myEstimate$par, length.out = boolUpdates)
    
    count = 1
    for (nExper in 1:dim(simResults)[1]) {
        for (nSig in 1:dim(simResults)[2]) {
            
            yOut = splineStore[[count]](xCoords)
            yInter[nExper, nSig, ] = yOut
            count = count + 1
        }
    }
    
    diff <- (simResults - yInter)
    r <- diff^2
    deviationPen <- sum(r[!is.na(r)])/nTimes
    NAPen <- NAFac * length(which(is.na(simResults)))
    nDataPts <- dim(CNOlist@signals[[1]])[1] * dim(CNOlist@signals[[1]])[2] * nTimes
    nInputs <- length(which(model$interMat == -1))
    
    # nInTot: number of inputs of expanded model nInputs: number of inputs of cut model
    sizePen <- (nDataPts * sizeFac * nInputs)/nInTot
    
    score <- deviationPen + NAPen + sizePen
    return(list(score = score, estimate = myEstimate$par, xCoords = xCoords, yInter = yInter, simResults = simResults))
    
}

