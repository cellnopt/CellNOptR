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
# $Id$

getFitDelay <- function(CNOlist, model, simList, indexList, sizeFac = 1e-04, NAFac = 1, nInTot, boolUpdates) {
    	    
    if ((class(CNOlist) == "CNOlist") == FALSE) {
        CNOlist = CellNOptR::CNOlist(CNOlist)
    }
      
    # dimensions, time points
    nTimes = length(CNOlist@timepoints)
    sigs = dim(CNOlist@signals[[1]])
    
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
	
   	ySilico = array(dim=c(dim(CNOlist@signals[[1]])[1], length(indexList$signals), boolUpdates))
   	xCoords = seq(0, by=boolUpdates/CNOlist@timepoints[length(CNOlist@timepoints)], length.out=boolUpdates)

   	count = 1
   	for(nExper in 1:dim(ySilico)[1]) {
    	for(nSig in 1:dim(ySilico)[2]) {			             
            yOut = splineStore[[count]](xCoords)
            ySilico[nExper,nSig,] = yOut;
            count = count+1
   		}
  	}
	
	################################################################################
	
	# TODO: neg data should be index list from model$reacID
	# of what strongWeak edges there are
	strongWeak = rep(0,length(model$reacID))
	negEdges = feedbackWrapper(model) # the index of edges that can be changed

	# optimization
	# what to optimize
	delayFinder <- function(optimStringF) {
		
	#	whatScale = optimStringF[1] # continuous with upper/lower bounds
		delayThresh = optimStringF[1:length(model$reacID)] # integer section of optimString
	#	strongWeak = c(0,0,0,0,1,0,0,0,0,0,0)
		strongWeak[negEdges] = optimStringF[(length(model$reacID)+1):length(optimStringF)] # binary section of optimString
		
		simResults = simulatorDelay(CNOlist, model, simList, indexList,
		boolUpdates=boolUpdates, delayThresh=delayThresh, strongWeak=strongWeak)
		simResults = convert2array(simResults, dim(CNOlist@signals[[1]])[1], length(model$namesSpecies), boolUpdates)
		simResults = simResults[,indexList$signals,]

  		errorVector = as.vector(ySilico) - as.vector(simResults)
    	sse = sum(errorVector^2)
    	return(sse)			
	}

	optimString = c(rep(1,length(model$reacID)), rep(0,length(negEdges)))
	problem = list(
		f=delayFinder,
		x_O = optimString,
		x_L=c(rep(0,length(model$reacID)), rep(0,length(negEdges))),
		x_U=c(rep(10,length(model$reacID)), rep(1,length(negEdges))),
		int_var=length(model$reacID) + length(negEdges)
	#	bin_var=length(model$reacID)
	)

	opts = list(maxtime=15, local_solver=0, dim_refset=6, ndiverse=50)
	est = essR(problem, opts)	

	################################################################################

	bestDelay = est$xbest[1:length(model$reacID)]
	bestSW = est$xbest[(length(model$reacID)+1):length(est$xbest)]
	strongWeak[negEdges] = bestSW

	simBest = simulatorDelay(CNOlist, model, simList, indexList,
	boolUpdates=boolUpdates, delayThresh=bestDelay, strongWeak=strongWeak)
	simBest = convert2array(simBest, dim(CNOlist@signals[[1]])[1], length(model$namesSpecies), boolUpdates)
	simBest = simBest[,indexList$signals,]
#	plotCNOlist(plotData(CNOlistPB, simBest))
	
	diff <- (simBest - ySilico)
    r <- diff^2
    deviationPen <- sum(r[!is.na(r)])/nTimes
    NAPen <- NAFac * length(which(is.na(simBest)))
    nDataPts <- dim(CNOlist@signals[[1]])[1] * dim(CNOlist@signals[[1]])[2] * nTimes
    nInputs <- length(which(model$interMat == -1))
    
    # nInTot: number of inputs of expanded model nInputs: number of inputs of cut model
    sizePen <- (nDataPts * sizeFac * nInputs)/nInTot
    
    score <- deviationPen + NAPen + sizePen
	return(list(score=score, estimate=est$xbest, xCoords=xCoords, yInter=ySilico, simResults=simBest))

}