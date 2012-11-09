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

getFitPause <- function(CNOlist, model, indexList, sizeFac = 1e-04, NAFac = 1, nInTot, boolUpdates,  
    lowerB, upperB) {
    	    
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
	
	
	################################################################################
	
	# TODO: neg data should be index list from model$reacID
	# of what strongWeak edges there are
	strongWeak = rep(0,length(model$reacID))
	negEdges = feedbackWrapper(model) # the index of edges that can be changed
	# optimization

	# what to optimize
	delayFinder <- function(optimString) {
		
		whatScale = optimString[1] # continuous with upper/lower bounds
		delayThresh = optimString[2:(length(model$reacID)+1)] # integer section of optimString
		strongWeak = c(0,0,0,0,1,0,0,0,0,0,0)
	#	strongWeak = optimString[(length(model$reacID)+2):length(optimString)] # binary section of optimString
		
		simResults = simulatorTest(CNOlist, model, simList, indexList,
		boolUpdates=boolUpdates, delayThresh=delayThresh, strongWeak=strongWeak)
		simResults = convert2array(simResults, dim(CNOlist@signals[[1]])[1], length(model$namesSpecies), boolUpdates)
		simResults = simResults[,indexList$signals,]

		ySilico = array(dim=dim(simResults))
		xCoords = seq(0, by=whatScale, length.out=dim(simResults)[3])
   		count = 1
    
    	for(nExper in 1:dim(simResults)[1]) {
    		for(nSig in 1:dim(simResults)[2]) {
    				             
            	yOut = splineStore[[count]](xCoords)
            	ySilico[nExper,nSig,] = yOut;
            	count = count+1
   			}
  		}
  	
  		errorVector = as.vector(ySilico) - as.vector(simResults)
    	sse = sum(errorVector^2)
    	return(sse)			
	}

	problem = list(
		f=delayFinder,
		x_L=c(0.9, rep(0,length(model$reacID))),
		x_U=c(9.9, rep(10,length(model$reacID))),
		int_var=length(model$reacID)
	#	bin_var=length(model$reacID)
	)
	opts = list(maxtime=15, local_solver=0, dim_refset=6, ndiverse=50)
	optimString = c(1, rep(1,length(model$reacID)))
	est = essR(problem, opts)	

	################################################################################

	bestScale = est$xbest[1]
	bestDelay = est$xbest[2:length(est$xbest)]

	simBest = simulatorTest(CNOlist, model, simList, indexList,
	boolUpdates=boolUpdates, delayThresh=bestDelay, strongWeak=strongWeak)
	simBest = convert2array(simBest, dim(CNOlist@signals[[1]])[1], length(model$namesSpecies), boolUpdates)
	simBest = simBest[,indexList$signals,]
	
#	plotCNOlist(plotData(CNOlistPB, simBest))


	yInter = array(dim = dim(simBest))
	xCoords = seq(0, by=bestTimeStep, length.out=boolUpdates)

	count = 1
	for(nExper in 1:dim(yInter)[1]) {
		for(nSig in 1:dim(yInter)[2]) { 
                
			yOut = spline.store[[count]](xCoords)
     		yFinal[nExper,nSig,] = yOut;
      		count = count+1
   		}
	}

    diff <- (simBest - yInter)
    r <- diff^2
    deviationPen <- sum(r[!is.na(r)])/nTimes
    NAPen <- NAFac * length(which(is.na(simResults)))
    nDataPts <- dim(CNOlist@signals[[1]])[1] * dim(CNOlist@signals[[1]])[2] * nTimes
    nInputs <- length(which(model$interMat == -1))
    
    # nInTot: number of inputs of expanded model nInputs: number of inputs of cut model
    sizePen <- (nDataPts * sizeFac * nInputs)/nInTot
    
    score <- deviationPen + NAPen + sizePen
	return(list(score=score, estimate=bestScale, bestDelay=bestDelay, xCoords=xCoords, yInter=yInter, simResults=simBestl))

}