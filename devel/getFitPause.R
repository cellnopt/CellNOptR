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
		delayThresh = # integer section of optimString
		strongWeak[negEdges] = # binary section of optimString
		
		simResults = simulatorTest(CNOlist, model, simList, indexList,
		boolUpdates=boolUpdates, delayThresh=delayThresh, strongWeak=strongWeak)
		simResults = convert2array(simResults, dim(CNOlist@signals[[1]])[1], length(model$namesSpecies), boolUpdates)
		simResults = simResults[,indexList$signals,]

		ySilico = array(dim=dim(simResults))
		numberPoints = dim(simResults)[3]
		xCoords = seq(0, by=whatScale, length.out=numberPoints)
   		count = 1
    
    	for(nExper in 1:dim(simResults)[1]) {
    		for(nSig in 1:dim(simResults)[2]) {
    				             
            	yOut = splines[[count]](xCoords)
            	ySilico[nExper,nSig,] = yOut;
            	count = count + 1
   			}
  		}
  	
  		errorVector = as.vector(ySilico) - as.vector(simResults)
    	sse = sum(errorVector^2)
    	return(sse)			
	}

	lenDelay = length(Model$reacID)
	len.strongWeak = length(optimString[(2+length(Model$reacID)):length(optimString)])
	problem = list(f=delayStateFinder, x_L=c(0.1, rep(0,len.delay)), x_U=c(9.9, rep(boolUpdates,len.delay)), int_var=length(Model$reacID))
	opts = list(maxtime=240, local_solver=0, dim_refset=6, ndiverse=50)
	est.1 = essR(problem, opts)	

	################################################################################

	optimString = c(1, rep(1,length(Model$reacID)))
	my.estimate = findDelayState(optimString=optimString, splines=spline.store)
	bestTimeStep = my.estimate$xbest[1]
	bestDelays = my.estimate$xbest[2:length(my.estimate$xbest)]

	yBool = simulatorPause(CNOlist, Model, SimList, indexList, boolUpdates, bestDelays)
	yBool = yBool[,indexList$signals,]
	yFinal = array(dim = dim(yBool))
	xCoords = seq(0, by=bestTimeStep, length.out=boolUpdates)
	count.2 = 1

	for(nExper in 1:dim(yBool)[1]) {
		for(nSig in 1:dim(yBool)[2]) { 
                
			yOut = spline.store[[count.2]](xCoords)
     		yFinal[nExper,nSig,] = yOut;
      		count.2 = count.2 + 1
   		}
	}

	Diff <- (yBool - yFinal)
	r <- Diff^2
	deviationPen <- sum(r[!is.na(r)])
	NAPen <- NAPenFac * length(which(is.na(yBool)))
	dims = dim(yFinal)
	nDataPts <- dims[1] * dims[2] * dims[3]
	nReac <- length(Model$reacID)
	nInputs <- length(which(Model$interMat == -1))
	sizePen <- (nDataPts * sizeFac * nInputs) / nReac
	score <- deviationPen + NAPen + sizePen
	return(list(score=score, estimate=bestTimeStep, bestDelays=bestDelays, xCoords=xCoords, yInter=yFinal, yBool= yBool))

}