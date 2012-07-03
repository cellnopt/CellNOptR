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

getFitTimeScale <- function(simList, CNOlist, model, indexList, sizeFac=0.0001, NApenFac=1, boolUpdates, timeSplit="early", divTime=NULL, simResultsT1, lowerB=lowerB, upperB=upperB) {

	if(is.null(divTime)) {divTime = CNOlist$timeSignals[length(CNOlist$timeSignals)]}
	if(timeSplit=="early") {
		timeExper = CNOlist$timeSignals[CNOlist$timeSignals <= divTime]
		timeIndex = which(CNOlist$timeSignals <= divTime)
		boolUpdates = boolUpdates[1]
	} else if(timeSplit=="late") {
		timeExper = CNOlist$timeSignals[CNOlist$timeSignals > divTime]
		timeIndex = which(CNOlist$timeSignals > divTime)
		boolUpdates = boolUpdates[2]
	}	

	# make sure timeExper[1] = 0 for fitting
	timeExper = timeExper - timeExper[1]
	
	# dimensions, time points
	times = length(timeExper)
	sigs = dim(CNOlist$valueSignals[[1]])

	# simulator
	if(timeSplit=="early") {
		yBool = simulatorTimeScale(CNOlist, model, simList, indexList, boolUpdates)
	} else if(timeSplit=="late") {
		yBool = simulatorTimeScaleT2(simResultsT1[,,dim(simResultsT1)[3]], CNOlist, model, simList, indexList, boolUpdates)
	}
	
	yBool = yBool[,indexList$signals,]

	##### IN LOOP #####

	splineStore = list()
	splineAdd = 1

	for (nExper in 1:dim(CNOlist$valueSignals[[1]])[1]) {
		for (nSig in 1:dim(CNOlist$valueSignals[[1]])[2]) {
			yTest = c()
			for (a in timeIndex) {
				yTest =c(yTest, CNOlist$valueSignals[[a]][nExper, nSig])
			}

			if(!is.na(yTest[1])) {
				cS = splinefun(timeExper, yTest)
				splineStore[splineAdd] = list(cS)
			} else {
				cS = splinefun(timeExper,rep(0,times))
				splineStore[splineAdd] = list(cS)
			}
		splineAdd = splineAdd + 1
		}
	}

	# optimization
	findTimeScale <- function(yB, splines) {

		# what to optimize
		taufinder <- function(whatScale) {
	
			ySilico = array(dim=dim(yB))
			numberPoints = dim(yB)[3]
			xCoords = seq(0,by=whatScale,lengthOut=numberPoints)
   			count.1 = 1
    
    		for(nExper in 1:dim(yB)[1]) {
    			for(nSig in 1:dim(yB)[2]) { 
                
            		yOut = splines[[count.1]](xCoords)
            		ySilico[nExper,nSig,] = yOut;
            		count.1 = count.1 + 1
   				}
  			}
  	
  			ErrorVector = as.vector(ySilico) - as.vector(yB)
    		sse = sum(ErrorVector^2)
    		return(sse)	
		}

		seed1 = 0.99
		est1 = optim(seed1, taufinder, method="L-BFGS-B", lower=lowerB, upper=upperB)
	}

	myEstimate = findTimeScale(yBool, splineStore)
	yFinal = array(dim = dim(yBool))
	xCoords = seq(0,by=myEstimate$par,lengthOut=boolUpdates)
	
	count2 = 1
	for(nExper in 1:dim(yBool)[1]) {
		for(nSig in 1:dim(yBool)[2]) { 
                
			yOut = splineStore[[count2]](xCoords)
     		yFinal[nExper,nSig,] = yOut;
      		count2 = count2 + 1
   		}
	}

	diff <- (yBool - yFinal) # * my.weights
	r <- diff^2
	deviationPen <- sum(r[!is.na(r)])
	NApen <- NApenFac * length(which(is.na(yBool)))
	dims = dim(yFinal)
	nDataPts <- dims[1] * dims[2] * dims[3]
	nReac <- length(model$reacID)
	nInputs <- length(which(model$interMat == -1))
	sizePen <- (nDataPts * sizeFac * nInputs) / nReac
	score <- deviationPen + NApen + sizePen
	return(list(score=score, estimate=myEstimate$par, xCoords=xCoords, yInter=yFinal, yBool=yBool))

}
