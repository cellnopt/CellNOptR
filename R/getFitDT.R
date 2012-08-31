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
# $Id$

getFitDT <- function(
  simResults,
  CNOlist,
  model,
  indexList,
  timePoint=c("t1","t2"),
  sizeFac=0.0001,
  NAFac=1,
  nInTot,
  simResultsT0=NA,
# additional arguments over getFit
  boolUpdates = 20,
  upperB=,
  lowerB=,

) {

# ignore for the moment and return to for TN
# 	if(is.null(divTime)) {divTime = CNOlist$timeSignals[length(CNOlist$timeSignals)]}
# 	if(timeSplit=="early") {
# 		timeExper = CNOlist$timeSignals[CNOlist$timeSignals <= divTime]
# 		timeIndex = which(CNOlist$timeSignals <= divTime)
# 		boolUpdates = boolUpdates[1]
# 	} else if(timeSplit=="late") {
# 		timeExper = CNOlist$timeSignals[CNOlist$timeSignals > divTime]
# 		timeIndex = which(CNOlist$timeSignals > divTime)
# 		boolUpdates = boolUpdates[2]
# 	}	
  
    if ((class(CNOlist)=="CNOlist")==FALSE){
      CNOlist = CellNOptR::CNOlist(CNOlist)
    }
	# make sure timeExper[1] = 0 for fitting
	timeExper = CNOlist$timeSignals - CNOlist$timeSignals[1]
	
	# dimensions, time points
	times = length(timeExper)
	sigs = dim(CNOlist$valueSignals[[1]])

	# simulator will need to add simulatorDT-TN here, ignore for now
# 	if(timeSplit=="early") {
# 		yBool = simulatorTimeScale(CNOlist, model, simList, indexList, boolUpdates)
# 	} else if(timeSplit=="late") {
# 		yBool = simulatorTimeScaleT2(simResultsT1[,,dim(simResultsT1)[3]], CNOlist, model, simList, indexList, boolUpdates)
# 	}
	# this is simResults
  #yBool = simulatorTimeScale(CNOlist, model, simList, indexList, boolUpdates)
	#yBool = yBool[,indexList$signals,]

	simResults<-simResults[,indexList$signals,]
  # combine simResults with T0?
  # add T0 only if interested in first pseudo steady state
  # simResultsT0 has been passed
  
	Diff0 <- simResultsT0[,indexList$signals]-CNOlist@signals[[1]]
	
  
  
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
			xCoords = seq(0,by=whatScale,length.out=numberPoints)
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
	xCoords = seq(0,by=myEstimate$par,length.out=boolUpdates)
	
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
