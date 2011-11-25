
getFitTimeScale <- function(SimList, CNOlist, Model, indexList, sizeFac=0.0001, NAPenFac=1, boolUpdates, timeSplit="early", divTime=NULL, SimResultsT1) {

	if(is.null(divTime)) {divTime = CNOlist$timeSignals[length(CNOlist$timeSignals)]}
	if(timeSplit=="early") {
		time.exper = CNOlist$timeSignals[CNOlist$timeSignals <= divTime]
		time.index = which(CNOlist$timeSignals <= divTime)
		boolUpdates = boolUpdates[1]
	} else if(timeSplit=="late") {
		time.exper = CNOlist$timeSignals[CNOlist$timeSignals > divTime]
		time.index = which(CNOlist$timeSignals > divTime)
		boolUpdates = boolUpdates[2]
	}	

	# make sure time.exper[1] = 0 for fitting
	time.exper = time.exper - time.exper[1]
	
	# dimensions, time points
	times = length(time.exper)
	sigs = dim(CNOlist$valueSignals[[1]])

	# simulator
	if(timeSplit=="early") {
		yBool = simulatorTimeScale(CNOlist, Model, SimList, indexList, boolUpdates)
	} else if(timeSplit=="late") {
		yBool = simulatorTimeScaleT2(SimResultsT1[,,dim(SimResultsT1)[3]], CNOlist, Model, SimList, indexList, boolUpdates)
	}
	
	yBool = yBool[,indexList$signals,]

	##### IN LOOP #####

	spline.store = list()
	spline.add = 1

	for (nExper in 1:dim(CNOlist$valueSignals[[1]])[1]) {
		for (nSig in 1:dim(CNOlist$valueSignals[[1]])[2]) {
			yTest = c()
			for (a in time.index) {
				yTest =c(yTest, CNOlist$valueSignals[[a]][nExper, nSig])
			}

			if(!is.na(yTest[1])) {
				cs = splinefun(time.exper, yTest)
				spline.store[spline.add] = list(cs)
			} else {
				cs = splinefun(time.exper,rep(0,times))
				spline.store[spline.add] = list(cs)
			}
		spline.add = spline.add + 1
		}
	}

	# optimization
	findTimeScale <- function(yB, splines) {

		# what to optimize
		taufinder <- function(what.scale) {
	
			ySilico = array(dim=dim(yB))
			number.points = dim(yB)[3]
			x.coords = seq(0,by=what.scale,length.out=number.points)
   			count.1 = 1
    
    		for(nExper in 1:dim(yB)[1]) {
    			for(nSig in 1:dim(yB)[2]) { 
                
            		yOut = splines[[count.1]](x.coords)
            		ySilico[nExper,nSig,] = yOut;
            		count.1 = count.1 + 1
   				}
  			}
  	
  			ErrorVector = as.vector(ySilico) - as.vector(yB)
    		sse = sum(ErrorVector^2)
    		return(sse)	
		}

		seed.1 = 0.99
		est.1 = optim(seed.1, taufinder, method="L-BFGS-B", lower=0.5, upper=1.5)
	}

	my.estimate = findTimeScale(yBool, spline.store)
	yFinal = array(dim = dim(yBool))
	xCoords = seq(0,by=my.estimate$par,length.out=boolUpdates)
	
	count.2 = 1
	for(nExper in 1:dim(yBool)[1]) {
		for(nSig in 1:dim(yBool)[2]) { 
                
			yOut = spline.store[[count.2]](xCoords)
     		yFinal[nExper,nSig,] = yOut;
      		count.2 = count.2 + 1
   		}
	}

	Diff <- (yBool - yFinal) # * my.weights
	r <- Diff^2
	deviationPen <- sum(r[!is.na(r)])
	NAPen <- NAPenFac * length(which(is.na(yBool)))
	dims = dim(yFinal)
	nDataPts <- dims[1] * dims[2] * dims[3]
	nReac <- length(Model$reacID)
	nInputs <- length(which(Model$interMat == -1))
	sizePen <- (nDataPts * sizeFac * nInputs) / nReac
	score <- deviationPen + NAPen + sizePen
	return(list(score=score, estimate=my.estimate$par, xCoords=xCoords, yInter=yFinal, yBool=yBool))

}