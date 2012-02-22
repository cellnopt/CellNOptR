# delay and state finder
# 20/10/2011

getFitPause <- function(SimList, CNOlist, Model, indexList, sizeFac=0.0001, NAPenFac=1, boolUpdates, timeSplit="early", divTime=NULL, SimResultsT1) {

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
	
	# dimensions, time.points
	times = length(time.exper)
	sigs = dim(CNOlist$valueSignals[[1]])

	################################################################################
	
	# fit the curves to the experimental data for each cell
	# uses splinefun

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
	
	################################################################################
	
	# here pick out the negative feedback that can be optimized
	# as either strong/weak with/without delay

	what.loops = feedbackFinder1(Model)
	neg.reacs = c()
	for(b in 1:length(what.loops)) {
		loop1 = what.loops[[b]]
		loop1 = loop1[loop1 > 0]
		loop1 = c(loop1, loop1[1])
		for(a in 1:length(loop1)-1) {
			lhs = loop1[a]
			rhs = loop1[a+1]
			lhs.reac = which(Model$interMat[lhs,] == -1)
			rhs.reac = which(Model$interMat[rhs,] == 1)
			reac = intersect(lhs.reac, rhs.reac)
			if(any(Model$notMat[,reac] == 1)) {
				neg.reacs = c(neg.reacs, reac)	
			}
		}
	}
	neg.reacs = unique(neg.reacs)
	negEdges = rep(0,length(Model$reacID))
	negEdges[neg.reacs] = 1
	strongWeak = rep(0,length(negEdges[negEdges == 1]))
	
	################################################################################

	# optimization

	findDelayState <- function(optimString, splines) {

		# what to optimize
		delayStateFinder <- function(optimString) {
		
			what.scale = optimString[1]
			delayThresh = optimString[2:(1+length(optimString))]
		#	strongWeak = optimString[(2+length(Model$reacID)):length(optimString)]
			yB = simulatorPause(CNOlist, Model, SimList, indexList, boolUpdates, delayThresh)
			yB = yB[,indexList$signals,]
		
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

		len.delay = length(Model$reacID)
		len.strongWeak = length(optimString[(2+length(Model$reacID)):length(optimString)])
		problem = list(f=delayStateFinder, x_L=c(0.5, rep(0,len.delay)), x_U=c(1.5, rep(boolUpdates,len.delay)), int_var=length(Model$reacID))
		opts = list(maxtime=240, local_solver=0, dim_refset=6, ndiverse=50)
		est.1 = essR(problem, opts)	
	}

	################################################################################

	optimString = c(1, rep(1,length(Model$reacID)))
	my.estimate = findDelayState(optimString=optimString, splines=spline.store)
	bestTimeStep = my.estimate$xbest[1]
	bestDelays = my.estimate$xbest[2:length(my.estimate$xbest)]

	yBool = simulatorTotal3(CNOlist, Model, SimList, indexList, boolUpdates, bestDelays)
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