# boolsync simulator
simulatorTimeScale <- function(CNOlist, Model, SimList, indexList, boolUpdates) {

	nSp <- dim(Model$interMat)[1]
	nReacs <- dim(Model$interMat)[2]	
	nCond <- dim(CNOlist$valueStimuli)[1]
	yBool = array(dim=c(nCond, nSp, boolUpdates))
	yBool[,,1] = 0

	if(is.null(dim(Model$interMat))) { 
		nSp <- length(Model$interMat)
		nReacs <- 1
	}

	# this holds, for each sp, how many reactions have that sp as output
	# + other function definitions here
	endIx <- rep(NA,nSp)
		for(i in 1:nSp){
			endIx[i] <- length(which(SimList$maxIx == i))
		}
	
	filltempCube <- function(x) {
		cMatrix <- matrix(data=x, nrow=nReacs, ncol=nCond)
		cVector <- apply(cMatrix, 1, function(x){return(x)})
		return(cVector)
	}
		
	minNA <- function(x) {
		if(all(is.na(x))) {
			return(NA)
		} else {
			return(min(x, na.rm=TRUE))
		}
	}
				
	compOR <- function(x) {
		if(all(is.na(x[which(SimList$maxIx == s)]))) {
			res <- NA
		} else {
			res <- max(x[which(SimList$maxIx == s)], na.rm=TRUE)
		}
		return(res)
	}
	
	# this value is used to test the stop condition for difference between 2 iterations
	testVal <- 1E-3

	# create an initial values matrix	
	initValues <- matrix(data=NA, nrow=nCond, ncol=nSp)
	colnames(initValues) <- Model$namesSpecies

	# set the initial values of the stimuli
	initValues[,indexList$stimulated] <- CNOlist$valueStimuli

	# flip the inhibitors so that 0 = inhibited / 1 = noninhibited
	valueInhibitors <- 1-CNOlist$valueInhibitors
	valueInhibitors[which(valueInhibitors == 1)] <- NA

	# set the initial values of the inhibited species: 0 if inhibited, untouched if not inhibited
	initValues[,indexList$inhibited] <- valueInhibitors

	# initialise main loop
	newInput <- initValues

	# main loop

	for(count in 2:boolUpdates) {
		
		outputPrev <- newInput
		# this is now a 2 column matrix that has a column for each input (column in finalCube)
		# and a set of rows for each reac (where a set contains as many rows as conditions)
		# all concatenated into one long column
				
		if(nReacs > 1) {
			tempStore <- apply(SimList$finalCube, 2, function(x){return(outputPrev[,x])})
			tempIxNeg <- apply(SimList$ixNeg, 2, filltempCube)
			tempIgnore <- apply(SimList$ignoreCube, 2, filltempCube)
		} else {
			tempStore <- outputPrev[,SimList$finalCube]
			tempIxNeg <- matrix(SimList$ixNeg, nrow=nCond, ncol=length(SimList$ixNeg), byrow=TRUE)
			tempIgnore <- matrix(SimList$ignoreCube, nrow=nCond, ncol=length(SimList$ignoreCube), byrow=TRUE)
		}

		# set to NA the values that are "dummies", so they won't influence the min
		tempStore[tempIgnore] <- NA

		# Flip the values that enter with a negative sign
		tempStore[tempIxNeg] <- 1-tempStore[tempIxNeg]
	
		# compute all the ands by taking, for each gate, the min value across the inputs of that gate
		if(nReacs > 1) {
			
			outputCube <- apply(tempStore, 1, minNA)

			# outputCube is now a vector of length (nCond*nReacs) that contains the input of each reaction in
			# each condition, concatenated as such allcond4reac1, allcond4reac2, etc...			# this is transformed into a matrix with a column for each reac and a row for each cond		
			outputCube <- matrix(outputCube, nrow=nCond, ncol=nReacs)
			# go through each species, and if it has inputs, then take the max across those input reactions
			# i.e. compute the ORs
		
			for(s in 1:nSp){
				if(endIx[s] != 0){
					newInput[,s] <- apply(outputCube, 1, compOR)
				}
			}
		
		} else {
			outputCube <- ifelse(all(is.na(tempStore)), NA, min(tempStore,na.rm=TRUE))
			newInput[,SimList$maxIx] <- outputCube
		}
	
		# reset the inhibitors and stimuli
		for(stim in 1:length(indexList$stimulated)) {
			stimM <- cbind(CNOlist$valueStimuli[,stim], newInput[,indexList$stimulated[stim]])
			maxNA <- function(x) {
				return(max(x, na.rm=TRUE))
			}
			stimV <- apply(stimM, 1, maxNA)
			newInput[,indexList$stimulated[stim]] <- stimV
		}
		
		valueInhibitors <- 1-CNOlist$valueInhibitors
		newInput[,indexList$inhibited] <- valueInhibitors * newInput[,indexList$inhibited]
		# replace NAs with zeros to avoid having the NA penalty applying to unconnected species
		newInput[is.na(newInput)] <- 0
		outputPrev[is.na(outputPrev)] <- 0
		yBool[,,count] = newInput		
	}

	return(yBool)
}

