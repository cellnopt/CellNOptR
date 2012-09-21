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

simulatorTimeScaleT2 <- function(simResultsT1, CNOlist, model, simList, indexList, boolUpdates) {
		
	nSp <- dim(model$interMat)[1]
	nReacs <- dim(model$interMat)[2]
	nCond <- dim(CNOlist$valueStimuli)[1]
	yBoolT2 = array(dim=c(nCond, nSp, boolUpdates))
	
	if(is.null(dim(model$interMat))) { 
		nSp <- length(model$interMat)
		nReacs <- 1
	}
	
	# this holds, for each sp, how many reactions have that sp as output
	# + other functions here
	
	endIx <- rep(NA, nSp)
	for(i in 1:nSp) {
		endIx[i] <- length(which(simList$maxIx == i))
	}
		
	fillTempCube <- function(x) {
		cMatrix <- matrix(data=x, nrow=dim(simList$ixNeg)[1], ncol=nCond)
		cVector <- apply(cMatrix, 1, function(x){return(x)})
		return(cVector)
	}
	
	minNA <- function(x) {
		if(all(is.na(x))) {
			return(NA)
		} else {
			return(min(x,na.rm=TRUE))
		}
	}
	
	compOr <- function(x) {
		if(all(is.na(x[which(simList$maxIx == s)]))){
			res <- NA
		} else {
			res <- max(x[which(simList$maxIx == s)], na.rm=TRUE)
		}
		return(res)
	}
	
	maxNA <- function(x) {
		return(max(x, na.rm=TRUE))
	}
	
	# create an initial values matrix	
	initValues <- simResultsT1

	# initialise main loop
	newInput <- initValues

	# first iteration
	outputPrev <- newInput
	tempStore <- apply(simList$finalCube, 2, function(x){return(outputPrev[,x])})
	tempIxNeg <- apply(simList$ixNeg, 2, fillTempCube)
	tempIgnore <- apply(simList$ignoreCube, 2, fillTempCube)
	tempStore[tempIgnore] <- NA
	tempStore[tempIxNeg] <- 1-tempStore[tempIxNeg]

	outputCube <- apply(tempStore, 1, minNA)
	outputCube <- matrix(outputCube, nrow=nCond, ncol=nReacs)
	
	# rewrite anything that comes to a node that also receives a t2 branch,ie set it to the same 
	# as our t2 reac for those reacs, so that they won't influence the OR
	reacsT2 <- which(model$times == 2)
	if(length(reacsT2) > 0) {
		for(i in reacsT2) {
			outNode <- which(model$interMat[,i] > 0)
			reacs2Overwrite <- which(model$interMat[outNode,] > 0)
			if(length(reacs2Overwrite) != 0) {
				for(n in 1:length(reacs2Overwrite)) {
					outputCube[,reacs2Overwrite[n]] <- outputCube[,i]
				}
			}	
		}
	}	
	
	for(s in 1:nSp) {
		if(endIx[s] != 0) {
			newInput[,s] <- apply(outputCube, 1, compOr)
		}
	}	
	
	for(stim in 1:length(indexList$stimulated)) {
		stimM <- cbind(CNOlist$valueStimuli[,stim], newInput[,indexList$stimulated[stim]])
		stimV <- apply(stimM, 1, maxNA)
		newInput[,indexList$stimulated[stim]] <- stimV
	}
	
	valueInhibitors <- 1-CNOlist$valueInhibitors
	newInput[,indexList$inhibited] <- valueInhibitors*newInput[,indexList$inhibited]
	newInput[is.na(newInput)] <- 0
	outputPrev[is.na(outputPrev)] <- 0
	firstIter <- newInput
	yBoolT2[,,1] = firstIter

	##### main loop #####
	
	for(countLoop in 2:boolUpdates) {
		
		outputPrev <- newInput
		
		# this is now a 2 columns matrix that has a column for each input (column in finalCube)
		# and a set of rows for each reac (where a set contains as many rows as conditions)
		# all concatenated into one long column
		tempStore <- apply(simList$finalCube, 2, function(x){return(outputPrev[,x])})
		tempIxNeg <- apply(simList$ixNeg, 2, fillTempCube)
		tempIgnore <- apply(simList$ignoreCube, 2, fillTempCube)
		
		# set to NA the values that are "dummies", so they won't influence the min
		tempStore[tempIgnore] <- NA

		# flip the values that enter with a negative sign
		tempStore[tempIxNeg] <- 1-tempStore[tempIxNeg]
		
		# compute all the ands by taking, for each gate, the min value across the inputs of that gate
		outputCube <- apply(tempStore, 1, minNA)

		# outputCube is now a vector of length (nCond*nReacs) that contains the input of each reaction in
		# each condition, concatenated as such allcond4reac1,allcond4reac2,etc..
		# this is transformed into a matrix with a column for each reac and a row for each cond
		outputCube <- matrix(outputCube, nrow=nCond,ncol=nReacs)

		# go through each species, and if it has inputs, then take the max across those input reactions
		# i.e. compute the ORs
		for(s in 1:nSp) {
			if(endIx[s] != 0) {
				compOr <- function(x) {
					if(all(is.na(x[which(simList$maxIx == s)]))) {
						res <- NA
					} else {
						res <- max(x[which(simList$maxIx == s)], na.rm=TRUE)
					}
					return(res)
				}
				newInput[,s] <- apply(outputCube, 1, compOr)
			}
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
		newInput[,indexList$inhibited] <- valueInhibitors*newInput[,indexList$inhibited]
		
		# set all the nodes that are targets of a t2 reaction to the state that they had at the first iteration
		if(length(reacsT2) != 0) {
			t2reacs <- model$interMat[,reacsT2]
			for(r in 1:length(reacsT2)) {
				if(length(reacsT2) == 1) {
					target <- which(t2reacs > 0)
				} else {
					target <- which(t2reacs[,r] > 0)
				}
				newInput[,target] = firstIter[,target]
			}
		}

		# replace NAs with zeros to avoid having the NA penalty applying to unconnected species
		newInput[is.na(newInput)] <- 0
		outputPrev[is.na(outputPrev)] <- 0
		yBoolT2[,,countLoop] = newInput
	}
	
	return(yBoolT2)
}

