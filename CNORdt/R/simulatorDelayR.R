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

simulatorDelayR <- function(CNOlist, model, simList, indexList, boolUpdates, delayThresh, strongWeak) {

	nSp <- dim(model$interMat)[1] # number of species in model
	nReacs <- dim(model$interMat)[2] # number of reactions
	nCond <- dim(CNOlist$valueStimuli)[1] # number of conditions

	yBool = array(dim=c(nCond, nSp, boolUpdates)) # an array to store the model simulation
	colnames(yBool) = model$namesSpecies
	yBool[,,1] = 0 # initial conditions for simulation

  # if the model has only 1 reaction
	if(is.null(dim(model$interMat))) { 
		nSp <- length(model$interMat)
		nReacs <- 1
	}
		
	# for each species, how many times is it an output
  endIx <- rep(NA,nSp)
	for(i in 1:nSp){
		endIx[i] <- length(which(simList$maxIx == i))
	}
		
  # functions 
	compOR <- function(x){
		if(all(is.na(x[which(simList$maxIx == s)]))){
			res <- NA
		} else {
			res <- max(x[which(simList$maxIx == s)], na.rm=TRUE)
		}
			return(res)
	}
	
	minNA <- function(x){
		if(all(is.na(x))){
			return(NA)
		} else {
			return(min(x, na.rm=TRUE))
		}
	}
	
	maxNA <- function(x){
		return(max(x, na.rm=TRUE))
	}
	
	filltempCube <- function(x) {
		cMatrix <- matrix(data=x, nrow=nReacs, ncol=nCond)
		cVector <- apply(cMatrix, 1, function(x){return(x)})
		return(cVector)
	}
	
	# create an initial values matrix	
	initValues <- matrix(data=NA, nrow=nCond, ncol=nSp)
	colnames(initValues) <- model$namesSpecies

	# set the initial values of the stimuli
	initValues[,indexList$stimulated] <- CNOlist$valueStimuli

	# flip the inhibitors so that 0 = inhibited / 1 = noninhibited
	valueInhibitors <- 1-CNOlist$valueInhibitors
	valueInhibitors[which(valueInhibitors == 1)] <- NA
	
  # set the initial values of the inhibited species: 0 if inhibited, untouched if not inhibited
	initValues[,indexList$inhibited] <- valueInhibitors
	
	# set everything else = 0 (necessary for time-course data)
	initValues[is.na(initValues)] = NA

	# initialise main loop
	newInput <- initValues

	# time delay
	delayThreshTot = rep(delayThresh, 1, each=nCond) # expand the delay vector across all conditions
	strongWeakTot = rep(strongWeak, 1, each=nCond) # expand the strongWeak vector across all conditions

	delayCount = which(delayThreshTot > 0) # index of above where delay is present
	strongWeakCount = which(strongWeakTot==1)	# index of above where a strong edge is present
	
	# make a matrix to store outputCubes
	allCubes = matrix(NA, nrow=nReacs*nCond, ncol=boolUpdates)
	rownames(allCubes) = rep(model$reacID,1,each=nCond)
	allCubes[,1]=0
    
	# main loop
	for(countBool in 2:boolUpdates) {
		
		outputPrev <- newInput # values of species for each condition copied to t-1
				
		if(nReacs > 1) {
			# tempStore is a matrix of ncol = the maximum inputs across all reactions
      # and nrow = reactions repeated by number of conditions
      tempStore <- apply(simList$finalCube, 2, function(x){return(outputPrev[,x])})
			tempIxNeg <- apply(simList$ixNeg, 2, filltempCube)
			tempIgnore <- apply(simList$ignoreCube, 2, filltempCube)
		} else {
			tempStore <- outputPrev[,simList$finalCube]
			tempIxNeg <- matrix(simList$ixNeg, nrow=nCond, ncol=length(simList$ixNeg), byrow=TRUE)
			tempIgnore <- matrix(simList$ignoreCube, nrow=nCond, ncol=length(simList$ignoreCube), byrow=TRUE)
		}

		tempStore[tempIgnore] <- NA # update tempStore by tempIgnore e.g. where max. inputs = 2 across all reactions
    # and there is only 1 input in a particular reaction

		# flip the values that enter with a negative sign
		tempStore[tempIxNeg] <- 1-tempStore[tempIxNeg]
	    
		# compute all the ands by taking, for each gate, the min value across the inputs of that gate
		if(nReacs > 1) {
			
			outputCube <- apply(tempStore, 1, minNA) # find the minimum per row - this is the output of the AND gates
      # if the reaction is not an AND gate, it will just take the single value
      # as the other value(s) will be NA as a result of the 'tempStore[tempIgnore]' command above
			      
      # delays
			outputOn = which(!is.na(outputCube)) # which reactions (across conditions) have produced an output
			delayOn = intersect(outputOn, delayCount) # of these reactions, is there a delay?
			strongWeakOn = intersect(outputOn, strongWeakCount) # or, is this a strong reaction?

			if(length(delayOn)) { # if there is a delay on activated edges
				for(a in delayOn) { # for each delay
					toAdd = countBool:(countBool+delayThreshTot[a]) # the indices in allCubes to edit
					if(toAdd[length(toAdd)] > boolUpdates) {
						toAdd = countBool:boolUpdates # truncate this indices if it goes past boolUpdates
					}				
					futureData = c(rep(0.5,delayThreshTot[a]), outputCube[a])[1:length(toAdd)] # add a delay 0.5 and then the output
					futureCheck = which(is.na(allCubes[a, toAdd])) # make sure it doesn't overlap with past edits
					if(length(futureCheck)) {
						allCubes[a, toAdd][futureCheck] = futureData[futureCheck]
					}
				}
			}
			
      # where allCubes has not been edited, add outputCube to current column
			as.normal = which(is.na(allCubes[,countBool]))
			allCubes[as.normal,countBool] = outputCube[as.normal]
    
			# check for strong edges
			if(length(strongWeakOn)) { # if there are strong edges...
				for(b in strongWeakOn) { # for each 'stronge edge' index
					if(!is.na(allCubes[b,countBool])) { # if this entry for this time is not na
					  
            whatReac = ceiling(b/nCond) # what reactions are involved
					  whatCond = (b %% nCond) # what conditions are involved
            if(!whatCond) whatCond = nCond
					  output1 = which(model$interMat[,whatReac] > 0) # what are the output species
					  reacsToFreeze = which(model$interMat[output1,] > 0) # what other reactions have these output species
					  range1 = countBool:boolUpdates # what is the (time) range to override?
            bPlus = (reacsToFreeze-1)*nCond + whatCond # what other interactions are affected
					  bPlus = setdiff(bPlus,b)
            
            print(paste("Round:",countBool))
            print(paste("Reaction:",whatReac))
					  print(paste("Condition:",whatCond))
            
            toChange = which(!allCubes[b,range1]==0.5 | is.na(allCubes[b,range1])) # don't overwrite the delay
            allCubes[b,range1[toChange]] = outputCube[b]
				    if(length(bPlus)) { # change the other interactions as well
				      allCubes[bPlus,range1[toChange]] = outputCube[b]
				    }  
            
						strongWeakCount = setdiff(strongWeakCount, b) # remove 'b' from list - stays 'strong' for whole simulation
					}
				}	
			}
		
			# reform outputCube (the same as current (time) column of allCubes as a matrix
			outputCube <- matrix(allCubes[,countBool], nrow=nCond, ncol=nReacs)
		
      # compute the OR gates and calculate 'newInput'
			for(s in 1:nSp) { # for each species
				if(endIx[s] != 0) {	# if it is an output
					newInput[,s] <- apply(outputCube, 1, compOR) # find the maximum across reactions where the species is an output
					for(p in 1:length(newInput[,s])) { # for all conditions for each species
						if(!is.na(newInput[p,s])) { # if there is data
							if(newInput[p,s] == 0.5) { # if there is a delay
								newInput[p,s] = outputPrev[p,s] # set the species value to that of t-1
							}	
						}
					}
				}
			}
			      
		} else { # there is only 1 reaction TODO fix delay code here as well
			outputCube <- ifelse(all(is.na(tempStore)), NA, min(tempStore,na.rm=TRUE))
			newInput[,simList$maxIx] <- outputCube
		}
	
    # NOTE the stimuli here have values, different from c version
    # where 'new_input' is reinitialized so at this point (in c version) any
    # species that is not an output (i.e. stimuli) is NA
    
		# reset the inhibitors and stimuli
		for(stim in 1:length(indexList$stimulated)) {
			stimM <- cbind(CNOlist$valueStimuli[,stim], newInput[,indexList$stimulated[stim]])
			stimV <- apply(stimM, 1, maxNA)
			newInput[,indexList$stimulated[stim]] <- stimV
		}
		
		valueInhibitors <- 1-CNOlist$valueInhibitors
		newInput[,indexList$inhibited] <- valueInhibitors * newInput[,indexList$inhibited]
    
		# replace NAs with zeros to avoid having the NA penalty applying to unconnected species
		readout <- newInput
		readout[is.na(readout)] <- 0
		yBool[,,countBool] = readout
	}

	return(list(yBool, allCubes))

}
