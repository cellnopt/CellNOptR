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

gaBinaryTimeScaleT2 <- function(CNOlist, model, simList, indexList, bStringT1, simResT1, sizeFac=0.0001, NApenFac=1, popSize=50, pMutation=0.5, maxTime=60, maxGens=500, stallGenMax=100, selPress=1.2, elitism=5, relTol=0.1, verbose=TRUE, boolUpdates, divTime, lowerB=lowerB, upperB=upperB) {

	# find the bits to optimise
	bits2optimise <- which(bStringT1==0)

	# initialise
	bLength <- length(bits2optimise)
	pop <- round(matrix(runif(bLength*(popSize)), nrow=(popSize),ncol=bLength))
	bestBit <- pop[1,]
	bestObj <- Inf
	stop <- FALSE
	obj <- rep(0,popSize)
	g <- 0
	stallGen <- 0
	res<-rbind(c(g,bestObj, toString(bestBit), stallGen, Inf, Inf, toString(bestBit), 0), c(g,bestObj, toString(bestBit), stallGen, Inf, Inf, toString(bestBit), 0))
	colnames(res) <- c("Generation", "Best_score", "Best_bitString", "Stall_Generation", "Avg_score_Gen", "Best_score_Gen", "Best_bit_Gen", "Iter_time")
	popTol <- rep(NA, bLength)
	popTolScores <- NA

	# function that produces the score for a specific bitstring
	getObj <-  function(x) {
	
		bitString <- bStringT1
		bitString[which(bStringT1 == 0)] <- x
		bStringtimes <- bStringT1
		bStringtimes[which(bStringT1 == 0)] <- x*2
		
		# cut the model according to bitstring	
		modelCut <- model
		modelCut$interMat <- modelCut$interMat[,as.logical(bitString)]
		modelCut$notMat <- modelCut$notMat[,as.logical(bitString)]
		modelCut$reacID <- modelCut$reacID[as.logical(bitString)]
		modelCut$times <- bStringtimes[which(bStringtimes != 0)]
		simListCut <- simList
		simListCut$finalCube <- simListCut$finalCube[as.logical(bitString),]
		simListCut$ixNeg <- simListCut$ixNeg[as.logical(bitString),]
		simListCut$ignoreCube <- simListCut$ignoreCube[as.logical(bitString),]
		simListCut$maxIx <- simListCut$maxIx[as.logical(bitString)]
	
		# compute the score	
		getFitData <- getFitTimeScale(simList=simListCut, CNOlist=CNOlist, model=modelCut, indexList=indexList, boolUpdates, sizeFac=sizeFac, NApenFac=NApenFac, timeSplit="late", divTime, simresultsT1=simResT1, lowerB=lowerB, upperB=upperB)
		score = getFitData$score
		# ***** FIX nDataP *****
		nDataP <- sum(!is.na(CNOlist$valueSignals[[2]]))
		score <- score/nDataP
		return(score)
	}

	# loop
	t0 <- Sys.time()
	t <- t0
	
	while(!stop){
		
		# compute the scores
		scores <- apply(pop, 1, getObj)
		
		# fitness assignment: ranking, linear
		rankP <- order(scores, decreasing=TRUE)
		pop <- pop[rankP,]
		scores <- scores[rankP]
		fitness <- 2-selPress+(2*(selPress-1)*(c(1:popSize)-1)/(popSize-1))
		
		# selection:stochastic uniform sampling 
		wheel1 <- cumsum(fitness/sum(fitness))
		breaks <- runif(1)*1/popSize
		breaks <- c(breaks, breaks+((1:(popSize-1)))/popSize)
		sel <- rep(1, popSize)
		for(i in 1:length(breaks)) {
			sel[i] <- which(wheel1>breaks[i])[1]
		}
		
		# intermediate generation
		pop2 <- pop[sel,]
		pSize2 <- dim(pop2)[1]		
		pSize3 <- popSize-elitism
		
		# recombination: uniform: each bit has a 0.5 prob of being inherited from each parent
		mates <- cbind(ceiling(runif(pSize3)*pSize2), ceiling(runif(pSize3)*pSize2))
		
		# this holds the probability, for each bit, to be inherited from parent 1 (if TRUE) or 2 (if FALSE)
		inhBit <- matrix(runif((pSize3*bLength)), nrow=pSize3, ncol=bLength)
		inhBit <- inhBit < 0.5
		pop3par1 <- pop2[mates[,1],]
		pop3par2 <- pop2[mates[,2],]
		pop3 <- pop3par2
		pop3[inhBit] <- pop3par1[inhBit]
		
		# mutation
		mutProba <- matrix(runif((pSize3*bLength)), nrow=pSize3, ncol=bLength)
		mutProba <- (mutProba < (pMutation/bLength))
		pop3[mutProba] <- 1-pop3[mutProba]
		
		# compute stats
		t <- c(t, Sys.time())
		g <- g+1
		thisGenBest <- scores[length(scores)]
		thisGenBestBit <- pop[length(scores),]
		if(is.na(thisGenBest)){
			thisGenBest <- min(scores, na.rm=TRUE)
			thisGenBestBit <- pop[which(scores == thisGenBest)[1],]
		}
		if(thisGenBest < bestObj){
			bestObj <- thisGenBest
			bestBit <- thisGenBestBit
			stallGen <- 0
		} else {
			stallGen <- stallGen+1
		}
		
		resThisGen <- c(g, bestObj, toString(bestBit), stallGen, (mean(scores,na.rm=TRUE)), thisGenBest, toString(thisGenBestBit), as.numeric((t[length(t)]-t[length(t)-1]), units="secs"))		
		names(resThisGen) <- c("Generation", "Best_score", "Best_bitString", "Stall_Generation", "Avg_score_Gen", "Best_score_Gen", "Best_bit_Gen", "Iter_time")
		if(verbose) print(resThisGen)
		res <- rbind(res,resThisGen)
		
		# check stopping criteria
		criteria <- c((stallGen > stallGenMax), (as.numeric((t[length(t)]-t[1]), units="secs") > maxTime), (g > maxGens))
		if(any(criteria)) stop <- TRUE
		
		# check for bitstrings that are within the tolerance of the best bitstring
		tolScore <- scores[length(scores)]*0.1
		tolBS <- which(scores < scores[length(scores)]+tolScore)
		if(length(tolBS) > 0) {
			popTol <- rbind(popTol,pop[tolBS,])
			popTolScores <- c(popTolScores,scores[tolBS])
		}
		if(elitism > 0) {
			pop <- rbind(pop3, pop[(popSize-elitism+1):popSize,])
		} else {
			pop <- pop3
		}
	}
	
	popTol <- popTol[-1,]
	popTolScores <- popTolScores[-1]
	tolBS <- which(popTolScores < scores[length(scores)]+tolScore)
	popTol <- popTol[tolBS,]
	popTolScores <- popTolScores[tolBS]
	popTolT <- cbind(popTol, popTolScores)
	popTolT <- unique(popTolT, MARGIN=1)
	if(!is.null(dim(popTolT))) { 
		popTol <- popTolT[,1:(dim(popTolT)[2]-1)]
		popTolScores <- popTolT[,dim(popTolT)[2]]
	} else {
		popTol <- popTolT[1:(length(popTolT)-1)]
		popTolScores <- popTolT[length(popTolT)]
	}	
	
	res <- res[3:dim(res)[1],]	
	rownames(res) <- NULL
	return(list(bString=bestBit, results=res, stringsTol=popTol, stringsTolScores=popTolScores))

}
