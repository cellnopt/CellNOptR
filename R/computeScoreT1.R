#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EMBL - European Bioinformatics Institute
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id$

#Function that computes the score of a specific bitstring
# todo: this is similar to wha is done in gaBinaryT1. need to do the same for T2
computeScoreT1<-function(CNOlist, model, bString, simList=NULL, indexList=NULL, 
    sizeFac=0.0001, NAFac=1, timeIndex=2){
    # simList and indexList are computed inside this function. 
    # However, for back-compatibility, we keep the arguments so that if
    # provided, we can still use them.

    # uncomment if required for debugging.
    # if (timeIndex<2){ stop("timeIndex must be >=2")}
    # if (timeIndex>length(CNOlist@timeSignals)){ 
    #      stop(paste("timeIndex must be <= ", length(CNOlist@timeSignals),sep=" "))
    # }


    if (is.null(simList)==TRUE){
        simList = prep4sim(model)
    }
    if (is.null(indexList)==TRUE){
        indexList = indexFinder(CNOlist, model)
    }

	# TC oct 2012. cutModel is a function. It performs the cut to select interMat, 
    # reacID, nameSpecies, and notMat. Here, we need only need the 3  and possibly only 2, so let us 
    # just copy and paste the code. This will save 10% of computational time
    # modelCut = cutModel(model, bString)
    bs = as.logical(bString)
    modelCut <- list()
    modelCut$interMat <- model$interMat[, bs]
    modelCut$reacID <- model$reacID[bs]


    simListCut <- cutSimList(simList, bString)


    # Compute the simulated results
    nStimuli = length(indexList$stimuli)
    nInhibitors <- length(indexList$inhibited)
    nCond <- dim(CNOlist@stimuli)[1]
    nReacs <- length(modelCut$reacID)
    nSpecies <- length(model$namesSpecies) # this is correct. No need to get modelCut$namesSpecies 
    nMaxInputs <- dim(simListCut$finalCube)[2]
    nStimuli <- dim(CNOlist@stimuli)[2]

    # simList matrices. C code must handle the matrix indices carefully.
    # This is faster than transforming into a vector as in the previous code.
    finalCube = as.integer(simListCut$finalCube-1)
    ixNeg = as.integer(simListCut$ixNeg)
    ignoreCube = as.integer(simListCut$ignoreCube)
    maxIx = as.integer(simListCut$maxIx-1)

    # index. convertion from R to C indices convention.
    indexSignals <- as.integer(indexList$signals-1)
    indexStimuli <- as.integer(indexList$stimulated-1)
    indexInhibitors <- as.integer(indexList$inhibited-1)
    nSignals <- length(indexSignals)


    # cnolist
    valueInhibitors <- as.integer(CNOlist@inhibitors)
    valueStimuli <- as.integer((CNOlist@stimuli))

	simResults = .Call("simulatorT1", nStimuli, nInhibitors,
		nCond, nReacs, nSpecies, nSignals, nMaxInputs,
        finalCube, ixNeg, ignoreCube, maxIx,
        indexSignals, indexStimuli, indexInhibitors, valueInhibitors,
        valueStimuli, as.integer(1))

    simResultsT0 = .Call("simulatorT1", nStimuli, nInhibitors,
        nCond, nReacs, nSpecies, nSignals, nMaxInputs,
        finalCube, ixNeg, ignoreCube, maxIx,
        indexSignals, indexStimuli, indexInhibitors, 
        valueInhibitors, valueStimuli, as.integer(0))


    # this step is commented in the C code. Does not seem to work properly. Try
    # test_simulateTN.R for instance.
    #simResultsT0 = simResultsT0[, indexList$signals]
    #simResults = simResults[, indexList$signals]
    nInTot = length(which(model$interMat == -1))

    #Compute the score
    mode = 1 # 1 for TRUE: takes into account T0
    Score = .Call("getFit", 
        nCond,
        nSignals, 
        nReacs,
        nSpecies,
        sum(bs),
        nInTot,
        as.real(simResultsT0), 
        as.real(simResults), 
        as.real(CNOlist@signals[[1]]), #T0 data
        as.real(CNOlist@signals[[timeIndex]]), #T1 data
        as.integer(modelCut$interMat),
		as.real(sizeFac),
		as.real(NAFac), as.integer(mode)
		)[[1]]

  return(Score)
}
