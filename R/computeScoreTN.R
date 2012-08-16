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

# Function that computes the score of a specific bitstring
# Although it is very similar to computeScoreT1, there are enough differences to
# have a different function.
computeScoreTN<-function(CNOlist, model, simList=NULL, indexList=NULL, simResPrev, bStringPrev, 
    bStringNext, timeIndex=3, sizeFac=0.0001, NAFac=1){

    # by default same behaviour as computeScoreT2
    # timeIndex=3 stands for T2 by default.
    #timeIndex = timeIndex # i.e., "tN"
    if (is.null(simList)==TRUE){
        simList = prep4sim(model)
    }
    if (is.null(indexList)==TRUE){
        indexList = indexFinder(CNOlist, model)
    }

    bitString <- bStringPrev
    bitString[which(bStringPrev == 0)] <- bStringNext
    bStringTimes = bStringPrev
    bStringTimes[which(bStringPrev == 0)] <- bStringNext * (timeIndex-1)


    modelCut = cutModel(model, bitString)
    modelCut$times <- bStringTimes[which(bStringTimes != 0)]
    simListCut <- cutSimList(simList, bitString)

    # Compute the simulated results
    simResults <- simulatorTN(
        simResultsPrev=simResPrev,
        CNOlist=CNOlist,
        model=modelCut,
        simList=simListCut,
        indexList=indexList,
        timeIndex=timeIndex)

    #Compute the score
    Score <- getFit(
        simResults=simResults,
        CNOlist=CNOlist,
        model=modelCut,
        indexList=indexList,
        timePoint=timeIndex,
        sizeFac=sizeFac,
        NAFac=NAFac,
        nInTot=length(which(model$interMat == -1)),
		simResultsT0=NA)

  nDataP <- sum(!is.na(CNOlist$valueSignals[[2]]))
  Score <- Score/nDataP


  return(Score)
}
