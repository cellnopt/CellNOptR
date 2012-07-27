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
computeScoreT2<-function(CNOlist, model, simList, indexList, simResT1, bStringT1,
    bStringT2, sizeFac=0.0001, NAFac=1){

	timeIndex = 3 # i.e., "t2"

    bitString <- bStringT1
    bitString[which(bStringT1 == 0)] <- bStringT2
    bStringTimes <- bStringT1
    bStringTimes[which(bStringT1 == 0)] <- bStringT2 * (timeIndex-1)

    modelCut = cutModel(model, bitString)
    modelCut$times <- bStringTimes[which(bStringTimes != 0)]
    simListCut <- cutSimList(simList, bitString)

    # Compute the simulated results
    simResults <- simulatorT2(
        simResT1,
        CNOlist=CNOlist,
        model=modelCut,
        simList=simListCut,
        indexList=indexList)

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
