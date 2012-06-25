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
computeScoreT2<-function(CNOlist, Model, SimList, indexList, SimResT1, bStringT1,
    bStringT2, sizeFac=0.0001, NAFac=1){

    bitString <- bStringT1
    bitString[which(bStringT1 == 0)]<- bStringT2
    bStringTimes <- bStringT1
    bStringTimes[which(bStringT1 == 0)]<-bStringT2*2

    ModelCut = cutModel(Model, bitString)
    ModelCut$times <- bStringTimes[which(bStringTimes != 0)]
    SimListCut <- cutSimList(SimList, bitString)


    # We may want to to use the T0 information.
    SimResults <- simulatorT2(
        SimResultsT1=SimResT1,
        CNOlist=CNOlist, 
        Model=ModelCut,
        SimList=SimListCut,
        indexList=indexList)

    #Compute the score
    Score <- getFit(
        SimResults=SimResults,
        CNOlist=CNOlist,
        Model=ModelCut,
        indexList=indexList,
        timePoint="t2",
        sizeFac=sizeFac,
        NAFac=NAFac,
        nInTot=length(which(Model$interMat == -1)))

  nDataP <- sum(!is.na(CNOlist$valueSignals[[2]]))
  Score <- Score/nDataP

  return(Score)
}
