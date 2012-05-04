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

#Function that computes the score of a specific bitstring
# todo: this is similar to wha is done in gaBinaryT1. need to do the same for T2
# todo: timepoints and nafac are hardcoded.
computeScore<-function(CNOlist,Model,SimList,indexList,x){

  bitString<-x
  
  #cut the model according to bitstring
  ModelCut<-Model
  ModelCut$interMat<-ModelCut$interMat[,as.logical(bitString)]
  ModelCut$notMat<-ModelCut$notMat[,as.logical(bitString)]
  ModelCut$reacID<-ModelCut$reacID[as.logical(bitString)]
  
  SimListCut<-cutSimList(SimList, bitString)
  
  #compute the simulated results
  SimResults<-simulatorT1(
    CNOlist=CNOlist,
    Model=ModelCut,
    SimList=SimListCut,
    indexList=indexList)
  
  # We may want to to use the T0 information.
  SimResultsT0<-simulatorT0(
    CNOlist=CNOlist,
    Model=ModelCut,
    SimList=SimListCut,
    indexList=indexList)
  
  #Compute the score
  Score<-getFit(
    SimResults=SimResults,
    SimResultsT0=SimResultsT0,
    CNOlist=CNOlist,
    Model=ModelCut,
    indexList=indexList,
    timePoint="t1",
    sizeFac=0.0001,
    NAFac=1,
    nInTot=length(which(Model$interMat == -1)))
  nDataP<-sum(!is.na(CNOlist$valueSignals[[2]]))
  Score<-Score/nDataP
  
  return(Score)
}
