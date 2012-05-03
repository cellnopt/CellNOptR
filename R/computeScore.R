#Function that produces the score for a specific bitstring
getObj<-function(CNOlist,Model,SimList,IndexList,x){

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