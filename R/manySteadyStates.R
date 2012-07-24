#takes as an input the CUT (NCNOcutCompExp and indicesNCNOcutCompExp) model

manySteadyStates <-function(
  CNOlist,
  model,
  simList,
  indexList,
  sizeFac=0.0001,
  NAFac=1,
  popSize=50,
  Pmutation=0.5,
  maxTime=60,
  maxGens=500,
  stallGenMax=100,
  selPress=1.2,
  elitism=5, 
  relTol=0.1,
  verbose=FALSE,
  priorBitString=NULL){

initBstring<-rep(1,length(model$reacID))

Opt<-list()
bString<-list()
simRes<-list()
T1opt<-gaBinaryT1(CNOlist=CNOlist,
                  model=model,
                  simList=simList,
                  indexList=indexList,
                  initBstring=initBstring,
                  stallGenMax=stallGenMax,
                  sizeFac=sizeFac,
                  NAFac=NAFac,
                  popSize=popSize,
                  pMutation=Pmutation,
                  maxTime=maxTime,
                  maxGens=maxGens,
                  selPress=selPress,
                  elitism=elitism,
                  relTol=relTol,
                  verbose=verbose,
                  priorBitString=priorBitString)

Opt[[2]]<-T1opt
simT1<-simulateT1(CNOlist=CNOlist,
                  model=model,
                  bStringT1=T1opt$bString,
                  simList=simList,
                  indexList=indexList)
simRes[[2]]<-simT1
currBstring<-T1opt$bString    
bStringTimes<-T1opt$bString

for(i in 3:length(CNOlist$valueSignals)){
  timeIndex<-i
  Opt[[i]]<-gaBinaryTN(CNOlist=CNOlist,
                    model=model,
                    simList=simList,
                    indexList=indexList,
                    bStringPrev=currBstring,
                    simResPrev=simRes[[i-1]],
                    timeIndex=i,
                    bStringTimes=bStringTimes,
                    stallGenMax=stallGenMax,
                    maxTime=maxTime,
                    sizeFac=sizeFac,
                    NAFac=NAFac,
                    popSize=popSize,
                    pMutation=Pmutation,
                    maxGens=maxGens,
                    selPress=selPress,
                    elitism=elitism,
                    relTol=relTol,
                    verbose=verbose,
                    priorBitString=priorBitString)
  
  currBstring[which(currBstring==0)]<-Opt[[i]]$bString
  bStringTimes[which(bStringTimes==0)]<-Opt[[i]]$bString*(timeIndex-1)
  
  simRes[[i]]<-simulateTN(CNOlist,model,simRes[[i-1]],currBstring,bStringTimes,simList,indexList,timeIndex)
  
}
                  return(list(bString=currBstring,
                              bStringTimes=bStringTimes,
                              Opt=Opt,
                              simRes=simRes))               
}

