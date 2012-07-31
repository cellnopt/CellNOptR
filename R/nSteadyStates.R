#A script to run multiple timepoints (after T1 has already been run)
#Assume want to have a steady state for every input time point (if not the case, can just subset CNOlist)
nSteadyStates <-function(
  simResultsT1,
  bStringT1,
  CNOlist,
  Model,
  SimList,
  indexList,
  sizeFac=0.0001,
  NAFac=1,
  PopSize=50,
  Pmutation=0.5,
  MaxTime=60,
  maxGens=500,
  StallGenMax=100,
  SelPress=1.2,
  elitism=5, 
  RelTol=0.1){
  
  
  #initialize Opt array
  Opt<-list()
  #initialze a bString array
  bString<-list()
  bString[[2]]<-bStringT1
  #initialize a SimResults array
  simRes<-list()
  simRes[[2]]<-simResultsT1
  
  # if(length(CNOlist$valueSignals)<3)
  # {
  # stop("There is not more than one time point to optimize")
  # }
  bStringPrev<-bStringT1
  for(i in 3:length(CNOlist$valueSignals))
  {
    Opt[[i]]<-gaBinaryTN(CNOlist=CNOlist,
                         Model=Model,
                         indexList=indexList,
                         bStringT1=bString[[i-1]],
                         SimResT1=simRes[[i-1]],
                         SimList=SimList,
                         timePoint=i,
                         verbose = FALSE,
                         sizeFac=0.0001,
                         NAFac=1,
                         PopSize=50,
                         Pmutation=0.5,
                         MaxTime=60,
                         maxGens=500,
                         StallGenMax=100,
                         SelPress=1.2,
                         elitism=5, 
                         RelTol=0.1)
    bString[[i]]<-Opt[[i]]$bString
#     bitString2 <- bStringPrev
#     bitString2[which(bStringPrev == 0)] <- bString[[i]]
    simRes[[i]]<-simulateTN(CNOlist=CNOlist,
                                Model=Model,
                                bStringT2=bString[[i]],
                                bStringT1=bString[[i-1]],
                                SimList=SimList,
                                indexList=indexList,
                                SimResT1=simRes[[i-1]])
    #                           TimePoint=(i-1))
    bStringPrev<-bitString2
  }
  
  return(list(bString=bString,
              Opt=Opt,
              simRes=simRes))
}
