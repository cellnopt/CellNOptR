#A script to run multiple timepoints (after T1 has already been run)
#Assume want to have a steady state for every input time point (if not the case, can just subset CNOlist)
nSteadyStates <-function(
  simResultsT1,
  bStringT1,
  CNOlist,
  model,
  simList,
  indexList,
  sizeFac=0.0001,
  NAFac=1,
  popSize=50,
  pMutation=0.5,
  maxTime=60,
  maxGens=500,
  stallGenMax=100,
  selPress=1.2,
  elitism=5, 
  relTol=0.1){
  
  
    #initialize Opt array
    Opt<-list()

    #initialze a bString array
    bStrings<-list()
    bStrings[[1]]<-bStringT1

    #initialize a SimResults array
    simRes<-list()
    simRes[[1]]<-simResultsT1
  

    bStringPrev<-bStringT1
    for(i in 3:length(CNOlist$valueSignals))
    {
      Opt[[i]]<-gaBinaryTN(CNOlist=CNOlist,
                         model=Model,
                         bStrings=bitStrings,
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

        bitStrings[[i-1]] = Opt[[i-1]]$bString
        simRes[[i]]<-simulateTN(CNOlist,model,bStrings)
  }
  
  return(list(bString=bString, Opt=Opt, simRes=simRes))
}
