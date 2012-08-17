#takes as an input the CUT (NCNOcutCompExp and indicesNCNOcutCompExp) model

manySteadyStates <-function(
  CNOlist,
  model,
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
    simRes<-list()
    bitStrings  = list()

    T1opt<-gaBinaryT1(CNOlist=CNOlist,
                  model=model,
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

    Opt[[1]]<-T1opt
    simT1<-simulateTN(CNOlist=CNOlist, model=model, bStrings=list(T1opt$bString))
    simRes[[1]]<-simT1
    bitStrings[[1]] = T1opt$bString


    for(i in 3:length(CNOlist$valueSignals)){
        Opt[[i-1]]<-gaBinaryTN(CNOlist=CNOlist,
                    model=model,
                    bStringPrev=bitStrings,
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

        bitStrings[[i-1]] = Opt[[i-1]]$bString

        simRes[[i]]<-simulateTN(CNOlist,model,bStrings)

    }
    return(list(bStrings=bStrings, Opt=Opt, simRes=simRes))
}

