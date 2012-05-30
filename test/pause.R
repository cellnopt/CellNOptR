library(CellNOptR)
library(CNORdt)

source("../CNORdt/R/simulatorPause.R")
source("../CNORdt/R/getFitPause.R")

load("CNOlistDelay.RData")
modelDelay = readSif("toyDelay.sif")

indexDelay = indexFinder(CNOlist=CNOlistDelay, Model=modelDelay, verbose=T)
simDelay = prep4Sim(modelDelay)
delayThresh = rep(0,6)

# add feedback finder here

source("../CNORdt/R/feedbackFinder1.R")
source("../CNORdt/R/feedbackWrapper.R")

negEdges = feedbackWrapper(modelDelay)

# simulate with delays
source("../CNORdt/R/simulatorPause.R")
sim1 = simulatorPause(CNOlistDelay, modelDelay, simDelay, indexDelay,
boolUpdates=30, delayThresh, strongWeak=negEdges)
sim1 = sim1[,indexDelay$signals,]


plotCNOlist(plotData(CNOlistDelay, sim1))


# make list from bool.sim
valueSignals = list()
for(a in 1:dim(sim1)[3]) {
	valueSignals = c(valueSignals, list(sim1[,,a]))
}
CNOlistBool = CNOlistDelay
CNOlistBool$valueSignals = valueSignals
CNOlistBool$timeSignals = 0:(boolUpdates-1)
plotCNOlist(CNOlistBool)



sim1 = system.time(replicate(100,simulatorPause(CNOlistDelay, modelDelay, simDelay, indexDelay,
boolUpdates=10, delayThresh)))
sim1 = system.time(simulatorTimeScale(CNOlistDelay, modelDelay, simDelay, indexDelay,boolUpdates=10))
system.time(replicate(100,simulatorT1(CNOlistDelay, modelDelay, simDelay, indexDelay)))

plotData <- function(CNOlist, data) {

    valueSignals = list()
    outSim = data
    if(is.list(outSim)) {
        valueSignals = outSim
        t = length(outSim)
    } else {   
        for(a in 1:dim(outSim)[3]) {
            valueSignals = c(valueSignals, list(outSim[,,a]))
        }
        t = dim(outSim)[3]
    }   
    CNOlistSim = CNOlist
    CNOlistSim$valueSignals = valueSignals
    CNOlistSim$timeSignals = 0:(t-1)

    return(CNOlistSim)   

}





optDelay = list()
optDelay$bString = c(1,0,0,0,1,0,1,1,0,1,1,1,1,1,1,0,0,0,1,1)
cutData = cutModel(ModelCutCompressExpand,SimList=fields4Sim,bitString=optDelay$bString)
simDelay = cutData[[2]]
ModelDelay = cutData[[1]]

indexDelay = indexFinder(CNOlist=CNOlist, Model=ModelDelay, verbose=T)


load("objects/CNOlist.RData")
load("objects/ModelDelay.RData")
load("objects/simDelay.RData")
load("objects/indexDelay.RData")
source("objects/feedbackFinder1.R")

delayThresh = rep(0,length(ModelDelay$reacID))
boolUpdates = 30
delayThresh[c(11,12)] = c(10,5)

source("../CNORdt/R/simulatorPause.R")
source("../CNORdt/R/getFitPause.R")
gfPause = getFitPause(SimList=simDelay, CNOlist=CNOlist, Model=ModelDelay, indexList=indexDelay, sizeFac=0.0001, NAPenFac=1, boolUpdates=boolUpdates)

outSim = simulatorPause(CNOlist, ModelDelay, simDelay, indexDelay, boolUpdates=30, delayThresh)
#outSim = simulatorTimeScale(CNOlist, ModelDelay, simDelay, indexDelay, boolUpdates=20)

outSim = outSim[,indexDelay$signals,]

# make list from bool.sim
valueSignals = list()
for(a in 1:dim(outSim)[3]) {
	valueSignals = c(valueSignals, list(outSim[,,a]))
}
CNOlistBool = CNOlist
CNOlistBool$valueSignals = valueSignals
CNOlistBool$timeSignals = 0:(boolUpdates-1)
plotCNOlist(CNOlistBool)

#####

SimList= simDelay
Model = ModelDelay
indexList = indexDelay
sizeFac=0.0001
NAPenFac=1
timeSplit="early"
divTime=NULL
