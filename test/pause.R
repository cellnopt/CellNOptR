library(CellNOptR)
library(CNORdiscreteTime)
library(CNORode)
library(CNOR.CFL)

# load model and cnolist



optDelay = list()
optDelay$bString = c(0,1,1,1,1,0,0,0,1,0,0,1,1,1,1,1,1,0,0,0,1,0)
cutData = cutModel(ModelCutCompressExpand,SimList=fields4Sim,bitString=optDelay$bString)
simDelay = cutData[[2]]
ModelDelay = cutData[[1]]
indexDelay = indexFinder(CNOlist=CNOlist, Model=ModelDelay, verbose=T)
delayThresh = rep(0,length(ModelDelay$reacID))
boolUpdates = 30
delayThresh = c(0, 8 ,0 ,2 ,0 ,20, 0 ,20, 6 ,0, 0, 3)

source("../../../trunk/CNOR_dt/R/simulatorPause.R")
source("../../../trunk/CNOR_dt/R/getFitPause.R")

gfPause = getFitPause(SimList=simDelay, CNOlist=CNOlist, Model=ModelDelay, indexList=indexDelay, sizeFac=0.0001, NAPenFac=1, boolUpdates=boolUpdates)

outSim = simulatorPause(CNOlist, ModelDelay, simDelay, indexDelay, boolUpdates=20, delayThresh)
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












source("functionEdits/cutAndPlotResultsT1.R")
source("functionEdits/cutAndPlotResultsT2.R")
source("functionEdits/cutAndPlotResultsCFL.R")
source("functionEdits/plotOptimResultsNew.R")
source("functionEdits/plotOptimResultsNewT2.R")
source("functionEdits/plotOptimResultsNewCFL.R")
source("functionEdits/plotCNOlistOver.R")
source("functionEdits/cutAndPlotResultsTimeScaleT2.R")
source("functionEdits/plotOptimResultsTimeScale.R")
source("functionEdits/getFitTimeScale.R")

# add model and CNOlist:
Model = readSif("model/ModelV2PKN.sif")

# CNOlist: data after noise addition and time point selection
# CNOlistInSilico: data before noise addition and time point selection 
load("data/CNOlist.Rdata")
load("data/CNOlistInSilico.Rdata")

# visualize data (figure 3)
pdf(file="./results/figure3.pdf",width=14.5,height=11)
plotCNOlistOver(CNOlist=CNOlistInSilico,data2=CNOlist)
dev.off()


########## PROCESSING ##########


# processing the model (indexing, compression and expansion)
# index the stimuli, readouts and inhibitors
indexOrig <- indexFinder(CNOlist=CNOlist, Model=Model, verbose=T)

# find the indexes of the non-observables and the non-controllable species
indexNONC <- findNONC(Model=Model, indexes=indexOrig, verbose=T)

# cut the data according to 'indexNONC'
ModelCut <- cutNONC(Model=Model, NONCindexes=indexNONC)

# find the indexes again as model may have changed
indexNONCcut <- indexFinder(CNOlist=CNOlist, Model=ModelCut)

# compress the model
ModelCutCompress <- compressModel(Model=ModelCut, indexes=indexNONCcut)

# find the indexes again
indexNONCcutComp <- indexFinder(CNOlist=CNOlist, Model=ModelCutCompress)

# expand
ModelCutCompressExpand <- expandGates(Model=ModelCutCompress)

# extract information for simulation
fields4Sim <- prep4Sim(Model=ModelCutCompressExpand)
initBstring <- rep(1, length(ModelCutCompressExpand$reacID))


########## /PROCESSING/ ##########


# STEADY STATE
# what time point is 'steady state' in the data?

t = 10
CNOlistSS = CNOlist
tIndex = which(CNOlistSS$timeSignals == t)
# make a new CNOlist with a single measurement time point
CNOlistSS$timeSignals = c(0,t)
CNOlistSS$valueSignals = list(CNOlist$valueSignals[[1]], CNOlist$valueSignals[[tIndex]])

optPart1 <- gaBinaryT1(CNOlist=CNOlistSS, Model=ModelCutCompressExpand,
SimList=fields4Sim, indexList=indexNONCcutComp,
initBstring=initBstring, verbose = TRUE, MaxTime=180)

pdf(file="./results/figure4.pdf",width=14.5,height=11)
cutAndPlotResultsT1(Model=ModelCutCompressExpand, bString=optPart1$bString, SimList=fields4Sim,
CNOlist=CNOlistSS, indexList=indexNONCcutComp, plotPDF=FALSE, CNOlistDT=CNOlist)
dev.off()


#  2 STEADY STATE

t = c(10,30)
CNOlistSS2 = CNOlist
tIndex = which(CNOlistSS2$timeSignals == t[1])
tIndex[2] = which(CNOlistSS2$timeSignals == t[2])
# make a new CNOlist with 2 time points
CNOlistSS2$timeSignals = c(0,t)
CNOlistSS2$valueSignals = list(CNOlist$valueSignals[[1]], CNOlist$valueSignals[[tIndex[1]]], CNOlist$valueSignals[[tIndex[2]]])

optPart1B1 <- gaBinaryT1(CNOlist=CNOlistSS2, Model=ModelCutCompressExpand,
SimList=fields4Sim, indexList=indexNONCcutComp,
initBstring=initBstring, verbose=TRUE, MaxTime=180)

# optimise T2
SimT1 <- simulateT1(CNOlist=CNOlistSS2, Model=ModelCutCompressExpand, bStringT1=optPart1B1$bString,
SimList=fields4Sim, indexList=indexNONCcutComp)

optPart1B2 <- gaBinaryT2(CNOlist=CNOlistSS2, Model=ModelCutCompressExpand, SimList=fields4Sim,
indexList=indexNONCcutComp, bStringT1=optPart1B1$bString, SimResT1=SimT1)

pdf(file="./results/figure5.pdf",width=14.5,height=11)
cutAndPlotResultsT2(Model=ModelCutCompressExpand, bStringT1=optPart1B1$bString, bStringT2=optPart1B2$bString, SimList=fields4Sim,
CNOlist=CNOlistSS2, indexList=indexNONCcutComp, CNOlistDT=CNOlist)
dev.off()


# 2 STEADY STATE - TIME COURSE

optPart3 <- gaBinaryTimeScale(CNOlist=CNOlist, Model=ModelCutCompressExpand,
SimList=fields4Sim, indexList=indexNONCcutComp,
initBstring=initBstring, verbose=TRUE, boolUpdates=c(10,20), divTime=10, MaxTime=180)

# simulate the above model to get the starting point for second steady state
dataStartPoint = cutModel(Model=ModelCutCompressExpand, SimList=fields4Sim, bitString=optPart3$bString)
simT1 = simulatorTimeScale(CNOlist=CNOlist, Model=dataStartPoint[[1]], SimList=dataStartPoint[[2]],
indexList=indexNONCcutComp, boolUpdates=10)

# optimize the 'second' network then stitch the results together
optPart4 <- gaBinaryTimeScaleT2(CNOlist=CNOlist, Model=ModelCutCompressExpand,
SimList=fields4Sim, indexList=indexNONCcutComp,
bStringT1=optPart3$bString, SimResT1=simT1, verbose=TRUE, boolUpdates=c(10,20), divTime=10)

# visualize total result
pdf(file="./results/figure6.pdf",width=14.5,height=11)
cutAndPlotResultsTimeScaleT2(Model=ModelCutCompressExpand, bStringT1=optPart3$bString, bStringT2=optPart4$bString,
SimList=fields4Sim, CNOlist=CNOlist, indexList=indexNONCcutComp, boolUpdates=c(10,20),
divTime=10)
dev.off()


# ODE - HERE WE UNCOMPRESS THE FEEDBACK

# take the best topology according to the other formalisms
optPart5 = list()
optPart5$bString = c(0,1,1,1,1,0,0,0,1,0,0,1,1,1,1,1,1,0,0,0,1,1) 
ModelODE = cutModel(ModelCutCompressExpand,SimList=fields4Sim,bitString=optPart5$bString)[[1]]
indicesODE <- indexFinder(CNOlist=CNOlist, Model=ModelODE) 

paramsInit = createLBodeContPars(ModelODE,random=TRUE)
paramsOpt = parEstimationLBodeSSm(CNOlist, ModelODE, paramsInit,
indices=indicesODE, maxtime=2000, ndiverse=100, dim_refset=10, maxStepSize=0.05)

source("functions/CNORdiscreteTimeSandBox/plotOptimResultsTimeScale.R")
source("functions/tutorialEdits/plotLBodeFitness.R")
pdf(file="out.pdf",width=14.5,height=11)
plotLBodeFitness(CNOlist, ode_parameters=paramsOpt, ModelODE, indices=indicesODE, maxStepSize=0.05)
dev.off()


# CONSTRAINED FUZZ LOGIC

# run the wrapper first
optCFL = CNORwrapFuzzy(Data=CNOlistSS,Model=Model)

a=2
CFLbestModel = optCFL$RedRef[[a]]$RefModel$RefinedModel
CFLindex = indexFinder(CNOlistSS, CFLbestModel)
CFLsimList = optCFL$RedRef[[a]]$RefModel$RefinedSimList
simCFL = simFuzzyT1(CNOlistSS, CFLbestModel, CFLsimList, CFLindex)
iSignals = match(CNOlistSS$namesSignals,colnames(simCFL))
simCFL = simCFL[,iSignals]

source("functions/tutorialEdits/plotOptimResultsNewCFL.R")
pdf(file="out.pdf",width=14.5,height=11)
plotOptimResultsNewCFL(SimResults=simCFL, xCoords=NULL, CNOlist=CNOlist)
dev.off()
