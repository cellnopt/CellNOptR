load_all("../CNORdt")
library(MEIGOR)


# DATA AND MODEL ----------------------------------------------------------

data(CNOlistPB, package="CNORdt")
data(modelPB, package="CNORdt")
cnolist = CNOlistPB
modelInit = preprocessing(cnolist, modelPB)
# the dt optimized bit string
bitString = c(1,0,0,0,1,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0)
model = cutModel(modelInit, bitString)
simList = prep4sim(model)
indexList = indexFinder(cnolist, model)

# load c file
dyn.load("simulatorDelay.so")


# TEST SIMULATOR ----------------------------------------------------------

boolUpdates = 30
delayThresh = rep(0,length(model$reacID))
negEdges = rep(0,length(model$reacID))
negEdges[5] = 1
delayThresh[5] = 5
delayThresh[9] = 3

# simulate with delays, r version
# TODO, some mistakes in this... (e.g. [10,6], [1,1:3])
simD1 = simulatorDelayR(cnolist, model, simList, indexList,
boolUpdates=boolUpdates, delayThresh=delayThresh, strongWeak=negEdges)
simD1 = simD1[[1]][,indexList$signals,]
plotCNOlist(plotData(cnolist, simD1))

# c version
simT0 = simulatorT0(cnolist, model, simList, indexList)
simD2 = simulatorDelay(cnolist, model, simList, indexList,
boolUpdates=boolUpdates, delayThresh=delayThresh, strongWeak=negEdges)
simD2 = convert2array(simD2, 10, 12, boolUpdates)
simD2 = simD2[,indexList$signals,]
plotCNOlist(plotData(cnolist, simD2))

# compare DT
simDT = simulatorDT(cnolist, model, simList, indexList, boolUpdates=boolUpdates, prevSim=NULL)
simDT = convert2array(simDT, 10, 12, boolUpdates)
simDT = simDT[,indexList$signals,]
plotCNOlist(plotData(cnolist, simDT))

res1 = system.time(replicate(1000, simulatorDelay(cnolist, model, simList, indexList,
boolUpdates=boolUpdates, delayThresh=delayThresh, strongWeak=negEdges)))
res2 = system.time(replicate(1000, simulatorDelayR(cnolist, model, simList, indexList,
boolUpdates=boolUpdates, delayThresh=delayThresh, strongWeak=negEdges)))
res2[1] / res1[1]


# TEST OPTIMIZATION -------------------------------------------------------

fit1 = getFitDelay(cnolist, model, simList, indexList, boolUpdates=boolUpdates, nInTot=length(which(modelInit$interMat==-1)))

# check computeScoreDelay
score1 = computeScoreDelay(CNOlist=cnolist, model=modelInit, bString=bitString, boolUpdates=30)
score2 = computeScoreDT(cnolist, modelInit, bString=bitString, boolUpdates=30, lowerB=0.8, upperB=10)

# check cutAndPlotResultsDelay
cutAndPlotResultsDelay(
  model=modelInit,
  CNOlist=cnolist,
  bString=bitString,
  plotPDF=FALSE,
  boolUpdates=30
)






# check the improvement in score over getFitDT
boolUpdates = 30
simResults <- simulatorDT(
  CNOlist=CNOlistPB,
  model=model,
  simList=simList,
  indices=indexList,
  boolUpdates=boolUpdates
)
simResults = convert2array(simResults, dim(CNOlistPB$valueSignals[[1]])[1],
length(model$namesSpecies), boolUpdates)

fit2 <- getFitDT(
  simResults=simResults,
  CNOlist=CNOlistPB,
  model=model,
  indexList=indexList,
  boolUpdates=boolUpdates,
  lowerB=0.8,
  upperB=10,
  nInTot=length(which(modelPB$interMat == -1))
)

fit2$score/fit1$score
