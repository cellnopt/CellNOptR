library(CellNOptR)
library(CNORdt)

data(CNOlistPB, package="CNORdt")
data(ToyModelPB, package="CNORdt")

# processing the model (indexing, compression and expansion)
# index the stimuli, readouts and inhibitors
indexOrig <- indexFinder(CNOlist=CNOlist, model=Model, verbose=T)

# find the indexes of the non-observables and the non-controllable species
indexNONC <- findNONC(model=Model, indexes=indexOrig, verbose=T)

# cut the data according to 'indexNONC'
ModelCut <- cutNONC(model=Model, NONCindexes=indexNONC)

# find the indexes again as model may have changed
indexNONCcut <- indexFinder(CNOlist=CNOlist, model=ModelCut)

# compress the model
ModelCutCompress <- compressModel(model=ModelCut, indexes=indexNONCcut)

# find the indexes again
indexNONCcutComp <- indexFinder(CNOlist=CNOlist, model=ModelCutCompress)

# expand
ModelCutCompressExpand <- expandGates(model=ModelCutCompress)

# extract information for simulation
fields4Sim <- prep4sim(model=ModelCutCompressExpand)
initBstring <- rep(1, length(modelCutCompressExpand$reacID))

opt1 <- gaBinaryTimeScale(CNOlist=CNOlist, model=ModelCutCompressExpand,
simList=fields4Sim, indexList=indexNONCcutComp, initBstring=initBstring,
verbose=TRUE, boolUpdates=c(10,20), divTime=10, MaxTime=100, lowerB=0.8, upperB=10)

cutAndPlotResultsTimeScale(model=ModelCutCompressExpand, bString=opt1$bString, simList=fields4Sim,
CNOlist=CNOlist, indexList=indexNONCcutComp, boolUpdates=c(10,20), divTime=10, lowerB=0.8, upperB=10)

dataStartPoint = cutModel(model=ModelCutCompressExpand,
SimList=fields4Sim, bitString=opt1$bString)
simT1 = simulatorTimeScale(CNOlist=CNOlist,
model=dataStartPoint[[1]], simList=dataStartPoint[[2]],
indexList=indexNONCcutComp, boolUpdates=10)

opt2 <- gaBinaryTimeScaleT2(CNOlist=CNOlist,
model=ModelCutCompressExpand, simList=fields4Sim,
indexList=indexNONCcutComp, bStringT1=opt1$bString,
simResT1=simT1, verbose=TRUE, boolUpdates=c(10,20),
divTime=10, lowerB=0.8, upperB=10)

cutAndPlotResultsTimeScaleT2(model=ModelCutCompressExpand,
bStringT1=opt1$bString, bStringT2=opt2$bString, simList=fields4Sim,
CNOlist=CNOlist, indexList=indexNONCcutComp, boolUpdates=c(10,20),
divTime=10, lowerB=0.8, upperB=10)

writeScaffold(
	modelComprExpanded=ModelCutCompressExpand,
	optimResT1=opt1,
	optimResT2=opt2,
	modelOriginal=Model,
	CNOlist=CNOlist
)
	
writeNetwork(
	modelOriginal=Model,
	modelComprExpanded=ModelCutCompressExpand,
	optimResT1=opt1,
	optimResT2=opt2,
	CNOlist=CNOlist
)
