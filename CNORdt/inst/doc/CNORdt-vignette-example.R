library(CellNOptR)
library(CNORdt)

data(CNOlistPB, package="CNORdt")
data(ToyModelPB, package="CNORdt")

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
resECNOlist <- residualError(CNOlist)
fields4Sim <- prep4Sim(Model=ModelCutCompressExpand)
initBstring <- rep(1, length(ModelCutCompressExpand$reacID))

opt1 <- gaBinaryTimeScale(CNOlist=CNOlist, Model=ModelCutCompressExpand,
SimList=fields4Sim, indexList=indexNONCcutComp, initBstring=initBstring,
verbose=TRUE, boolUpdates=c(10,20), divTime=10, MaxTime=100, lowerB=0.8, upperB=10)

cutAndPlotResultsTimeScale(Model=ModelCutCompressExpand, bString=opt1$bString, SimList=fields4Sim,
CNOlist=CNOlist, indexList=indexNONCcutComp, boolUpdates=c(10,20), divTime=10, lowerB=0.8, upperB=10)

dataStartPoint = cutModel(Model=ModelCutCompressExpand,
SimList=fields4Sim, bitString=opt1$bString)
simT1 = simulatorTimeScale(CNOlist=CNOlist,
Model=dataStartPoint[[1]], SimList=dataStartPoint[[2]],
indexList=indexNONCcutComp, boolUpdates=10)

opt2 <- gaBinaryTimeScaleT2(CNOlist=CNOlist,
Model=ModelCutCompressExpand, SimList=fields4Sim,
indexList=indexNONCcutComp, bStringT1=opt1$bString,
SimResT1=simT1, verbose=TRUE, boolUpdates=c(10,20),
divTime=10, lowerB=0.8, upperB=10)

cutAndPlotResultsTimeScaleT2(Model=ModelCutCompressExpand,
bStringT1=opt1$bString, bStringT2=opt2$bString, SimList=fields4Sim,
CNOlist=CNOlist, indexList=indexNONCcutComp, boolUpdates=c(10,20),
divTime=10, lowerB=0.8, upperB=10)

writeScaffold(
	ModelComprExpanded=ModelCutCompressExpand,
	optimResT1=opt1,
	optimResT2=opt2,
	ModelOriginal=Model,
	CNOlist=CNOlist
)
	
writeNetwork(
	ModelOriginal=Model,
	ModelComprExpanded=ModelCutCompressExpand,
	optimResT1=opt1,
	optimResT2=opt2,
	CNOlist=CNOlist
)
