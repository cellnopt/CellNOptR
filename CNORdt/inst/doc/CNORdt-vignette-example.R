library(CellNOptR)
library(CNORdt)

# data
data(CNOlistPB, package="CNORdt")
data(modelPB, package="CNORdt")

# pre-process model
model = preprocessing(CNOlistPB, modelPB)

# optimise
initBstring <- rep(1, length(model$reacID))

opt1 <- gaBinaryDT(CNOlist=CNOlistPB, model=model, initBstring=initBstring,
verbose=TRUE, boolUpdates=20, maxTime=60, lowerB=0.8, upperB=10)

cutAndPlotResultsDT(
	model=model,
	CNOlist=CNOlistPB,
	bString=opt1$bString,
	plotPDF=FALSE,
	boolUpdates=10,
	lowerB=0.8,
	upperB=10
)

writeScaffold(
	modelComprExpanded=model,
	optimResT1=opt1,
	optimResT2=NA,
	modelOriginal=modelPB,
	CNOlist=CNOlistPB
)

writeNetwork(
	modelOriginal=modelPB,
	modelComprExpanded=model,
	optimResT1=opt1,
	CNOlist=CNOlistPB
)
