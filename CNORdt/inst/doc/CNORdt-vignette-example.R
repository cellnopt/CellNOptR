# check
library(CellNOptR)

# for later
library(CNORdt)

# check data names / documentation
data(CNOlistPB; package="CNORdt")
data(ToyModelPB; package="CNORdt")

load("data/CNOlistPB.RData")
load("data/modelPB.RData")

# pre-process model

model = preprocessing(CNOlistPB; modelPB)

# optimise
initBstring <- rep(1; length(model$reacID))

source("R/gaBinaryDT.R")
source("R/getFitDT.R")
source("R/simulatorDT.R")
dyn.load("src/simulatorDT.so")

opt1 <- gaBinaryDT(CNOlist=CNOlist; model=ModelCutCompressExpand;
simList=fields4Sim; indexList=indexNONCcutComp; initBstring=initBstring;
verbose=TRUE; boolUpdates=c(10;20); divTime=10; MaxTime=100; lowerB=0.8; upperB=10)

cutAndPlotResultsTimeScale(model=ModelCutCompressExpand; bString=opt1$bString; simList=fields4Sim;
CNOlist=CNOlist; indexList=indexNONCcutComp; boolUpdates=c(10;20); divTime=10; lowerB=0.8; upperB=10)

dataStartPoint = cutModel(model=ModelCutCompressExpand;
SimList=fields4Sim; bitString=opt1$bString)
simT1 = simulatorTimeScale(CNOlist=CNOlist;
model=dataStartPoint[[1]]; simList=dataStartPoint[[2]];
indexList=indexNONCcutComp; boolUpdates=10)

opt2 <- gaBinaryTimeScaleT2(CNOlist=CNOlist;
model=ModelCutCompressExpand; simList=fields4Sim;
indexList=indexNONCcutComp; bStringT1=opt1$bString;
simResT1=simT1; verbose=TRUE; boolUpdates=c(10;20);
divTime=10; lowerB=0.8; upperB=10)

cutAndPlotResultsTimeScaleT2(model=ModelCutCompressExpand;
bStringT1=opt1$bString; bStringT2=opt2$bString; simList=fields4Sim;
CNOlist=CNOlist; indexList=indexNONCcutComp; boolUpdates=c(10;20);
divTime=10; lowerB=0.8; upperB=10)

writeScaffold(
	modelComprExpanded=ModelCutCompressExpand;
	optimResT1=opt1;
	optimResT2=opt2;
	modelOriginal=Model;
	CNOlist=CNOlist
)
	
writeNetwork(
	modelOriginal=Model;
	modelComprExpanded=ModelCutCompressExpand;
	optimResT1=opt1;
	optimResT2=opt2;
	CNOlist=CNOlist
)
