library(CellNOptR)

# for later
library(CNORdt)

# check data names / documentation
data(CNOlistPB, package="CNORdt")
data(ToyModelPB, package="CNORdt")

load("data/CNOlistPB.RData")
load("data/modelPB.RData")

# pre-process model
model = preprocessing(CNOlistPB, modelPB)

# optimise
initBstring <- rep(1, length(model$reacID))

source("R/gaBinaryDT.R")
source("R/getFitDT.R")
source("R/simulatorDT.R")
source("R/computeScoreDT.R")
dyn.load("src/simulatorDT.so")


opt1 <- gaBinaryDT(CNOlist=CNOlistPB, model=model, initBstring=initBstring,
verbose=TRUE, boolUpdates=10, maxTime=30, lowerB=0.8, upperB=10)

source("R/cutAndPlotResultsDT.R")
cutAndPlotResultsDT(
	model=model,
	CNOlist=CNOlistPB,
	bString=opt1$bString,
	plotPDF=FALSE,boolUpdates=10, lowerB=0.8, upperB=10)

