cutAndPlotResultsTimeScale <- function(Model, bString, SimList, CNOlist, indexList, boolUpdates=boolUpdates, divTime=NULL) {
	
	Modelcut <- Model
	Modelcut$interMat <- Modelcut$interMat[,as.logical(bString)]
	Modelcut$notMat <- Modelcut$notMat[,as.logical(bString)]
	Modelcut$reacID <- Modelcut$reacID[as.logical(bString)]
	SimListcut <- SimList
	SimListcut$finalCube <- SimListcut$finalCube[as.logical(bString),]
	SimListcut$ixNeg <- SimListcut$ixNeg[as.logical(bString),]
	SimListcut$ignoreCube <- SimListcut$ignoreCube[as.logical(bString),]
	SimListcut$maxIx <- SimListcut$maxIx[as.logical(bString)]
	
	boolUpdates = boolUpdates[1]
	SimRes <- simulatorTimeScale(CNOlist=CNOlist, Model=Modelcut, SimList=SimListcut, indexList=indexList, boolUpdates=boolUpdates)
	SimRes = SimRes[,indexList$signals,]
	getFitData <- getFitTimeScale(SimList=SimListcut, CNOlist=CNOlist, Model=Modelcut, indexList=indexList, boolUpdates=boolUpdates, divTime=divTime)
	AU <- getFitData$estimate
	x.coords <- getFitData$x.coords
	expResults <- CNOlist$valueSignals
	
	plotOptimResultsTimeScale(SimResults=SimRes, yInterpol=getFitData$yInter, xCoords=x.coords, expResults=expResults, times=CNOlist$timeSignals, namesCues=CNOlist$namesCues, namesSignals=CNOlist$namesSignals, valueCues=CNOlist$valueCues)

}

