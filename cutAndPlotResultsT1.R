cutAndPlotResultsT1<-function(
	Model,
	bString,
	SimList,
	CNOlist,
	indexList,
	plotPDF=FALSE){
	
	Modelcut<-Model
	Modelcut$interMat<-as.matrix(Modelcut$interMat[,as.logical(bString)])
	Modelcut$notMat<-as.matrix(Modelcut$notMat[,as.logical(bString)])
	Modelcut$reacID<-as.vector(Modelcut$reacID[as.logical(bString)])
	SimListcut<-SimList
	SimListcut$finalCube<-as.matrix(SimListcut$finalCube[as.logical(bString),])
	SimListcut$ixNeg<-as.matrix(SimListcut$ixNeg[as.logical(bString),])
	SimListcut$ignoreCube<-as.matrix(SimListcut$ignoreCube[as.logical(bString),])
	SimListcut$maxIx<-as.vector(SimListcut$maxIx[as.logical(bString)])
		
	Sim<-simulatorT1(CNOlist=CNOlist,Model=Modelcut,SimList=SimListcut,indexList=indexList)
	
	SimRes<-Sim[,indexList$signals]
	SimResults<-list(t0=matrix(data=0,nrow=dim(SimRes)[1],ncol=dim(SimRes)[2]),t1=SimRes)
	expResults<-list(t0=CNOlist$valueSignals[[1]],t1=CNOlist$valueSignals[[2]])
	
	plotOptimResults(
		SimResults=SimResults,
		expResults=expResults,
		times=CNOlist$timeSignals[1:2],
		namesCues=CNOlist$namesCues,
		namesSignals=CNOlist$namesSignals,
		valueCues=CNOlist$valueCues)
		
	if(plotPDF == TRUE){
		plotOptimResultsPDF(
			SimResults=SimResults,
			expResults=expResults,
			times=CNOlist$timeSignals[1:2],
			filename=paste(deparse(substitute(Model)),"SimResultsT1.pdf",sep=""),
			namesCues=CNOlist$namesCues,
			namesSignals=CNOlist$namesSignals,
			valueCues=CNOlist$valueCues)
		}
	}

