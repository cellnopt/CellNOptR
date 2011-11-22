cutAndPlotResultsT2 <-function(
	Model,
	bStringT1,
	bStringT2,
	SimList,
	CNOlist,
	indexList,
	plotPDF=FALSE, 
	tag=NULL){
	

	
	#Simulate T1
		#prepare the model (i.e. cut)
		
	Modelcut<-Model
	Modelcut$interMat<-Modelcut$interMat[,as.logical(bStringT1)]
	Modelcut$notMat<-Modelcut$notMat[,as.logical(bStringT1)]
	Modelcut$reacID<-Modelcut$reacID[as.logical(bStringT1)]
	
	SimListcut<-SimList
	SimListcut$finalCube<-SimListcut$finalCube[as.logical(bStringT1),]
	SimListcut$ixNeg<-SimListcut$ixNeg[as.logical(bStringT1),]
	SimListcut$ignoreCube<-SimListcut$ignoreCube[as.logical(bStringT1),]
	SimListcut$maxIx<-SimListcut$maxIx[as.logical(bStringT1)]
	
		#simulate
		
	SimT1<-simulatorT1(
		CNOlist=CNOlist,
		Model=Modelcut,
		SimList=SimListcut,
		indexList=indexList)
	SimResT1<-SimT1[,indexList$signals]
	
	#Simulate T2
	
		#Prepare the model
	bitString2<-bStringT1
	bitString2[which(bStringT1 == 0)]<-bStringT2
	BStimes<-bStringT1
	BStimes[which(bStringT1 == 0)]<-bStringT2*2
	
	Modelcut<-Model
	Modelcut$interMat<-Modelcut$interMat[,as.logical(bitString2)]
	Modelcut$notMat<-Modelcut$notMat[,as.logical(bitString2)]
	Modelcut$reacID<-Modelcut$reacID[as.logical(bitString2)]
	Modelcut$times<-BStimes[which(BStimes != 0)]
	
	SimListcut<-SimList
	SimListcut$finalCube<-SimListcut$finalCube[as.logical(bitString2),]
	SimListcut$ixNeg<-SimListcut$ixNeg[as.logical(bitString2),]
	SimListcut$ignoreCube<-SimListcut$ignoreCube[as.logical(bitString2),]
	SimListcut$maxIx<-SimListcut$maxIx[as.logical(bitString2)]
	
		#Simulate
	SimT2<-simulatorT2(
		SimResultst1=SimT1,
		CNOlist=CNOlist,
		Model=Modelcut,
		SimList=SimListcut,
		indexList=indexList)
	SimResT2<-SimT2[,indexList$signals]
	
	#Put it all together
	
	SimResults<-list(
		t0=matrix(data=0,nrow=dim(SimResT1)[1],ncol=dim(SimResT1)[2]),
		t1=SimResT1,
		t2=SimResT2)
		
	expResults<-list(
		t0=CNOlist$valueSignals[[1]],
		t1=CNOlist$valueSignals[[2]],
		t2=CNOlist$valueSignals[[3]])
		
	plotOptimResults(
		SimResults=SimResults,
		expResults=expResults,
		times=CNOlist$timeSignals[1:3],
		namesCues=CNOlist$namesCues,
		namesSignals=CNOlist$namesSignals,
		valueCues=CNOlist$valueCues)
		
	if(plotPDF == TRUE){
		if ( is.null(tag)){
			filename <- paste(deparse(substitute(Model)),"SimResultsT1T2.pdf",sep="")
		}
		else{
            filename<-paste(tag, "SimResultsT1T2.pdf", sep="_")
		}

		plotOptimResultsPDF(
			SimResults=SimResults,
			expResults=expResults,
			times=CNOlist$timeSignals[1:3],
			filename=filename,
			namesCues=CNOlist$namesCues,
			namesSignals=CNOlist$namesSignals,
			valueCues=CNOlist$valueCues)
		}
		
	}

