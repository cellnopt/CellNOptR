cutAndPlotResultsT1<-function(
	Model,
	bString,
	SimList,
	CNOlist,
	indexList,
	plotPDF=FALSE,
    tag=NULL,
    show=TRUE){
	
	Modelcut<-Model
	Modelcut$interMat<-Modelcut$interMat[,as.logical(bString)]
	Modelcut$notMat<-Modelcut$notMat[,as.logical(bString)]
	Modelcut$reacID<-Modelcut$reacID[as.logical(bString)]
	
	SimListcut<-SimList
	SimListcut$finalCube<-SimListcut$finalCube[as.logical(bString),]
	SimListcut$ixNeg<-SimListcut$ixNeg[as.logical(bString),]
	SimListcut$ignoreCube<-SimListcut$ignoreCube[as.logical(bString),]
	SimListcut$maxIx<-SimListcut$maxIx[as.logical(bString)]
	
	Sim<-simulatorT1(CNOlist=CNOlist,Model=Modelcut,SimList=SimListcut,indexList=indexList)
	
	SimRes<-Sim[,indexList$signals]
	SimResults<-list(t0=matrix(data=0,nrow=dim(SimRes)[1],ncol=dim(SimRes)[2]),t1=SimRes)
	expResults<-list(t0=CNOlist$valueSignals[[1]],t1=CNOlist$valueSignals[[2]])
	
    if (show == TRUE){
    	plotOptimResults(
	    	SimResults=SimResults,
		    expResults=expResults,
    		times=CNOlist$timeSignals[1:2],
	    	namesCues=CNOlist$namesCues,
		    namesSignals=CNOlist$namesSignals,
    		valueCues=CNOlist$valueCues)
	}
	if(plotPDF == TRUE){
        if ( is.null(tag)){
               filename<-paste(deparse(substitute(Model)), "SimResultsT1.pdf", sep="")
        }
        else{
            filename<-paste(tag, "SimResultsT1.pdf", sep="_")
        }
		plotOptimResultsPDF(
			SimResults=SimResults,
			expResults=expResults,
			times=CNOlist$timeSignals[1:2],
			filename=filename,
			namesCues=CNOlist$namesCues,
			namesSignals=CNOlist$namesSignals,
			valueCues=CNOlist$valueCues)
		}
	}

