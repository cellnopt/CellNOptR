#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv2 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-2.0.html
#
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/software.html
#
##############################################################################
# $Id$

cutAndPlotResultsT2 <-function(
	Model,
	bStringT1,
	bStringT2,
	SimList,
	CNOlist,
	indexList,
	plotPDF=FALSE, 
	tag=NULL,
	tPt=CNOlist$timeSignals[2:3]) {

	# simulate T1
	# prepare the model (i.e. cut)
		
	Modelcut <- Model
	Modelcut$interMat <- Modelcut$interMat[,as.logical(bStringT1)]
	Modelcut$notMat <- Modelcut$notMat[,as.logical(bStringT1)]
	Modelcut$reacID <- Modelcut$reacID[as.logical(bStringT1)]

    SimListCut <- cutSimList(SimList, bStringT1)
	
	# simulate
	SimT1 <- simulatorT1(
		CNOlist=CNOlist,
		Model=Modelcut,
		SimList=SimListCut,
		indexList=indexList)
	SimResT1 <- as.matrix(SimT1[,indexList$signals])
	
	# simulate T2
	
	# prepare the model
	bitString2 <- bStringT1
	bitString2[which(bStringT1 == 0)] <- bStringT2
	BStimes <- bStringT1
	BStimes[which(bStringT1 == 0)] <- bStringT2*2
	
	Modelcut <- Model
	Modelcut$interMat <- Modelcut$interMat[,as.logical(bitString2)]
	Modelcut$notMat <- Modelcut$notMat[,as.logical(bitString2)]
	Modelcut$reacID <- Modelcut$reacID[as.logical(bitString2)]
	Modelcut$times <- BStimes[which(BStimes != 0)]
	
    SimListCut <- cutSimList(SimList, bitString2)
	
	# simulate
	SimT2 <- simulatorT2(
		SimResultst1=SimT1,
		CNOlist=CNOlist,
		Model=Modelcut,
		SimList=SimListCut,
		indexList=indexList)
	SimResT2 <- SimT2[,indexList$signals]
	
	# put it all together

    Sim0 <- simulatorT0(CNOlist=CNOlist,Model=Modelcut,SimList=SimListCut,indexList=indexList)
    SimResT0 <- as.matrix(Sim0[,indexList$signals])
	
	SimResults <- list(
		t0=SimResT0,
		t1=SimResT1,
		t2=SimResT2)

		
	plotOptimResultsPan(
		SimResults=SimResults,
		CNOlist=CNOlist,
		formalism="ss2",
		tPt=tPt
	)
		
	if(plotPDF == TRUE) {
		if(is.null(tag)) {
			filename <- paste(deparse(substitute(Model)),"SimResultsT1T2.pdf",sep="")
		}
		else {
			filename<-paste(tag, "SimResultsT1T2.pdf", sep="_")
		}

		plotOptimResultsPan(
			SimResults=SimResults,
			CNOlist=CNOlist,
			formalism="ss2",
			tPt=tPt,
			pdfFileName=filename,
			pdf=TRUE
		)
	}		
}

