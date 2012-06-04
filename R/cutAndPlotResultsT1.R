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

cutAndPlotResultsT1 <- function(
	Model,
	bString,
	SimList,
	CNOlist,
	indexList,
	plotPDF=FALSE,
    tag=NULL,
    show=TRUE,
    tPt=CNOlist$timeSignals[2]) {

	Modelcut <- Model
	Modelcut$interMat <- Modelcut$interMat[,as.logical(bString)]
	Modelcut$notMat <- Modelcut$notMat[,as.logical(bString)]
	Modelcut$reacID <- Modelcut$reacID[as.logical(bString)]

	SimListCut<-cutSimList(SimList,bString)

	Sim <- simulatorT1(CNOlist=CNOlist,Model=Modelcut,SimList=SimListCut,indexList=indexList)
	SimRes <- as.matrix(Sim[,indexList$signals])

    # former code when t0 was not taken into account (everything set to zero)
	#SimResults <- list(t0=matrix(data=0,nrow=dim(SimRes)[1],ncol=dim(SimRes)[2]),t1=SimRes)

	# new code
	Sim0 <- simulatorT0(CNOlist=CNOlist,Model=Modelcut,SimList=SimListCut,indexList=indexList)
	SimRes0 <- as.matrix(Sim0[,indexList$signals])
	SimResults <- list(t0=SimRes0,t1=SimRes)
	dev.new()
    if(show==TRUE) {
    	plotOptimResultsPan(
			SimResults=SimResults,
			CNOlist=CNOlist,
			formalism="ss1",
			tPt=tPt
		)
	}
	if(plotPDF == TRUE) {
		if(is.null(tag)) {
			filename <- paste(deparse(substitute(Model)),"SimResultsT1.pdf",sep="")
        } else {
			filename <- paste(tag,"SimResultsT1.pdf",sep="_")
        }
		plotOptimResultsPan(
			SimResults=SimResults,
			CNOlist=CNOlist,
			pdf=TRUE,
			formalism="ss1",
			pdfFileName=filename,
			tPt=tPt
		)
	}
}

