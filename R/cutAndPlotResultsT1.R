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


    SimListCut<-cutSimList(SimList, bString)

	
	Sim<-simulatorT1(CNOlist=CNOlist,Model=Modelcut,SimList=SimListCut,indexList=indexList)
	
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

