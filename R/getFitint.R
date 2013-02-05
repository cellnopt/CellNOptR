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
getFitint<-function(
	SimResults,
	CNOlist,
	model,
	indexList,
	timePoint=c("t1","t2"),
	sizeFac=0.0001,
	LinkPen,
	NAFac=1,
	nInTot){
	
	SimResults<-SimResults[,indexList$signals]
	
	if(timePoint == "t1") tPt<-2
	if(timePoint == "t2") tPt<-3
	
	Diff<-SimResults-CNOlist$valueSignals[[tPt]]
	r<-Diff^2
	
	deviationPen<-sum(r[!is.na(r)])
	
	NAPen<-NAFac*length(which(is.na(SimResults)))
	
	nDataPts<-dim(CNOlist$valueSignals[[tPt]])[1]*dim(CNOlist$valueSignals[[tPt]])[2]
	
	#nInputs<-length(which(model$interMat == -1))
			
	#nInputsIntegr<-length(which(model$interMat[,model$indexIntegr] == -1))
				
	#sizePen<-((nDataPts*sizeFac)/nInTot)*(nInputs+(integrFac-1)*nInputsIntegr)
	sizePen<-((nDataPts*sizeFac)/nInTot)*LinkPen
		
	score<-deviationPen+NAPen+sizePen
	
	return(score)
	
	}

