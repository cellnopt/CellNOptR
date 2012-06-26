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

cutAndPlotResultsT1 <- function(model, bString, simList, CNOlist, indexList,
    plotPDF=FALSE, tag=NULL, show=TRUE, tPt=CNOlist$timeSignals[2])
{

    modelCut <- cutModel(model, bString)

    simListCut<-cutSimList(simList,bString)

    Sim <- simulatorT1(CNOlist=CNOlist,model=modelCut,simList=simListCut,indexList=indexList)
    simRes <- as.matrix(Sim[,indexList$signals])

    # former code when t0 was not taken into account (everything set to zero)
    #simResults <- list(t0=matrix(data=0,nrow=dim(simRes)[1],ncol=dim(simRes)[2]),t1=simRes)

    # new code
    Sim0 <- simulatorT0(CNOlist=CNOlist,model=modelCut,simList=simListCut,indexList=indexList)
    simRes0 <- as.matrix(Sim0[,indexList$signals])
    simResults <- list(t0=simRes0,t1=simRes)
    #dev.new()
    if(show==TRUE) {
        plotOptimResultsPan(
            simResults=simResults,
            CNOlist=CNOlist,
            formalism="ss1",
            tPt=tPt
        )
    }
    if(plotPDF == TRUE) {
        if(is.null(tag)) {
            filename <- paste(deparse(substitute(model)),"SimResultsT1.pdf",sep="")
        } else {
            filename <- paste(tag,"SimResultsT1.pdf",sep="_")
        }
        plotOptimResultsPan(
            simResults=simResults,
            CNOlist=CNOlist,
            pdf=TRUE,
            formalism="ss1",
            pdfFileName=filename,
            tPt=tPt
        )
    }
}

