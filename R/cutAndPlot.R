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

cutAndPlot <- function(model, CNOlist, bStringT1, bStringT2=NULL, 
    plotPDF=FALSE, tag=NULL, show=TRUE, tPt=2)
{

    simList <- prep4sim(model)
    indexList <- indexFinder(CNOlist=CNOlist,model=model)

    # if tPt nothing to plot
    if (tPt==1){
        stop("noting to do with time data at time 0")
    }

    # if tPt=2 (default), call cutAndPlotResultsT1
    if (tPt == 2){
       cutAndPlotResultsT1(model=model, bString=bStringT1, simList=simList, 
            CNOlist=CNOlist, indexList=indexList, plotPDF=plotPDF, tag=tag, 
            show=show)
    }
    # if tPt=2 (default), call cutAndPlotResultsT2
    if (tPt==3){
       cutAndPlotResultsT2(
         model=model,
         bStringT1=bStringT1,
         bStringT2=bStringT2,
         simList=simList,
         CNOlist=CNOlist,
         indexList=indexList,
         plotPDF=plotPDF,
         tag=tag) 
    }

    if (tPt>=3){
       cutAndPlotResultsTN(
         model=model,
         bStringT1=bStringT1,
         bStringT2=bStringT2,
         simList=simList,
         CNOlist=CNOlist,
         indexList=indexList,
         plotPDF=plotPDF,
         tag=tag) 
    }
    # if tPt=2 (default), call cutAndPlotResultsT2


}
