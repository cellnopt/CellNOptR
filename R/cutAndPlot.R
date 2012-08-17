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

cutAndPlot <- function(CNOlist, model, bStrings=NULL, plotPDF=FALSE, tag=NULL, show=TRUE)
{

    # bitStrings must be a list of bitString (T1, T2, ...TN)
    tPt = length(bStrings)+1

    simList <- prep4sim(model)
    indexList <- indexFinder(CNOlist=CNOlist,model=model)

    # if tPt nothing to plot
    if (tPt==1){
        stop("noting to do with time data at time 0")
    }

    # if tPt=2 (default), call cutAndPlotResultsT1
    if (tPt == 2){
       cutAndPlotResultsT1(CNOlist=CNOlist, model=model, bString=bStrings[[1]], simList=simList, 
             indexList=indexList, plotPDF=plotPDF, tag=tag, show=show)
    }

    # if tPt=2 (default), call cutAndPlotResultsT2
    if (tPt==3){ # use TN variant now.
       #print("Entering cutAndPlotResultsT2")
       cutAndPlotResultsT2(CNOlist, model=model, bStringT1=bStrings[[1]], 
         bStringT2=bStrings[[2]],  simList=simList, indexList=indexList,  plotPDF=plotPDF, tag=tag) 
    }

    if (tPt>3){
       #print("Entering cutAndPlotResultsTN")
       cutAndPlotResultsTN(
         CNOlist=CNOlist,
         model=model,
         bStrings=bStrings,
         simList=simList,
         indexList=indexList,
         plotPDF=plotPDF,
         tag=tag) 
    }
    # if tPt=2 (default), call cutAndPlotResultsT2


}
