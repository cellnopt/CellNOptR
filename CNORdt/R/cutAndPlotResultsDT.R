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
# $Id: $

cutAndPlotResultsDT <- function(model, bString, simList=NULL, CNOlist, indexList=NULL,
 plotPDF=FALSE, tag=NULL, plotParams=list(maxrow=10),boolUpdates=boolUpdates, lowerB=lowerB, upperB=upperB,sizeFac = 1e-04, NAFac = 1)
{
	
    if ((class(CNOlist)=="CNOlist")==FALSE){
        CNOlist = CellNOptR::CNOlist(CNOlist)
    } 
    if (is.null(simList)==TRUE){
        simList = prep4sim(model)
    }
    if (is.null(indexList)==TRUE){
        indexList = indexFinder(CNOlist, model)
    }

    if ("maxrow" %in% names(plotParams) == FALSE){
        plotParams$maxrow = 10
    }

    # keep simList and indxList for back compatibility ?
    modelCut <- cutModel(model, bString)
    simListCut <- cutSimList(simList, bString)

    # t0
    Sim0 <- simulatorT0(CNOlist=CNOlist, model=modelCut, simList=simListCut, indexList=indexList)

    Sim <- simulatorDT(CNOlist=CNOlist, model=modelCut, simList=simListCut, indices=indexList, boolUpdates=boolUpdates, prevSim = Sim0) 
    
    simResults = convert2array(Sim, dim(CNOlist@signals[[1]])[1], length(model$namesSpecies), boolUpdates)

    optimResults <- getFitDT(simResults = simResults, CNOlist = CNOlist, model = modelCut, indexList = indexList, sizeFac = sizeFac, 
        NAFac = NAFac, nInTot = length(which(model$interMat == -1)), boolUpdates,  lowerB = lowerB, upperB = upperB)
	    simResults <- simResults[, indexList$signals, ]
  if(dim(CNOlist@signals[[1]])[1] == 1) {
    	simResults = array(simResults,
    	c(dim(CNOlist@signals[[1]])[1], length(indexList$signals), boolUpdates)
    	)
    }
    dim1 = dim(CNOlist@signals[[1]])[1]
    dim2 = dim(CNOlist@signals[[1]])[2]

    CNOlistSet = list()
    simResultsSet = list()

    if(dim1 > plotParams$maxrow) { #|| dim2 > 10) {

        par1 = ceiling(dim1/plotParams$maxrow)
        div1 = ceiling(dim1/par1)
        #par2 = ceiling(dim2/plotParams$maxrow)
        #div2 = ceiling(dim2/par2)

        count1 = 1
        for(a in 1:par1) {
            CNOdiv = CNOlist
            simDiv = simResults
            finalN = div1 * a
            if(finalN > dim1) {finalN = dim1}
            CNOdiv@cues = CNOdiv@cues[count1:finalN,]
            CNOdiv@stimuli = CNOdiv@stimuli[count1:finalN,]
            CNOdiv@inhibitors = CNOdiv@inhibitors[count1:finalN,]
            for(b in 1:length(CNOdiv@signals)) {
                CNOdiv@signals[[b]] = CNOdiv@signals[[b]][count1:finalN,]
            }
            for(d in 1:length(simDiv)) {
                simDiv[[d]] = simDiv[[d]][count1:finalN,]
            }
            count1 = count1 + div1
            CNOlistSet = c(CNOlistSet, list(CNOdiv))
            simResultsSet = c(simResultsSet, list(simDiv))
        }
    } else {

        CNOlistSet = list(CNOlist)
        simResultsSet = list(simResults)
    }

    outputFilenames = list()
    for(f in 1:length(CNOlistSet)) {

        plotOptimResultsPan(
            simResults=simResultsSet[[f]],
            CNOlist=CNOlistSet[[f]],
        	yInterpol=optimResults$yInter,
			xCoords=optimResults$xCoords,
			formalism="dt",
			tPt=CNOlist@timepoints,
            plotParams=plotParams
            )

        if(plotPDF == TRUE) {
            if(is.null(tag)) {
                filename <- paste("SimResultsT1_", f, ".pdf", sep="")
            } else {
                filename <- paste(tag,"SimResultsT1",f,".pdf",sep="_")
            }
            plotOptimResultsPan(
                 simResults=simResultsSet[[f]],
            CNOlist=CNOlistSet[[f]],
        	yInterpol=optimResults$yInter,
			xCoords=optimResults$xCoords,
			formalism="dt",
			tPt=CNOlist@timepoints,
            plotParams=plotParams,
                            pdf=TRUE,
                pdfFileName=filename

            )
            outputFilenames[[f]] = filename
        }
    }
    return(outputFilenames)
}
