# This file is part of the CNO software
# 
# Copyright (c) 2011-2013 - EBI
# 
# File author(s): CNO developers (cno-dev@ebi.ac.uk)
# 
# Distributed under the GPLv3 License.  See accompanying file LICENSE.txt or copy at
# http://www.gnu.org/licenses/gpl-3.0.html
# 
# CNO website: http://www.cellnopt.org
# 
# $Id$

computeScoreDT <- function(CNOlist, model, bString, simList = NULL, indexList = NULL, sizeFac = 1e-04, NAFac = 1, boolUpdates, 
    lowerB = lowerB, upperB = upperB) {
    
    # simList and indexList are computed inside this function
    # for back-compatibility, we keep the arguments so
    # that if provided, we can still use them
    if (is.null(simList) == TRUE) {
        simList = prep4sim(model)
    }
    if (is.null(indexList) == TRUE) {
        indexList = indexFinder(CNOlist, model)
    }
    
        if ((class(CNOlist) == "CNOlist") == FALSE) {
        CNOlist = CellNOptR::CNOlist(CNOlist)
    }
    
    modelCut = cutModel(model, bString)
    simListCut <- cutSimList(simList, bString)
    
    # compute the simulated results
    simResultsT0 <- simulatorT0(CNOlist = CNOlist, model = modelCut, simList = simListCut, indexList = indexList)
    
    simResults <- simulatorDT(CNOlist = CNOlist, model = modelCut, simList = simListCut, indices = indexList, boolUpdates = boolUpdates, 
        prevSim = simResultsT0)
    simResults = convert2array(simResults, dim(CNOlist@signals[[1]])[1], length(model$namesSpecies), boolUpdates)
    
    # compute the score
    optimResults <- getFitDT(simResults = simResults, CNOlist = CNOlist, model = modelCut, indexList = indexList, sizeFac = sizeFac, 
        NAFac = NAFac, nInTot = length(which(model$interMat == -1)), boolUpdates, lowerB = lowerB, upperB = upperB)   

    nDataP <- sum(!is.na(unlist(CNOlist@signals)))
    optimResults$score <- optimResults$score/nDataP
       
    return(optimResults$score)
}
