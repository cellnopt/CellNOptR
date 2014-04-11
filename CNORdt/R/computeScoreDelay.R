# This file is part of the CNO software
# 
# Copyright (c) 2011-2012 - EBI
# 
# File author(s): CNO developers (cno-dev@ebi.ac.uk)
# 
# Distributed under the GPLv2 License.  See accompanying file LICENSE.txt or copy at
# http://www.gnu.org/licenses/gpl-2.0.html
# 
# CNO website: http://www.ebi.ac.uk/saezrodriguez/software.html
# 
# $Id$

computeScoreDelay <- function(CNOlist, model, bString, simList = NULL, indexList = NULL, sizeFac = 1e-04, NAFac = 1, boolUpdates) {
    
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
   
    # compute the score
    optimResults <- getFitPause(CNOlist = CNOlist, model = modelCut, indexList = indexList, simList=simListCut, sizeFac = sizeFac, 
        NAFac = NAFac, nInTot = length(which(model$interMat == -1)), boolUpdates)   
	print(optimResults$score)
    nDataP <- sum(!is.na(unlist(CNOlist@signals)))
    optimResults$score <- optimResults$score/nDataP
       
    return(optimResults$score)
}
