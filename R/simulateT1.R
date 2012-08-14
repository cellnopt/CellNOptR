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
simulateT1<-function(CNOlist, model, bStringT1, simList=NULL, indexList=NULL){
    # simList and indexList are computed inside this function. 
    # However, for back-compatibility, we keep the arguments so that if
    # provided, we can still use them.


    # cut the model
    modelCut <- cutModel(model, bStringT1)
    if (is.NULL(simList)==TRUE){
        simList = prep4sim(model)
    }
    if (is.NULL(indexList)==TRUE){
        indexList = indexFinder(CNOlist, model)
    }

    # cut the model
    newSimList = cutSimList(simList, bStringT1)

    # compute the results
    simRes <- simulatorT1(CNOlist=CNOlist, model=modelCut, simList=newSimList, 
        indexList=indexList)

    return(simRes)
}
