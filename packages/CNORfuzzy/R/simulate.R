#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI - Massachusetts Institute of Technology
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv2 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-2.0.html
#
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/software/cno
#
##############################################################################
simulate<-function(postRefThresh=NULL, allFinalMSEs, allRes)
{

# initialize some things
    ReferenceRes = allRes[[1]]
    CNOlist=ReferenceRes$paramsList$data

     if ((class(CNOlist)=="CNOlist")==FALSE){
        CNOlist = CellNOptR::CNOlist(CNOlist)
    }


    indexList<-indexFinder(CNOlist=CNOlist,model=ReferenceRes$processedModel,verbose=ReferenceRes$paramsList$verbose)

    SimResAll = array(NA,dim(CNOlist@signals[[1]])[[1]]*dim(CNOlist@signals[[1]])[[2]]*length(allRes))
    dim(SimResAll) = c(dim(CNOlist@signals[[1]])[[1]],dim(CNOlist@signals[[1]])[[2]],length(allRes))
    # for each result, pick the best refined model, simulate it, and save the simulation results
    if(ReferenceRes$paramsList$doRefinement){
        for (eachRes in 1:length(allRes)){
            currIX = max(which(allFinalMSEs[eachRes,]-allFinalMSEs[eachRes,2]<=postRefThresh))
            currModel = allRes[[eachRes]]$redRef[[currIX]]$refModel$refinedModel
            currSimList = allRes[[eachRes]]$redRef[[currIX]]$refModel$refinedSimList
            SimCurr <- simFuzzyT1(CNOlist=CNOlist,currModel,currSimList)
            SimResAll[,,eachRes]=SimCurr[,indexList$signals]
        }
    }
    else{
        for (eachRes in 1:length(allRes)){
            currModel = allRes[[eachRes]]$unRef$model
            currSimList = allRes[[eachRes]]$unRef$simList
            SimCurr <- simFuzzyT1(CNOlist=CNOlist,currModel,currSimList)
            SimResAll[,,eachRes]=SimCurr[,indexList$signals]
        }
    }

# average the simulation results
    SimResMean = matrix(NA,dim(CNOlist@signals[[1]])[[1]],dim(CNOlist@signals[[1]])[[2]])
    for (eachSignal in 1:dim(CNOlist@signals[[1]])[[2]]){
        if (length(allRes)>1){
            SimResMean[,eachSignal]=apply(SimResAll[,eachSignal,],1,mean)
        }
        else{
            # if only one element in allRes, mean is useless, just copy the
            # data.
            SimResMean[,eachSignal]=SimResAll[,eachSignal,]
        }
    }

    simResults<-list(t0=matrix(data=0,nrow=dim(SimResAll)[1],ncol=dim(SimResAll)[2]),t1=SimResMean)

    return(simResults)
}

