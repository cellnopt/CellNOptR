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
#
# File author(s): T. Cokelaer based on gaDiscrete from M.K. Morris
computeScoreFuzzy <- function(CNOlist, model, simList=NULL, indexList=NULL,
    paramsList, intString=NULL, sizeFac=0.0001,NAFac=1){
    #initialise

    if ((class(CNOlist)=="CNOlist")==FALSE){
        CNOlist = CellNOptR::CNOlist(CNOlist)
    } 


    if (is.null(indexList)==TRUE){
        indexList = indexFinder(CNOlist, model, verbose=FALSE)
    }
    if (is.null(simList)==TRUE){
        simList = prep4simFuzzy(model, paramsList, verbose=FALSE)
    }

    if (is.null(intString)==TRUE){
        intString <- (sample.int(dim(paramsList$type2Funs)[1],
            (simList$numType1+simList$numType2),replace=TRUE)) - 1
    }


    nInTot = length(which(model$interMat==-1))

    #cut the model according to bitstring
    ModelCut<-model
    bitCube = matrix(data = 0,nrow=dim(simList$finalCube)[1],ncol=dim(simList$finalCube)[2])
    # currently not using params2train parameter in R version
    bitCube[simList$reshapeType1] = intString[1:simList$numType1]
    bitCube[simList$reshapeType2] = intString[(simList$numType1+1):length(intString)]
    maxString = apply(bitCube,1,max)
    bitString = maxString != 0
    ModelCut$interMat<-ModelCut$interMat[,bitString]
    ModelCut$notMat<-ModelCut$notMat[,bitString]
    ModelCut$reacID<-ModelCut$reacID[bitString]

    type1Vals = intString[1:simList$numType1];
    type2Vals = intString[(simList$numType1+1):length(intString)];

    simListCut<-simList

    for (i in 1:dim(paramsList$type1Funs)[1]) {
        #Set Type 1
        simListCut$kCube[simList$reshapeType1[type1Vals == i]] = paramsList$type1Funs[i,3];
        simListCut$nCube[simList$reshapeType1[type1Vals == i]] = paramsList$type1Funs[i,2];
        simListCut$gCube[simList$reshapeType1[type1Vals == i]] = paramsList$type1Funs[i,1];
        #Set Type 2
        simListCut$kCube[simList$reshapeType2[type2Vals == i]] = paramsList$type2Funs[i,3];
        simListCut$nCube[simList$reshapeType2[type2Vals == i]] = paramsList$type2Funs[i,2];
        simListCut$gCube[simList$reshapeType2[type2Vals == i]] = paramsList$type2Funs[i,1];
    }
    simListCut$finalCube<-simListCut$finalCube[bitString,]
    simListCut$ixNeg<-simListCut$ixNeg[bitString,]
    simListCut$ignoreCube<-simListCut$ignoreCube[bitString,]
    simListCut$gCube<-simListCut$gCube[bitString,]
    simListCut$nCube<-simListCut$nCube[bitString,]
    simListCut$kCube<-simListCut$kCube[bitString,]
    simListCut$maxIx<-simListCut$maxIx[bitString]

    #compute the simulated results
    SimResults<-simFuzzyT1(CNOlist=CNOlist,model=ModelCut,simList=simListCut)

    #Compute the score
    Score<-getFit(simResults=SimResults,CNOlist=CNOlist,model=ModelCut,indexList=indexList,timePoint="t1",sizeFac=sizeFac,NAFac=NAFac,nInTot=nInTot)
    nDataP<-sum(!is.na(CNOlist@signals[[2]]))
    Score<-Score/nDataP

    return(Score)
}

