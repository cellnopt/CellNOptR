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
interpretDiscreteGA = function(model,paramsList,intString,bitString=NULL){

    simList<-prep4simFuzzy(model=model,paramsList=paramsList, verbose=FALSE)



    #need to add cutting the type cube and adjusting numType1 and numType2
    # first determine what the bitstring should be if it isn not provided
    ModelCut<-model

    if (is.null(bitString)) {
        bitCube = matrix(data = 0,nrow=dim(simList$finalCube)[1],ncol=dim(simList$finalCube)[2])
        # currently not using params2train parameter in R version
        bitCube[simList$reshapeType1] = intString[1:simList$numType1]
        bitCube[simList$reshapeType2] = intString[(simList$numType1+1):length(intString)]
        maxString = apply(bitCube,1,max)
        bitString = maxString != 0
    } 
	else {
        bitString = bitString
    }


    # now cut model accordingly
    ModelCut$interMat<-ModelCut$interMat[,bitString]
    ModelCut$notMat<-ModelCut$notMat[,bitString]
    ModelCut$reacID<-ModelCut$reacID[bitString]
    type1Vals = intString[1:simList$numType1]
    type2Vals = intString[(simList$numType1+1):length(intString)]
    SimListCut<-simList   
    for (i in 1:dim(paramsList$type1Funs)[1]) {
        #Set Type 1
        SimListCut$kCube[simList$reshapeType1[type1Vals == i]] = paramsList$type1Funs[i,3]
        SimListCut$nCube[simList$reshapeType1[type1Vals == i]] = paramsList$type1Funs[i,2]
        SimListCut$gCube[simList$reshapeType1[type1Vals == i]] = paramsList$type1Funs[i,1]
        #Set Type 2
        SimListCut$kCube[simList$reshapeType2[type2Vals == i]] = paramsList$type2Funs[i,3]
        SimListCut$nCube[simList$reshapeType2[type2Vals == i]] = paramsList$type2Funs[i,2]
        SimListCut$gCube[simList$reshapeType2[type2Vals == i]] = paramsList$type2Funs[i,1]
    }

    SimListCut$finalCube<-SimListCut$finalCube[bitString,]
    SimListCut$ixNeg<-SimListCut$ixNeg[bitString,]
    SimListCut$ignoreCube<-SimListCut$ignoreCube[bitString,]
    SimListCut$gCube<-SimListCut$gCube[bitString,]
    SimListCut$nCube<-SimListCut$nCube[bitString,]
    SimListCut$kCube<-SimListCut$kCube[bitString,]
    SimListCut$typeCube<-SimListCut$typeCube[bitString,]
    SimListCut$maxIx<-SimListCut$maxIx[bitString]

    SimListCut$reshapeIx = which(!SimListCut$ignoreCube)
    SimListCut$reshapeType1 = which(SimListCut$typeCube==1)
    SimListCut$reshapeType2 = which(SimListCut$typeCube==2)
    SimListCut$numType1 = length(SimListCut$reshapeType1)
    SimListCut$numType2 = length(SimListCut$reshapeType2)

    #now further cut for redundant logic
    refInit = which(bitString)
    reacs2Cut = removeRedundant(ModelCut)
    cutBitString = bitString
    cutBitString[refInit[reacs2Cut]] = FALSE

    # make a new cut model that does not contain redundant logic
    ModelCutCut = model
    ModelCutCut$interMat<-ModelCutCut$interMat[,cutBitString]
    ModelCutCut$notMat<-ModelCutCut$notMat[,cutBitString]
    ModelCutCut$reacID<-ModelCutCut$reacID[cutBitString]	
    SimListCutCut = simList

    for (i in 1:dim(paramsList$type1Funs)[1]) {
	    #Set Type 1
        SimListCutCut$kCube[SimListCutCut$reshapeType1[type1Vals == i]] = paramsList$type1Funs[i,3];
        SimListCutCut$nCube[SimListCutCut$reshapeType1[type1Vals == i]] = paramsList$type1Funs[i,2];
        SimListCutCut$gCube[SimListCutCut$reshapeType1[type1Vals == i]] = paramsList$type1Funs[i,1];
    	#Set Type 2
        SimListCutCut$kCube[SimListCutCut$reshapeType2[type2Vals == i]] = paramsList$type2Funs[i,3];
        SimListCutCut$nCube[SimListCutCut$reshapeType2[type2Vals == i]] = paramsList$type2Funs[i,2];
        SimListCutCut$gCube[SimListCutCut$reshapeType2[type2Vals == i]] = paramsList$type2Funs[i,1];
    }
    SimListCutCut$finalCube<-SimListCutCut$finalCube[cutBitString,]
    SimListCutCut$ixNeg<-SimListCutCut$ixNeg[cutBitString,]
    SimListCutCut$ignoreCube<-SimListCutCut$ignoreCube[cutBitString,]
    SimListCutCut$gCube<-SimListCutCut$gCube[cutBitString,]
    SimListCutCut$nCube<-SimListCutCut$nCube[cutBitString,]
    SimListCutCut$kCube<-SimListCutCut$kCube[cutBitString,]
    SimListCutCut$typeCube<-SimListCutCut$typeCube[cutBitString,]
    SimListCutCut$maxIx<-SimListCutCut$maxIx[cutBitString]

    SimListCutCut$reshapeIx = which(!SimListCutCut$ignoreCube)
    SimListCutCut$reshapeType1 = which(SimListCutCut$typeCube==1)
    SimListCutCut$reshapeType2 = which(SimListCutCut$typeCube==2)
    SimListCutCut$numType1 = length(SimListCutCut$reshapeType1)
    SimListCutCut$numType2 = length(SimListCutCut$reshapeType2)

    # returns everything
    return(list(model = ModelCut,simList=SimListCut,bitString=bitString,cutModel = ModelCutCut,cutSimList=SimListCutCut,cutBitString = cutBitString))
}
