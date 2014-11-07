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
reduceFuzzy = function(firstCutOff, CNOlist, model, res, params)
{
 if ((class(CNOlist)=="CNOlist")==FALSE){
        CNOlist = CellNOptR::CNOlist(CNOlist)
    } 


    simList = prep4simFuzzy(model, params, verbose=FALSE)
    indexList<-indexFinder(CNOlist=CNOlist,model=model, verbose=FALSE)

    # [finalBits finalParams finalScores ActuallyRemove ActuallyReplace] =
    # CNOReduceFuzzyOptInside(FirstCutOff, CNOProject, Res)
    #
    # Reduces model "inside CNOGetModel" according to reduction threshold
    # (FirstCutOff).
    #
    # contributed by Melody K Morris Jan. 2011

    nDP = sum(!is.na(CNOlist@signals[[2]]))

    finalScores = array(NaN,4)
    interpModel = interpretDiscreteGA(model=model, intString=res$bString, paramsList=params)
    BaseSimRes = simFuzzyT1(CNOlist,interpModel$cutModel,interpModel$cutSimList)
    # tocheck: why nInTot = 1 ?
    BaseMSE = getFit(BaseSimRes,CNOlist,interpModel$cutModel,indexList, timePoint="t1",sizeFac=0,NAFac=1,nInTot=1)



    BaseMSE = BaseMSE/nDP
    BitString = interpModel$cutBitString
    finalScores[1] = BaseMSE

    ReacsToTest = which(BitString);
    Reacs2Remove = array(NA,0);
    for (y in 1:length(ReacsToTest)) {
        bitTry = BitString
        bitTry[ReacsToTest[y]] = FALSE
    ######
        currCutModel<-model
        currCutModel$interMat<-currCutModel$interMat[,bitTry]
        currCutModel$notMat<-currCutModel$notMat[,bitTry]
        currCutModel$reacID<-currCutModel$reacID[bitTry]
        type1Vals = res$bString[1:simList$numType1];
        type2Vals = res$bString[(simList$numType1+1):length(res$bString)];
        currSimList<-simList
        for (i in 1:dim(params$type1Funs)[1]) {
            #Set Type 1
            currSimList$kCube[simList$reshapeType1[type1Vals == i]] = params$type1Funs[i,3];
            currSimList$nCube[simList$reshapeType1[type1Vals == i]] = params$type1Funs[i,2];
            currSimList$gCube[simList$reshapeType1[type1Vals == i]] = params$type1Funs[i,1];
            #Set Type 2
                currSimList$kCube[simList$reshapeType2[type2Vals == i]] = params$type2Funs[i,3];
            currSimList$nCube[simList$reshapeType2[type2Vals == i]] = params$type2Funs[i,2];
            currSimList$gCube[simList$reshapeType2[type2Vals == i]] = params$type2Funs[i,1];
        }
        # now cut sim list
        currSimList$kCube = currSimList$kCube[bitTry,]
        currSimList$nCube = currSimList$nCube[bitTry,]
        currSimList$gCube = currSimList$gCube[bitTry,]
        currSimList$ixNeg = currSimList$ixNeg[bitTry,]
        currSimList$finalCube= currSimList$finalCube[bitTry,]
        currSimList$ignoreCube= currSimList$ignoreCube[bitTry,]
        currSimList$maxIx= currSimList$maxIx[bitTry]
        # simulate and score
        SimResults = simFuzzyT1(CNOlist,currCutModel,currSimList)
        testMSE = getFit(SimResults,CNOlist,currCutModel,indexList,timePoint="t1",sizeFac=0,NAFac=1,nInTot=1)
        testMSE = testMSE/nDP

        if (testMSE <= (BaseMSE + firstCutOff)) {
            Reacs2Remove = rbind(Reacs2Remove,ReacsToTest[y]);
        }
    }
    Reacs2Remove = unique(Reacs2Remove)

    #THIS PART computes the MSE after removal
    if (length(Reacs2Remove) > 0) {
        notBadBits = BitString
        notBadBits[Reacs2Remove] = FALSE
        currCutModel<-model
        currCutModel$interMat<-currCutModel$interMat[,notBadBits]
        currCutModel$notMat<-currCutModel$notMat[,notBadBits]
        currCutModel$reacID<-currCutModel$reacID[notBadBits]
        type1Vals = res$bString[1:simList$numType1];
        type2Vals = res$bString[(simList$numType1+1):length(res$bString)];
        currSimList<-simList
        for (i in 1:dim(params$type1Funs)[1]) {
            #Set Type 1
            currSimList$kCube[simList$reshapeType1[type1Vals == i]] = params$type1Funs[i,3];
            currSimList$nCube[simList$reshapeType1[type1Vals == i]] = params$type1Funs[i,2];
            currSimList$gCube[simList$reshapeType1[type1Vals == i]] = params$type1Funs[i,1];
            #Set Type 2
            currSimList$kCube[simList$reshapeType2[type2Vals == i]] = params$type2Funs[i,3];
            currSimList$nCube[simList$reshapeType2[type2Vals == i]] = params$type2Funs[i,2];
            currSimList$gCube[simList$reshapeType2[type2Vals == i]] = params$type2Funs[i,1];
        }
        # now cut sim list
        currSimList$kCube = currSimList$kCube[notBadBits,]
        currSimList$nCube = currSimList$nCube[notBadBits,]
        currSimList$gCube = currSimList$gCube[notBadBits,]
        currSimList$ixNeg = currSimList$ixNeg[notBadBits,]
        currSimList$finalCube= currSimList$finalCube[notBadBits,]
        currSimList$ignoreCube= currSimList$ignoreCube[notBadBits,]
        currSimList$maxIx= currSimList$maxIx[notBadBits]
    # simulate and calculate score
        SimResults = simFuzzyT1(CNOlist,currCutModel,currSimList)
        testMSE = getFit(SimResults,CNOlist,currCutModel,indexList, timePoint="t1",sizeFac=0,NAFac=1,nInTot=1)
        testMSE = testMSE/nDP
        finalScores[2] = testMSE;
    }

    # Replacement initialization

    Param1Const = simList$typeCube
    Param1Const[simList$typeCube==2] = NaN
    Param1Const[is.nan(simList$typeCube)] = NaN
    Param1Const[simList$reshapeType1] = 1:simList$numType1

    Param2Const = simList$typeCube
    Param2Const[simList$typeCube==1] = NaN
    Param2Const[is.nan(simList$typeCube)] = NaN
    Param2Const[simList$reshapeType2] = (simList$numType1+1):(simList$numType1+simList$numType2)


    #THIS PART FINDS THE ONES TO CONSIDER REPLACING

    inputsReacs = apply(model$interMat == -1,2,sum)

    NewFamilyParams = array(NA,0)
    NewFamilyBits = array(NA,0)
    ReplacementVec = array(NA,0)
    OrGates = which(inputsReacs == 1)
    #species index of all or-ed species
    OrIns = array(NA,length(OrGates))
    OrOuts = array(NA,length(OrGates))
    for (eachOr in 1:length(OrGates)){
        OrIns[eachOr] = which(model$interMat[,OrGates[eachOr]]==-1)
        OrOuts[eachOr] = which(model$interMat[,OrGates[eachOr]]==1)
    }

    AndGates = which(inputsReacs > 1 & BitString)
    #donot consider gates that will be removed anyway
    if (length(Reacs2Remove) > 0){
        for (eachRem in 1:length(Reacs2Remove)){
            AndGates = AndGates[AndGates!=Reacs2Remove[eachRem]]
            }
    }

    if (!any(AndGates)) {
        ReplacementVec = array(NA,0)
        ActuallyReplace = array(NA,0)
    }
    else {
        for (t in 1:length(AndGates)) {
            #species index to look for to compare to this specific gate
            inputsTest = which(model$interMat[,AndGates[t]]==-1)
            currOutput = which(model$interMat[,AndGates[t]]==1)
            for (o in 1:length(inputsTest)) {
                currBits = BitString;
                currBits[AndGates[t]] = FALSE;
                currBits[OrGates[inputsTest[o]==OrIns]] = TRUE
                currCutModel = model
                currSimList = simList
                for (i in 1:dim(params$type1Funs)[1]) {
                    #Set Type 1
                    currSimList$kCube[simList$reshapeType1[type1Vals == i]] = params$type1Funs[i,3];
                    currSimList$nCube[simList$reshapeType1[type1Vals == i]] = params$type1Funs[i,2];
                    currSimList$gCube[simList$reshapeType1[type1Vals == i]] = params$type1Funs[i,1];
                    #Set Type 2
                    currSimList$kCube[simList$reshapeType2[type2Vals == i]] = params$type2Funs[i,3];
                    currSimList$nCube[simList$reshapeType2[type2Vals == i]] = params$type2Funs[i,2];
                    currSimList$gCube[simList$reshapeType2[type2Vals == i]] = params$type2Funs[i,1];
                }
            saveScores = array(0,dim(params$type1Funs)[1])
            for (y in 1:dim(params$type1Funs)[1]) {
                currTestModel = currCutModel
                currTestSimList = currSimList
                # overwrite parameters that you are thinking about replacing
                # to see if AND gate can be replaced with OR gate with some
                # other transfer function.
                if (currTestSimList$typeCube[OrGates[inputsTest[o]==OrIns & currOutput==OrOuts],1] == 1) {
                    currTestSimList$kCube[OrGates[inputsTest[o]==OrIns & currOutput==OrOuts],1] = params$type1Funs[y,3]
                    currTestSimList$nCube[OrGates[inputsTest[o]==OrIns & currOutput==OrOuts],1] = params$type1Funs[y,2]
                    currTestSimList$gCube[OrGates[inputsTest[o]==OrIns & currOutput==OrOuts],1] = params$type1Funs[y,1]
                }
                else if (currTestSimList$typeCube[OrGates[inputsTest[o]==OrIns & currOutput==OrOuts],1] == 2) {
                    currTestSimList$kCube[OrGates[inputsTest[o]==OrIns & currOutput==OrOuts],1] = params$type2Funs[y,3]
                    currTestSimList$nCube[OrGates[inputsTest[o]==OrIns & currOutput==OrOuts],1] = params$type2Funs[y,2]
                    currTestSimList$gCube[OrGates[inputsTest[o]==OrIns & currOutput==OrOuts],1] = params$type2Funs[y,1]
                }
                else {
                    print('Something went wrong with typing')
                }
                Param1Ix = Param1Const[OrGates[inputsTest[o]==OrIns & currOutput==OrOuts],1]
                Param2Ix = Param2Const[OrGates[inputsTest[o]==OrIns & currOutput==OrOuts],1]
                # no    w cut sim list
                currTestModel$interMat = currTestModel$interMat[,currBits]
                currTestModel$notMat<-currTestModel$notMat[,currBits]
                currTestModel$reacID<-currTestModel$reacID[currBits]
                currTestSimList$kCube = simList$kCube[currBits,]
                currTestSimList$nCube = simList$nCube[currBits,]
                currTestSimList$gCube = simList$gCube[currBits,]
                currTestSimList$ixNeg = simList$ixNeg[currBits,]
                currTestSimList$finalCube= simList$finalCube[currBits,]
                currTestSimList$ignoreCube= simList$ignoreCube[currBits,];
                currTestSimList$maxIx= simList$maxIx[currBits];
                # simulate and score

                SimResults = simFuzzyT1(CNOlist,currTestModel,currTestSimList)
                testMSE = getFit(SimResults,CNOlist,currTestModel,indexList, timePoint="t1",sizeFac=0,NAFac=1,nInTot=1)
                testMSE = testMSE/nDP
                saveScores[y] = testMSE
            }
            currScore = min(saveScores)
            bestAlt = which(saveScores==min(saveScores))
            # if one was good enough, figure out which is best and sub that in
            if (currScore <= BaseMSE + firstCutOff) {
                NewFamilyBits = rbind(NewFamilyBits,BitString)
                NewFamilyBits[dim(NewFamilyBits)[1],AndGates[t]] = FALSE
                NewFamilyBits[dim(NewFamilyBits)[1],OrGates[inputsTest[o]==OrIns & currOutput==OrOuts]] = TRUE
                NewFamilyParams= rbind(NewFamilyParams,res$bString)
                if (!is.nan(Param1Ix)) {
                    NewFamilyParams[dim(NewFamilyParams)[1],Param1Ix] = bestAlt
                    ParamSave = Param1Ix
                    }
                else if (!is.nan(Param2Ix)) {
                    NewFamilyParams[dim(NewFamilyParams)[1],Param2Ix] = bestAlt
                    ParamSave = Param2Ix
                    }
                else {
                    print('something went wrong')
                }
            ReplacementVec = rbind(ReplacementVec,cbind(AndGates[t], OrGates[inputsTest[o]==OrIns & currOutput==OrOuts], ParamSave, bestAlt, currScore))
            }
        }
    }
}

#THIS PART computes the MSE after replacement
    if (dim(ReplacementVec)[1]>0) {
        toReplace = unique(ReplacementVec[,1]);
        ReplaceList = matrix(NaN,nrow=length(toReplace),ncol=5)
        for (p in 1:length(toReplace)) {
            currentLines = which(ReplacementVec[,1]==toReplace[p])
            minIx = ReplacementVec[currentLines,5] == min(ReplacementVec[currentLines,5])
            ReplaceList[p,] = ReplacementVec[currentLines[minIx],]
        }
        testParams = res$bString;
        testParams[ReplaceList[,3]] = ReplaceList[,4]
        testBits = BitString
        testBits[ReplaceList[,1]] = FALSE
        testBits[ReplaceList[,2]] = TRUE
        interpReplace = interpretDiscreteGA(model=model,paramsList=params,intString=testParams,bitString=testBits)
        SimResults = simFuzzyT1(CNOlist,interpReplace$model,interpReplace$simList)
        testMSE = getFit(SimResults,CNOlist,interpReplace$model,indexList, timePoint="t1",sizeFac=0,NAFac=1,nInTot=1)
        testMSE = testMSE/nDP
        ActuallyReplace = ReplaceList
        finalScores[3] = testMSE;
    }
    else{
        ActuallyReplace = array(NA,0)
        finalScores[3] = finalScores[2]
    }

    finalBits = BitString
    finalParams = res$bString
    if (length(ActuallyReplace) > 0) {
        finalParams[ActuallyReplace[,3]] = ActuallyReplace[,4];
        finalBits[ActuallyReplace[,1]] = FALSE
        finalBits[ActuallyReplace[,2]] = TRUE
    }
    if (length(Reacs2Remove) > 0) {
        finalBits[Reacs2Remove] = FALSE;
    }
    interpAll = interpretDiscreteGA(model=model, 
        paramsList=params, intString=finalParams, bitString=finalBits)
    SimResults = simFuzzyT1(CNOlist,interpAll$model,interpAll$simList)
    testMSE = getFit(SimResults,CNOlist,interpAll$model,indexList, timePoint="t1",sizeFac=0,NAFac=1,nInTot=1)
    testMSE = testMSE/nDP
    finalScores[4] = testMSE

    return(list(redModel=interpAll$model, redSimList=interpAll$simList, bitString=finalBits, intString=finalParams, MSE=finalScores[4]))

}




