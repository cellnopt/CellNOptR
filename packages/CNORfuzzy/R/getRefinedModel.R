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
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/software.html
#
##############################################################################
getRefinedModel = function(res, CNOlist, cutModel, cutSimList, refParams) {
    # ResRef = CNOGetRefinedModelInside(Res,CNOProject,varargin)
    # Refines cFL model parameters "inside" CNOGetModel.
    # Contributed by Melody K Morris Jan 2011
    # TODO(tc): lower and upper bound are harcoded. In particular upper is 68.5098
    # should be 68.9058 from a quick solver test makde on the hill equation.
    if ((class(CNOlist)=="CNOlist")==FALSE){
        CNOlist = CellNOptR::CNOlist(CNOlist)
    } 


    indexList<-indexFinder(CNOlist=CNOlist,model=cutModel, verbose=FALSE)


    verbose = refParams$verbose
    numTF1 = dim(refParams$type1Funs)[1]
    numTF2 = dim(refParams$type2Funs)[1]
    if (numTF1==numTF2) {
        numTF = numTF1
    }
    else {
        print('Number of type 1 and type 2 transfer functions must be equal!  Minimum number used')
        numTF =  min(numTF1,numTF2);
    }

    #This loop determined what will be refined after the discrete genetic algorithm.
    if (verbose){
        print('Only Relevent parameters will be refined')
    }

    if (all(all(refParams$type1Funs[,2:3]==1)) || all(refParams$type1Funs[,2]==1.01)) {
        if (verbose){
            print('Type I functions are linear')
        }
        refParams$type1Case='Linear'
        refParams$type1NumParams=1
    }
    else if (all(refParams$type1Funs[,1] == 1)) {
        if (verbose){
            print('Type I functions are normalized Hill (these can also be non-normalized by setting refParams$tfun)')
        }
        refParams$type1Case='NormHill'
        refParams$type1NumParams=2
        }
    else {
        if (verbose){
            print('Type I functions are normalized hill with gain')
        }
        refParams$type1Case='Hill'
        refParams$type1NumParams=3
    }
    if (all(all(refParams$type2Funs[,2:3]==1)) || all(refParams$type2Funs[,2]==1.01)) {
        if (verbose){
            print('Type II functions are linear')
        }
        refParams$type2Case='Linear'
        refParams$type2NumParams=1
        }
    else if (all(refParams$type2Funs[,1] == 1)) {
        if (verbose){
            print('Type II functions are normalized Hill (these can also be non-normalized by setting refParams$tfun)')
        }
        refParams$type2Case='NormHill'
        refParams$type2NumParams=2
        }
    else {
        if (verbose){
            print('Type II functions are normalized hill with gain')
        }
        refParams$type2Case='Hill'
        refParams$type2NumParams=3
    }



# objective fucntion

objFunParams <- function(SetString) {
# sub in parameters depending on the type of transfer functions used
    if (refParams$type1Case == 'Linear') {
    type1ValsG = SetString[1:cutSimList$numType1]
    cutSimList$gCube[cutSimList$reshapeType1] = type1ValsG
    }
    else if (refParams$type1Case=='NormHill') {
    type1ValsK = SetString[1:cutSimList$numType1]
    type1ValsN = SetString[(cutSimList$numType1+1):(2*cutSimList$numType1)]
    cutSimList$kCube[cutSimList$reshapeType1] = type1ValsK
    cutSimList$nCube[cutSimList$reshapeType1] = type1ValsN
    }
    else {
    type1ValsK = SetString[1:cutSimList$numType1]
    type1ValsN = SetString[(cutSimList$numType1+1):(2*cutSimList$numType1)]
    type1ValsG = SetString[(2*cutSimList$numType1+1):(3*cutSimList$numType1)]
    cutSimList$kCube[cutSimList$reshapeType1] = type1ValsK
    cutSimList$nCube[cutSimList$reshapeType1] = type1ValsN
    cutSimList$gCube[cutSimList$reshapeType1] = type1ValsG
    }

    if (refParams$type2Case == 'Linear') {
    type2ValsG = SetString[(refParams$type1NumParams*cutSimList$numType1+1):(length(SetString))];
    cutSimList$gCube[cutSimList$reshapeType2] = type2ValsG
    }
    else if (refParams$type2Case=='NormHill') {
    type2ValsK = SetString[(refParams$type1NumParams*cutSimList$numType1+1):(refParams$type1NumParams*cutSimList$numType1+cutSimList$numType2)]
    type2ValsN = SetString[(refParams$type1NumParams*cutSimList$numType1+cutSimList$numType2+1):(length(SetString))]
    cutSimList$kCube[cutSimList$reshapeType2] = type2ValsK
    cutSimList$nCube[cutSimList$reshapeType2] = type2ValsN
    }
    else {
    type2ValsK = SetString[(refParams$type1NumParams*cutSimList$numType1+1):(refParams$type1NumParams*cutSimList$numType1+cutSimList$numType2)]
    type2ValsN = SetString[(refParams$type1NumParams*cutSimList$numType1+cutSimList$numType2+1):(refParams$type1NumParams*cutSimList$numType1+2*cutSimList$numType2)]
    type2ValsG = SetString[(refParams$type1NumParams*cutSimList$numType1+2*cutSimList$numType2+1):(length(SetString))]
    cutSimList$kCube[cutSimList$reshapeType2] = type2ValsK
    cutSimList$nCube[cutSimList$reshapeType2] = type2ValsN
    cutSimList$gCube[cutSimList$reshapeType2] = type2ValsG
    }
# simulate
    SimResults = simFuzzyT1(CNOlist=CNOlist, model=cutModel, simList=cutSimList )
# calculate score
    Score<-getFit(simResults=SimResults,CNOlist=CNOlist,model=cutModel,indexList=indexList,timePoint="t1",sizeFac=0,NAFac=refParams$NAFac,nInTot=1)
    nDataP<-sum(!is.na(CNOlist@signals[[2]]))
    Score<-Score/nDataP
    return(Score)
}

 #set lower and upper bounds
if (refParams$type1Case=='Linear') {
    SetStringType1 = cbind(t(cutSimList$gCube[cutSimList$reshapeType1]))
    lowBoundType1 = cbind(t(array(0,cutSimList$numType1)))
    upBoundType1 = cbind(t(array(1,cutSimList$numType1)))
    }
else if (refParams$type1Case=='NormHill') {
    SetStringType1 = cbind(t(cutSimList$kCube[cutSimList$reshapeType1]),t(cutSimList$nCube[cutSimList$reshapeType1]))
    lowBoundType1 = cbind(t(array(0.05,cutSimList$numType1)),t(array(1.01,cutSimList$numType1)))
    upBoundType1 = cbind(t(array(68.5098,cutSimList$numType1)),t(array(20,cutSimList$numType1)))
    }
else {
    SetStringType1 = cbind(t(cutSimList$kCube[cutSimList$reshapeType1]),t(cutSimList$nCube[cutSimList$reshapeType1]),t(cutSimList$gCube[cutSimList$reshapeType1]));
    lowBoundType1 = cbind(t(array(0.05,cutSimList$numType1)),t(array(1.01,cutSimList$numType1)),t(array(0,cutSimList$numType1)));
    upBoundType1 = cbind(t(array(68.5098,cutSimList$numType1)),t(array(20,cutSimList$numType1)),t(array(1,cutSimList$numType1)));
}

if (refParams$type2Case=='Linear') {
    SetStringType2 = cbind(t(cutSimList$gCube[cutSimList$reshapeType2]))
    lowBoundType2 =  cbind(t(array(0,cutSimList$numType2)))
    upBoundType2 = cbind(t(array(1,cutSimList$numType2)))
    }
else if (refParams$type2Case=='NormHill') {
    SetStringType2 = cbind(t(cutSimList$kCube[cutSimList$reshapeType2]),t(cutSimList$nCube[cutSimList$reshapeType2]))
    lowBoundType2 = cbind(t(array(0.05,cutSimList$numType2)),t(array(1.01,cutSimList$numType2)))
    upBoundType2 = cbind(t(array(68.5098,cutSimList$numType2)),t(array(20,cutSimList$numType2)))
    }
else {
    SetStringType2 = cbind(t(cutSimList$kCube[cutSimList$reshapeType2]),t(cutSimList$nCube[cutSimList$reshapeType2]),t(cutSimList$gCube[cutSimList$reshapeType2]))
    lowBoundType2 = cbind(t(array(0.05,cutSimList$numType2)),t(array(1.01,cutSimList$numType2)),t(array(0,cutSimList$numType2)))
    upBoundType2 = cbind(t(array(68.5098,cutSimList$numType2)),t(array(20,cutSimList$numType2)),t(array(1,cutSimList$numType2)))
}

SetString0 = cbind(SetStringType1,SetStringType2);
lowBound = cbind(lowBoundType1,lowBoundType2);
upBound = cbind(upBoundType1,upBoundType2);

# Get/Set the default or user parameters related to the optimisation.
if (is.null(refParams$optimisation$algorithm)==TRUE) {
    refParams$optimisation$algorithm='NLOPT_LN_SBPLX'
    }
if (is.null(refParams$optimisation$xtol_abs)==TRUE) {
    refParams$optimisation$xtol_abs=0.001
}
if (is.null(refParams$optimisation$maxeval)==TRUE) {
    refParams$optimisation$maxeval=10000
}
if (is.null(refParams$optimisation$maxtime)==TRUE) {
    refParams$optimisation$maxtime = 5*60
}


ConOptions=list(
    algorithm=refParams$optimisation$algorithm,
    xtol_abs=refParams$optimisation$xtol_abs,
    maxeval=refParams$optimisation$maxeval,
    maxtime=refParams$optimisation$maxtime)


LocalMinList = nloptr(SetString0,objFunParams,lb=lowBound,ub=upBound,opts=ConOptions);
finalSet = LocalMinList$solution
MSE = LocalMinList$objective

# now save everything
ResRef = list();
ResRef$refinedModel = cutModel
ResRef$refinedSimList = cutSimList

# plug in optimal paramters one last time.
if (refParams$type1Case=='Linear') {
    type1ValsG = finalSet[1:ResRef$refinedSimList$numType1]
    ResRef$refinedSimList$gCube[ResRef$refinedSimList$reshapeType1] = type1ValsG
    }
else if (refParams$type1Case=='NormHill') {
    type1ValsK = finalSet[1:ResRef$refinedSimList$numType1]
    type1ValsN = finalSet[(ResRef$refinedSimList$numType1+1):(2*ResRef$refinedSimList$numType1)]
    ResRef$refinedSimList$kCube[ResRef$refinedSimList$reshapeType1] = type1ValsK
    ResRef$refinedSimList$nCube[ResRef$refinedSimList$reshapeType1] = type1ValsN
    }
else {
    type1ValsK = finalSet[1:ResRef$refinedSimList$numType1]
    type1ValsN = finalSet[(ResRef$refinedSimList$numType1+1):(2*ResRef$refinedSimList$numType1)]
    type1ValsG = finalSet[(2*ResRef$refinedSimList$numType1+1):(3*ResRef$refinedSimList$numType1)]
    ResRef$refinedSimList$kCube[ResRef$refinedSimList$reshapeType1] = type1ValsK
    ResRef$refinedSimList$nCube[ResRef$refinedSimList$reshapeType1] = type1ValsN
    ResRef$refinedSimList$gCube[ResRef$refinedSimList$reshapeType1] = type1ValsG
}

if (refParams$type2Case == 'Linear') {
    type2ValsG = finalSet[(refParams$type1NumParams*ResRef$refinedSimList$numType1+1):length(finalSet)]
    ResRef$refinedSimList$gCube[ResRef$refinedSimList$reshapeType2] = type2ValsG
    }
else if (refParams$type2Case == 'NormHill') {
    type2ValsK = finalSet[(refParams$type1NumParams*ResRef$refinedSimList$numType1+1):(refParams$type1NumParams*ResRef$refinedSimList$numType1+ResRef$refinedSimList$numType2)]
    type2ValsN = finalSet[(refParams$type1NumParams*ResRef$refinedSimList$numType1+ResRef$refinedSimList$numType2+1):(length(finalSet))]
    ResRef$refinedSimList$kCube[ResRef$refinedSimList$reshapeType2] = type2ValsK
    ResRef$refinedSimList$nCube[ResRef$refinedSimList$reshapeType2] = type2ValsN
    }
else {
    type2ValsK = finalSet[(refParams$type1NumParams*ResRef$refinedSimList$numType1+1):(refParams$type1NumParams*ResRef$refinedSimList$numType1+ResRef$refinedSimList$numType2)]
    type2ValsN = finalSet[(refParams$type1NumParams*ResRef$refinedSimList$numType1+ResRef$refinedSimList$numType2+1):(refParams$type1NumParams*ResRef$refinedSimList$numType1+2*ResRef$refinedSimList$numType2)]
    type2ValsG = finalSet[(refParams$type1NumParams*ResRef$refinedSimList$numType1+2*ResRef$refinedSimList$numType2+1):(length(finalSet))]
    ResRef$refinedSimList$kCube[ResRef$refinedSimList$reshapeType2] = type2ValsK
    ResRef$refinedSimList$nCube[ResRef$refinedSimList$reshapeType2] = type2ValsN
    ResRef$refinedSimList$gCube[ResRef$refinedSimList$reshapeType2] = type2ValsG
}
return(list( refModel = ResRef,finalSet = finalSet, MSE = MSE))

}
