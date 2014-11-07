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
# $Id: test_computeScore.R 2803 2012-11-21 10:59:05Z cokelaer $
# This is a unit test that checks te output of the computeScore function


test_computeScore <- function(){

    library(CNORfuzzy)
    library(RUnit)

    # load some data
    data(CNOlistToy, package="CellNOptR")
    data(ToyModel, package="CellNOptR")

    data = CNOlistToy
    pknmodel = ToyModel


    # prepares the parameters
    paramsList<-defaultParametersFuzzy(data, pknmodel)
    paramsList$verbose = FALSE

    # performs some preprocessing
    model = preprocessing(data, pknmodel)
    indices = indexFinder(data, model)
    simList = prep4simFuzzy(model=model,  paramsList=paramsList)

#initBstring<-(sample.int(dim(paramsList$type2Funs)[1],(fields4Sim$numType1+fields4Sim$numType2),replace = TRUE)) - 1

#T1opt<-gaDiscreteT1(CNOlist=paramsList$data,
#        Model=preModel,
#        SimList=fields4Sim,
#        indexList=indices,
#        initBstring=initBstring,
#        paramsList = paramsList,
#        sizeFac=paramsList$sizeFac,
#        NAFac=paramsList$NAFac,
#        PopSize=paramsList$PopSize,
#        Pmutation=paramsList$Pmutation,
#        MaxTime=paramsList$MaxTime,
#        maxGens=paramsList$maxGens,
#        StallGenMax=paramsList$StallGenMax,
#        SelPress=paramsList$SelPress,
#        elitism=paramsList$elitism,
#        RelTol=paramsList$RelTol,
#        verbose=paramsList$verbose)

    # test a specific parameter list
    initBstring = c(5, 1, 5, 4, 3, 4, 5, 5, 1, 3, 5, 4, 5, 2, 6, 6, 2, 2, 6)
    t1 = Sys.time()
    # run 100 times (not needed but allow sto have an idea of the time spent in
    # computeScoreFuzzy.
    for (i in seq(1,100)){
        #res = computeScoreFuzzy(paramsList$data, model, simList,  indices, paramsList, initBstring)
        res = computeScoreFuzzy(paramsList$data, model, NULL, NULL, paramsList, initBstring)
    }
    t2 = Sys.time()
    t2-t1
    print(res)
    
    #checkEquals(res, 0.131629217146272)
    checkEquals(res, 0.1316292, tolerance=1e-6)
    
    
    # 30 seconds on linux 64 bits august 2012 (TC.) for N = 1000
}
