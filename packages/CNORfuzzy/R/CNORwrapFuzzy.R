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
# $id: CNORwrapFuzzy.R 1499 2012-06-19 14:58:06Z cokelaer $

#This function is a wrapper around the whole CNOR analysis with cFL, it performs the following steps:
#1.Plot the CNOlist
#2.Checks data to model compatibility
#3.Find the indices, in the model, of the species that are inh/stim/sign
#4.Find the indices of the non-osb/non-contr
#5.Cut the nonc off the model
#6.Recompute the indices
#7.Compress the model
#8.Recompute the indices
#9.Expand the gates
#10.Compute the residual error
#11.Prepare for simulation
#12.Optimisation t1

CNORwrapFuzzy<-function(data, model, paramsList=NULL, verbose=TRUE){
#if the paramsList is set to NA, it means that only default parameters have been provided
#so we are going to build the parameters list with those

    if ((class(data)=="CNOlist")==FALSE){
        warning("Your input data is not a CNOlist class. Trying to convert it.
Please use CNOlist (not makeCNOlist). Future version of CellNOptR and CNORfuzzy
may not support makeCNOlist output anymore. ")
        data = CellNOptR::CNOlist(data)
    }



    ticBeg=Sys.time()
    if(is.null(paramsList)==TRUE){
        # Get default parameters
        print("Warning: get default parameters")
        paramsList<-defaultParametersFuzzy(data, model)
        }
    else{
        # make sure that the paramsList Data and Model are set.
        paramsList$data = data
        paramsList$model = model
    }

    #2.Checks data to model compatibility
    checkSignals(CNOlist=paramsList$data,model=paramsList$model)
    #!!!! model (the PKN provided as an argument) is overwritten here
    model = preprocessing(paramsList$data, paramsList$model, verbose=FALSE)

    #. Optimisation t1
    if (verbose==TRUE){
        print('Begining Optimization')
    }

    # generate a random bistring of numbers between 0 and
    # dim(paramsList$type2Funs)[1] of length fields4Sim$numType1+fields4Sim$numType2
    # fields4Sim<-prep4simFuzzy(model=NONCcutCompExp,paramsList=paramsList)
    # initBstring<- (sample.int(dim(paramsList$type2Funs)[1],(fields4Sim$numType1+fields4Sim$numType2),replace=TRUE)) - 1

    t0 = Sys.time()

    T1opt<-gaDiscreteT1(CNOlist=paramsList$data,
        model=model,
        paramsList=paramsList,
        sizeFac=paramsList$sizeFac,
        NAFac=paramsList$NAFac,
        popSize=paramsList$popSize,
        pMutation=paramsList$pMutation,
        maxTime=paramsList$maxTime,
        maxGens=paramsList$maxGens,
        stallGenMax=paramsList$stallGenMax,
        selPress=paramsList$selPress,
        elitism=paramsList$elitism,
        relTol=paramsList$relTol,
        verbose=paramsList$verbose)
    t1 = Sys.time()

    if (verbose==TRUE){
        print(paste('Discrete GA Finished in: ', format(t1-t0), sep=""))
    }

    if (verbose==TRUE){
        print('Calling interpretDiscreteGA')
    }
    #13. Interpret and remove redundant logic gates
    interpModel=interpretDiscreteGA(model=model,paramsList=paramsList,intString=T1opt$bString)

    # start making variable
    Res = list()
    Res$t1opt = T1opt
    Res$intString = T1opt$bString
    Res$currBestDiscrete = T1opt$currBest
    Res$bit = interpModel$bitString
    Res$cutBit = interpModel$cutBitString
    Res$unRef = interpModel
    Res$paramsList = paramsList
    Res$processedModel=model

    #SBS add MSE field to Res$unRef -- it may be better enough to use Res$t1opt$currBest here instead
    SimResults<-simFuzzyT1(
      CNOlist=paramsList$data,
      model=Res$unRef$model,
      simList=Res$unRef$simList
      )

    Score<-getFit(
      simResults=SimResults,
      CNOlist=paramsList$data,
      model=Res$unRef$model,
      indexList=indexFinder(paramsList$data, model, verbose=FALSE),
      timePoint="t1",
      sizeFac=paramsList$sizeFac,
      NAFac=paramsList$NAFac,
      nInTot=length(which(Res$unRef$model$interMat==-1)))
    nDataP<-sum(!is.na(paramsList$data@signals[[2]]))
    Score<-Score/nDataP


    #Res$unRef$simResults<-SimResults #probably don't need this, but it's nice for debugging purposes to have
    Res$unRef$MSE <-Score
    # 14. Refine models returned from GA
    # only do refinement if wanted (because it is VERY slow for large problems)
    if (!is.null(paramsList$doRefinement) && paramsList$doRefinement){
        # 14.a refine the cut model
        if (verbose==TRUE){
            print('Calling first Refinement')
        }
        refModel1 = getRefinedModel(res=T1opt,
                                    CNOlist=paramsList$data,
                                    cutModel=interpModel$model,
                                    cutSimList=interpModel$simList,
                                    refParams=paramsList)
        refModel1$bit = interpModel$bitString
        t2 = Sys.time()
        if (verbose==TRUE){
            print(paste('...First Refinement Complete ', format(t2-t1), sep=""))
        }

        if (verbose==TRUE){
            print('Calling second Refinement')
        }
        #14.a refine the cut model with logically redundant gates removed
        refModel2 = getRefinedModel(res=T1opt,
                                    CNOlist=paramsList$data,
                                    cutModel=interpModel$cutModel,
                                    cutSimList=interpModel$cutSimList,
                                    refParams=paramsList)
        refModel2$bit = interpModel$cutBitString
        t3 = Sys.time()
        if (verbose==TRUE){
            print(paste('...Second Refinement Complete ', format(t3-t2), sep=""))
        }

        #15. Reduce Models
        t0 = Sys.time()
        RedRef = list(refModel1,refModel2)
        prevIntString = T1opt$bString
        prevBitString = interpModel$cutBitString
        for (eachT in 1:length(paramsList$redThresh)) {
            if (verbose==TRUE){
                print(paramsList$redThresh[eachT])
            }
            t4=Sys.time()
            if (verbose==TRUE){
                print(paste('Calling reduceFuzzy ', eachT, sep=""))
            }
            redModel = reduceFuzzy(firstCutOff=paramsList$redThresh[eachT],
                                   CNOlist=paramsList$data,
                                   model=model,
                                   res=T1opt,
                                   params=paramsList)


            t5=Sys.time()
            if (verbose==TRUE){
                print(paste('...done ', format(t5-t4), sep=""))
            }
# refine parameters of reduced models
            if (all(redModel$intString==prevIntString) & all(redModel$bitString==prevBitString)){
                if (verbose==TRUE){
                    print('Reduction did not change Model.  Copying previous refinement')
                }
                currRedRef = RedRef[[eachT+1]]
            }
            else{
                if (verbose==TRUE){
                    print('Calling Refinement')
                }
                currRedRef = getRefinedModel(res=T1opt,
                                             CNOlist= paramsList$data,
                                             cutModel=redModel$redModel,
                                             cutSimList=redModel$redSimList,
                                             refParams=paramsList)
                currRedRef$bit = redModel$bitString
                t6=Sys.time()
                if (verbose==TRUE){
                    print(paste('...done ', format(t6-t5), sep=""))
                }
                prevIntString = redModel$intString
                prevBitString = redModel$bitString

            }
            RedRef[[eachT+2]] = currRedRef
    }
        Res$redRef = RedRef
        t1 = Sys.time()
        if (verbose==TRUE){
            print(paste('RedRef Finished.  Total time RedRef ', format(t1-t0)))
        }
    }
    ticEnd=Sys.time()
    if (verbose==TRUE){
        print(paste('Total Time: ', format(ticEnd-ticBeg), sep=" "))
    }
    return(Res)
}
