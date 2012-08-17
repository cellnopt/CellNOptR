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

CNORbool<-function(CNOlist, model, paramsList=defaultParameters(),
    compression=TRUE, expansion=TRUE, cutNONC=TRUE, verbose=FALSE)
{

    cnolist = CNOlist

    if (is.character(cnolist)==TRUE){
        cnolist = makeCNOlist(readMIDAS(cnolist), subfield=FALSE)
    }
    if (is.character(model)==TRUE){
        model = readSIF(model)
    }

    paramsList$verbose = verbose
    Name = "CNORbool"
    #1.Plot the CNOlist
    plotCNOlist(cnolist)
    plotCNOlistPDF(cnolist, filename=paste(Name,"DataPlot.pdf",sep="_"))

    #2. Checks data to model compatibility
    checkSignals(cnolist, model)
    #3.Cut the nonc off the model
    #4.Compress the model
    #5.Expand the gates
    res = preprocessing(cnolist, model,
        compression=compression, expansion=expansion, cutNONC=cutNONC, verbose=verbose)

    #6.Compute the residual error
    resE<-residualError(cnolist)

    #7.Prepare for simulation

    #8.Optimisation t1
    bStrings = list()
    initBstring<-rep(1,length(res$model$reacID))
    print("Entering gaBinaryT1")
    T1opt<-gaBinaryT1(CNOlist=cnolist,
        model=res$model,
        initBstring=initBstring,
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
    bStrings[[1]] = T1opt$bString
    #9.Plot simulated and experimental results
    #cutAndPlot(cnolist, res$model, bStrings=bStrings,plotPDF=TRUE)

    #10.Plot the evolution of fit
    pdf(paste(Name,"evolFitT1.pdf",sep="_"))
    plotFit(optRes=T1opt)
    dev.off()
    plotFit(optRes=T1opt)

    #11.Optimise t2
    Times = 1
    T2opt<-NA # default value

    for (i in 3:length(cnolist$valueSignals)){
        Times = Times + 1
        print(paste("Entering gaBinaryTN (time=", Times, ")", sep=" "))

        TNopt<-gaBinaryTN(
            CNOlist=cnolist,
            model=res$model,
            bStrings=bStrings,
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

        bStrings[[Times]] = TNopt$bString
        print(TNopt$bString)

        pdf(paste(Name,"evolFitT", i, ".pdf",sep=""))
        plotFit(optRes=TNopt)
        dev.off()
        plotFit(optRes=TNopt)
        if (Times==2){
            T2opt = TNopt
        }
    }

    print("ga done")
    #13.Write the scaffold and PKN
    #and
    #14.Write the report
    writeScaffold(
        modelComprExpanded=res$model,
        optimResT1=T1opt,
        optimResT2=T2opt,
        modelOriginal=model,
        CNOlist=cnolist)
    writeNetwork(
        modelOriginal=model,
        modelComprExpanded=res$model,
        optimResT1=T1opt,
        optimResT2=T2opt,
        CNOlist=cnolist)

    if(Times==2){

        namesfiles<-list(
            dataPlot=paste(Name,"DataPlot.pdf",sep=""),
            evolFitT1=paste(Name,"evolFitT1.pdf",sep=""),
            evolFitT2=paste(Name,"evolFitT2.pdf",sep=""),
            simResults2="NCNOcutCompExpSimResultsT1T2.pdf",
            simResults1="NCNOcutCompExpSimResultsT1.pdf",
            scaffold="Scaffold.sif",
            scaffoldDot="Scaffold.dot",
            tscaffold="TimesScaffold.EA",
            wscaffold="weightsScaffold.EA",
            PKN="PKN.sif",
            PKNdot="PKN.dot",
            wPKN="TimesPKN.EA",
            nPKN="nodesPKN.NA")

    }
    else{

        namesfiles<-list(
            dataPlot=paste(Name,"DataPlot.pdf",sep=""),
            evolFitT1=paste(Name,"evolFitT1.pdf",sep=""),
            evolFitT2=NA,
            simResults2=NA,
            simResults1="NCNOcutCompExpSimResultsT1.pdf",
            scaffold="Scaffold.sif",
            scaffoldDot="Scaffold.dot",
            tscaffold="TimesScaffold.EA",
            wscaffold="weightsScaffold.EA",
            PKN="PKN.sif",
            PKNdot="PKN.dot",
            wPKN="TimesPKN.EA",
            nPKN="nodesPKN.NA")
    }

    #writeReport(
    #    modelOriginal=paramsList$model,
    #    modelOpt=NCNOcutCompExp,
    #    optimResT1=T1opt,
    #    optimResT2=T2opt,
    ##    CNOlist=paramsList$data,
     #   directory=Name,
     #   namesFiles=namesfiles,
     #   namesData=NamesData,
     #   resE=resE)
    return(list(model=res$model, bStrings=bStrings))
}
