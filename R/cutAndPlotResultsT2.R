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

cutAndPlotResultsT2 <-function(model, bStringT1, bStringT2, simList, CNOlist,
    indexList, plotPDF=FALSE, tag=NULL, tPt=CNOlist$timeSignals[2:3])
{

    # simulate T1
    # prepare the model (i.e. cut)

    modelcut <- model
    modelcut$interMat <- modelcut$interMat[,as.logical(bStringT1)]
    modelcut$notMat <- modelcut$notMat[,as.logical(bStringT1)]
    modelcut$reacID <- modelcut$reacID[as.logical(bStringT1)]

    simListCut <- cutSimList(simList, bStringT1)

    # simulate
    SimT1 <- simulatorT1(
        CNOlist=CNOlist,
        model=modelcut,
        simList=simListCut,
        indexList=indexList)
    simResT1 <- as.matrix(SimT1[,indexList$signals])

    # simulate T2

    # prepare the model
    bitString2 <- bStringT1
    bitString2[which(bStringT1 == 0)] <- bStringT2
    BStimes <- bStringT1
    BStimes[which(bStringT1 == 0)] <- bStringT2*2

    modelcut <- model
    modelcut$interMat <- modelcut$interMat[,as.logical(bitString2)]
    modelcut$notMat <- modelcut$notMat[,as.logical(bitString2)]
    modelcut$reacID <- modelcut$reacID[as.logical(bitString2)]
    modelcut$times <- BStimes[which(BStimes != 0)]

    simListCut <- cutSimList(simList, bitString2)

    # simulate
    SimT2 <- simulatorT2(
        simResultsT1=SimT1,
        CNOlist=CNOlist,
        model=modelcut,
        simList=simListCut,
        indexList=indexList)
    simResT2 <- SimT2[,indexList$signals]

    # put it all together

    Sim0 <- simulatorT0(CNOlist=CNOlist,model=modelcut,simList=simListCut,indexList=indexList)
    simResT0 <- as.matrix(Sim0[,indexList$signals])

    simResults <- list(
        t0=simResT0,
        t1=simResT1,
        t2=simResT2)


    plotOptimResultsPan(
        simResults=simResults,
        CNOlist=CNOlist,
        formalism="ss2",
        tPt=tPt
    )

    if(plotPDF == TRUE) {
        if(is.null(tag)) {
            filename <- paste(deparse(substitute(model)),"SimResultsT1T2.pdf",sep="")
        }
        else {
            filename<-paste(tag, "SimResultsT1T2.pdf", sep="_")
        }

        plotOptimResultsPan(
            simResults=simResults,
            CNOlist=CNOlist,
            formalism="ss2",
            tPt=tPt,
            pdfFileName=filename,
            pdf=TRUE
        )
    }
}

