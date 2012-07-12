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

cutAndPlotResultsTN <-function(
  model,
  res,
  simList,
  CNOlist,
  indexList,
  plotPDF=FALSE, 
  tag=NULL,
  tPt=CNOlist$timeSignals[2:length(CNOlist$timeSignals)]){
  

  # prepare the model (i.e. cut)
  bStringT1<-res$Opt[[2]]$bString
  modelCut <- model
  modelCut$interMat <- modelCut$interMat[,as.logical(bStringT1)]
  modelCut$notMat <- modelCut$notMat[,as.logical(bStringT1)]
  modelCut$reacID <- modelCut$reacID[as.logical(bStringT1)]
  
  simListCut <- cutSimList(simList, res$Opt[[2]]$bString)
  
  # simulate
  
  Sim0 <- simulatorT0(CNOlist=CNOlist,model=modelCut,simList=simListCut,indexList=indexList)
  SimResT0 <- as.matrix(Sim0[,indexList$signals])
  simResults<-list()
  simResults[[1]]<-SimResT0
  for(i in 2:length(res$simRes)){
    simResTN<-res$simRes[[i]]
    cutRes<-simResTN[,indexList$signals]
    simResults[[i]]<-cutRes
  
  }
  
#   SimResults <- list(
#     t0=SimResT0,
#     t1=SimResT1)
#   
#   
#   # simulate T2
#   bStringPrev<-bStringT1
#   # prepare the model
#   SimResPrev<-SimT1
#   for(i in 3:length(bStringTN)){
#     SimResTN<-Res$simRes[[i]]
#     cutRes<-SimResTN[,indexList$signals]
#     SimResults[[i]]<-cutRes
# #   bitStringCurr <- bStringPrev
# #   bitStringCurr[which(bStringPrev == 0)] <- bStringTN[[i]]
# #   BStimes <- bStringT1
# #   BStimes[which(bStringPrev == 0)] <- bStringTN[[i]]*(i-1)
# #   
# #   Modelcut <- Model
#   Modelcut$interMat <- Modelcut$interMat[,as.logical(bitStringCurr)]
#   Modelcut$notMat <- Modelcut$notMat[,as.logical(bitStringCurr)]
#   Modelcut$reacID <- Modelcut$reacID[as.logical(bitStringCurr)]
#   Modelcut$times <- BStimes[which(BStimes != 0)]
#   
#   SimListCut <- cutSimList(SimList, bitStringCurr)
#   
# #   # simulate
# #   simTN<-simulateTN(CNOlist=CNOlist,
# #                           Model=Model,
# #                           bStringTN=bStringTN[[i]],
# #                           SimList=SimListCut,
# #                           indexList=indexList,
# #                           SimResultsPrev=SimResPrev)
# #   
# #   
#   SimTN <- simulatorTN(
#     SimResultsPrev=SimResPrev,
#     CNOlist=CNOlist,
#     Model=Modelcut,
#     SimList=SimListCut,
#     indexList=indexList,
#     TimePoint=(i-1))
#   SimResPrev<-SimTN
#   SimResults[[i]] <- SimTN[,indexList$signals]
#   bStringPrev<-bitStringCurr
# }

  plotOptimResultsPan(
    simResults=simResults,
    CNOlist=CNOlist,
    formalism="ssN",
    tPt=tPt,
    #timePoints=length(tPt)
  )
  
  if(plotPDF == TRUE) {
    if(is.null(tag)) {
      filename <- paste(deparse(substitute(Model)),"SimResultsT1T2.pdf",sep="")
    }
    else {
      filename<-paste(tag, "SimResultsT1T2.pdf", sep="_")
    }
    
    plotOptimResultsPan(
      simResults=simResults,
      CNOlist=CNOlist,
      formalism="ssN",
      tPt=tPt,
      pdfFileName=filename,
      pdf=TRUE,
      #TimePoints=length(tPt)
    )
  }		
}
