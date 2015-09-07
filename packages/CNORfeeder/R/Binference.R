#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id$
Binference <-
function(CNOlist, mode="AIC", tempCheckOrders=10, maxIter=100, filename="BAYESIAN"){
  
  # library(catnet)
  requireNamespace("catnet")
	
  if ((class(CNOlist)=="CNOlist")==FALSE){
	CNOlist = CellNOptR::CNOlist(CNOlist)
  }
  
  ##STIMULI
  #data matrix
  valueStimuli<-CNOlist@stimuli
	
  #parturbations matrix
  fixedStimuli<-valueStimuli

  ##INHIBITORS
  namesInhibitors <- setdiff(colnames(CNOlist@inhibitors), colnames(CNOlist@signals[[1]]))
  valueInhibitors <- CNOlist@inhibitors[,namesInhibitors]
  # = 0 not inhibited (unknown value) - remains = 0 in the perturbations matrix (will be = NA in the data matrix)
  # = 1 inhibited - remains = 1 (off) in both the perturbations matrix and in the data matrix

  #parturbations matrix
  fixedInhibitors<-valueInhibitors
  colnames(fixedInhibitors)<-NULL

  #data matrix
  valueInhibitors[valueInhibitors==0]=2
  #valueInhibitors[valueInhibitors==0]=NA

  ##SIGNALS
  valueSignals<-CNOlist@signals[[2]]

  #parturbations matrix
  fixedSignals<-matrix(data=0,nrow=dim(valueSignals)[1],ncol=dim(valueSignals)[2])
  fixedSignals[is.na(valueSignals)]<-1

  #data matrix
  valueSignals[is.na(valueSignals)]<-0
  valueSignals[valueSignals>0.5]<-2
  valueSignals[valueSignals<0.5]<-1

  ## ALL
  #data matrix
  dataset<-cbind(valueStimuli, valueInhibitors, valueSignals)
  #parturbations matrix
  perturbations<-cbind(fixedStimuli, fixedInhibitors, fixedSignals)
  #perturbations[perturbations>0]<-1
  
  
  psamples<-t(dataset)
  perturbations<-t(perturbations)

  parentSizes<-rep(3,dim(psamples)[1])
  parentSizes[1:length(colnames(CNOlist@stimuli))]<-0

  nets <- catnet::cnSearchSA(data=psamples, perturbations=perturbations, maxParentSet=3,
                       parentSizes=parentSizes, maxComplexity=0,
                       parentsPool=NULL, fixedParents=NULL,
                       edgeProb=NULL, dirProb=NULL, selectMode = mode,
                       tempStart=1, tempCoolFact=0.1, tempCheckOrders=tempCheckOrders,
                       maxIter=maxIter, 
                       orderShuffles=1, stopDiff=0,
                       numThreads=2, priorSearch=NULL, echo=FALSE)


  AllLinks<-matrix(nrow=0,ncol=2)
  for (i in 1:length(nets@nets)){
    netTMP<-nets@nets[[i]]
    AllLinks<-rbind(AllLinks, catnet::cnMatEdges(netTMP))
  }
  AllLinksUn<-unique(AllLinks)

  
  freqVec<-vector()
  for (j in 1:dim(AllLinksUn)[1]){
    freqVec[j]<-length(intersect(which(AllLinks[,1]==AllLinksUn[j,1]), which(AllLinks[,2]==AllLinksUn[j,2])))
  }
  freqVec<-freqVec/length(nets@nets)
  AllLinksOK<-AllLinksUn[freqVec>0.1,]

  sif<-cbind(AllLinksOK[,1], rep(1,dim(AllLinksOK)[1]), AllLinksOK[,2])
  
  write.table(sif, file=paste(filename,".sif",sep=""), row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

  return(sif)
}
