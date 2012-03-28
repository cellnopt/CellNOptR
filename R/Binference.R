Binference <-
function(CNOlist, Model="AIC", tempCheckOrders=10, maxIter=100, filename="BAYESIAN"){
  
  library(catnet)
  
  ##STIMULI
  #data matrix
  valueStimuli<-CNOlist$valueStimuli+1
  colnames(valueStimuli)<-CNOlist$namesStimuli

  #parturbations matrix
  fixedStimuli<-valueStimuli

  ##INHIBITORS
  namesInhibitors<-setdiff(CNOlist$namesInhibitors, CNOlist$namesSignals)
  valueInhibitors<-CNOlist$valueInhibitors
  colnames(valueInhibitors)<-CNOlist$namesInhibitors
  valueInhibitors<-as.matrix(valueInhibitors[,namesInhibitors])
  colnames(valueInhibitors)<-namesInhibitors
  # = 0 not inhibited (unknown value) - remains = 0 in the perturbations matrix (will be = NA in the data matrix)
  # = 1 inhibited - remains = 1 (off) in both the perturbations matrix and in the data matrix

  #parturbations matrix
  fixedInhibitors<-valueInhibitors
  colnames(fixedInhibitors)<-NULL

  #data matrix
  valueInhibitors[valueInhibitors==0]=2
  #valueInhibitors[valueInhibitors==0]=NA

  ##SIGNALS
  valueSignals<-CNOlist$valueSignals[[2]]

  #parturbations matrix
  fixedSignals<-matrix(data=0,nrow=dim(valueSignals)[1],ncol=dim(valueSignals)[2])
  fixedSignals[is.na(valueSignals)]<-1

  #data matrix
  valueSignals[is.na(valueSignals)]<-0
  valueSignals[valueSignals>0.5]<-2
  valueSignals[valueSignals<0.5]<-1
  colnames(valueSignals)<-CNOlist$namesSignals

  ## ALL
  #data matrix
  dataset<-cbind(valueStimuli, valueInhibitors, valueSignals)
  #parturbations matrix
  perturbations<-cbind(fixedStimuli, fixedInhibitors, fixedSignals)
  #perturbations[perturbations>0]<-1
  
  
  psamples<-t(dataset)
  perturbations<-t(perturbations)

  parentSizes<-rep(3,dim(psamples)[1])
  parentSizes[1:length(CNOlist$namesStimuli)]<-0

  nets <- cnSearchSA(data=psamples, perturbations=perturbations, maxParentSet=3,
                       parentSizes=parentSizes, maxComplexity=0,
                       parentsPool=NULL, fixedParents=NULL,
                       edgeProb=NULL, dirProb=NULL, selectMode = "AIC",
                       tempStart=1, tempCoolFact=0.1, tempCheckOrders=tempCheckOrders,
                       maxIter=maxIter, 
                       orderShuffles=1, stopDiff=0,
                       numThreads=2, priorSearch=NULL, echo=FALSE)


  AllLinks<-matrix(nrow=0,ncol=2)
  for (i in 1:length(nets@nets)){
    netTMP<-nets@nets[[i]]
    AllLinks<-rbind(AllLinks, cnMatEdges(netTMP))
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
