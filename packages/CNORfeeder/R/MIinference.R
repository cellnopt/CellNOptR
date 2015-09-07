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
MIinference <-
function(CNOlist, method="ARACNE", PKNgraph=NULL, filename="ARACNE"){
  
  # library(minet)
  requireNamespace("minet")
	
  if ((class(CNOlist)=="CNOlist")==FALSE){
	  CNOlist = CellNOptR::CNOlist(CNOlist)
  }
	
  valueStimuli <- CNOlist@stimuli
  
  namesInhibitors <- setdiff(colnames(CNOlist@inhibitors), colnames(CNOlist@signals[[1]]))
  valueInhibitors <- CNOlist@inhibitors[,namesInhibitors]
  
  valueInhibitors<-1-valueInhibitors
  #valueInhibitors[valueInhibitors==1]<-NA
  
  valueSignals<-CNOlist@signals[[2]]
  valueSignals[is.na(valueSignals)]<-0
  
  dataset<-cbind(valueStimuli, valueInhibitors, valueSignals)

  mim<-minet::build.mim(dataset=dataset)

  if (method=="ARACNE"){
    net<-minet::aracne(mim, eps=0)
  }else if (method=="CLR"){
    net<-minet::clr(mim)
  }
  
  graph<-as(net, "graphNEL")
  sif<-graph2sif(graph, writeSif=FALSE)
  
  # delete all links going to stimuli
  sif<-sif[-which(sif[,3]%in%colnames(CNOlist@stimuli)),]
  
  # if a PKN is given as input, the directionality of the links is selected according to the PKN
  
  if (!is.null(PKNgraph)){
    ixRem<-vector()
    for (i in 1:dim(sif)[1]){
      node1<-sif[i,1]
      node2<-sif[i,3]
      ck<-searchLinkGraph(node1,node2,PKNgraph)
      if (ck==1){
        oppLink<-intersect(which(sif[,1]==node2), which(sif[,3]==node1))
        if (length(oppLink)>0){
          #print(c(node1, node2))
          ixRem<-c(ixRem,oppLink)
        }
      }
    }
    if (length(ixRem>0)){sif<-sif[-ixRem,]}
  }
  
  
  write.table(sif, file=paste(filename,".sif",sep=""), row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

  return(sif)
}
