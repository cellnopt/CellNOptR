upSignalGraph <-
function(signal,graph,stopNodes){
  
  gg<-graph
  
  nod<-nodes(gg)
  edg<-edges(gg)
  
  SigUp<-signal
  
  findUp<-function(x){
    tmp<-vector()
    for (i in 1:length(nod)){
      if (any(edg[[i]]==x)){
         tmp<-c(tmp,nod[i])
      }
    }
    return(tmp)
  }
  
  tmp<-findUp(signal)
  tmp<-setdiff(tmp,stopNodes)
  queue<-tmp
  
  while (length(queue)>0){
    tmp<-findUp(queue[1])
    tmp<-setdiff(tmp,stopNodes)
    queue<-c(queue, tmp)
    SigUp<-c(SigUp,queue[1])
    queue<-setdiff(queue,SigUp)
    #queue<-queue[-1]
  }

  if (length(grep("and",SigUp))>0){
    SigUp<-SigUp[-grep("and",SigUp)]
  }
  return(SigUp)
}