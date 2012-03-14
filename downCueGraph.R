downCueGraph <-
function(cue,graph,stopNodes){
  
  gg<-graph
  
  nod<-nodes(gg)
  edg<-edges(gg)
  
  CueDown<-cue
  
  if (any(which(nod==cue))){
    tmp<-edg[[which(nod==cue)]]
  }else{
    tmp<-vector()
  }
  
  tmp<-setdiff(tmp,stopNodes)
  queue<-tmp
  
  while (length(queue)>0){
    tmp<-edg[[which(nod==queue[1])]]
    tmp<-setdiff(tmp,stopNodes)
    queue<-c(queue, tmp)
    CueDown<-c(CueDown,queue[1])
    queue<-setdiff(queue,CueDown)
    #queue<-queue[-1]
  }

  if (length(grep("and",CueDown))>0){
 	CueDown<-CueDown[-grep("and",CueDown)]
  }
  return(CueDown)
}