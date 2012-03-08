graph2sif<-function(graph,writeSif=FALSE){
  
  ed<-edges(graph)
  no<-nodes(graph)
  sifFile<-matrix(ncol=3)
  
  for (i in 1:length(no)){
    node1<-no[i]
    for (j in 1:length(ed[[i]])){
      node2<-ed[[i]][j]
      sifFile<-rbind(sifFile,c(node1,1,node2))
    }
  }
  
  sifFile<-sifFile[-1,]
  
  if (writeSif==TRUE){
    write.table(sifFile, file="Graph.sif", row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
  }

  return(sifFile)
}