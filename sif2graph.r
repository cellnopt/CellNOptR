sif2graph<-function(sif){
  # build the unique vertices from the column 1 and 3 of the SIF file
  vertices = unique(c(as.character(sif[,1]), as.character(sif[,3])))
  # some aliases
  v1 = sif[,1]
  v2 = sif[,3]
  edges = sif[,2]

  l = length(vertices) - 1
  g <- new("graphNEL", nodes=vertices, edgemode="directed")
  weights = rep(1, l)
  for (i in 1:length(v1)){
    g <- addEdge(as.character(v1[i]), as.character(v2[i]), g, weights[i])
  }
  return(g)
  
}