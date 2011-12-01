SearchLink <-
function(node1,node2,Model,noNodes=vector(), ck=0){
  # looks in the model if there is a connection between node1 and node2
  # giving as output 1 if they are connected and 0 otherwise
  
  # 1. finds reactions in which the investigated node2 is the target node
  ix <- which(Model$interMat[rownames(Model$interMat)==node2] == 1)
  if (ck != 1 && length(ix) > 0) {
    # 2. finds which species are source nodes in those reactions
    index <- unique(as.matrix(which(Model$interMat[,ix]==-1, arr.ind=TRUE))[,1])
    # 3. if one of those species is node1, ck is set to 1 (we found the connection)
    if (sum(rownames(Model$interMat)[index] == node1) > 0) {
      ck <- 1
      return(ck)
    }
	else {
      for (i in 1:length(index)){
        if (ck != 1){
		  # 4. if one of the noNodes is met (e.g. an inhibitor that we now is not affecting
		  #    the protein) the search in that particoular path stops
		  if (sum(noNodes == rownames(Model$interMat)[index[i]]) > 0) {
				
		  }	
		  # 5. otherwise the searching procedure is repeated for each upstream node
		  else {
          node2 <- rownames(Model$interMat)[index[i]]
          ck <- SearchLink(node1 = node1, node2 = node2, Model = Model, noNodes=noNodes, ck = ck)
		  }
        }
      }
    }
  }
  return(ck)
}