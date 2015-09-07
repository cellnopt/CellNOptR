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
searchPPIigraph <- function(node1, node2, UniprotID, PPINigraph, noPKNnodes=TRUE){
	
  requireNamespace("igraph")
	
  #node1 is the ending node and is always a 1 value vector
  node1ID<-UniprotID[[which(names(UniprotID)==node1)]]
  #node2 is the starting node and can have more than 1 element (AND gates)
  node2IDs<-vector()
  for (i in 1:length(node2)){
    node2IDs<-c(node2IDs, UniprotID[[which(names(UniprotID)==node2[i])]])
  }
	
  if (noPKNnodes==TRUE){
    #we are looking for paths that do not pass through other nodes of the PKN
    #thus we fist list the nodes of the PKN with all possible Uniprot identifiers..
    PKNnodes<-rapply(UniprotID,c)
    #..exept for the nodes we are interested in..
    noNodes<-setdiff(PKNnodes, c(node1ID,node2IDs))
    #..and we subtract them from the PIN graph (grate a new graph only with the other nodes)
    okNodes<-setdiff(igraph::V(PPINigraph)$name,noNodes)
	# this gives the index of the nodes
	ixOkNodes<-match(okNodes,igraph::V(PPINigraph)$name)
    gg<-igraph::induced.subgraph(PPINigraph,ixOkNodes)
  }else{
    gg<-PPINigraph
  }
  
  ckmin<-rep(Inf,length(node2))
  path<-list()
  #node1 is always a 1 value vector
  #node2 can have more than 1 element
  node1ID<-UniprotID[[which(names(UniprotID)==node1)]]
  for (i in 1:length(node1ID)){
    for (j in 1:length(node2)){
      node2ID<-UniprotID[[which(names(UniprotID)==node2[j])]]
      for (j1 in 1:length(node2ID)){
        if ((node1ID[i] %in% igraph::V(gg)$name) && (node2ID[j1] %in% igraph::V(gg)$name)){
          ix<-igraph::get.shortest.paths(graph=gg,from=node1ID[i],to=node2ID[j1])[[1]]
          ck<-length(ix)
          ck[ck==0]<-Inf
          ckmin[j]<-min(ckmin[j],ck)
        }
      }
    }
  }
  #also the nodes are includes in the path, thus the real length of the path not-including the ending nodes is N-2
  ckmin<-ckmin-2
  #for and nodes the maximum distance from all the nodes is taken
  ckmin<-max(ckmin)
  return(ckmin)
}
