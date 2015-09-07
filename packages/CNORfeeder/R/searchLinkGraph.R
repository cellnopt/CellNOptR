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
searchLinkGraph <-
function(node1,node2,graph,noNodes=vector()){
  # looks in the graph if there is a connection between node1 and node2
  # giving as output 1 if they are connected and 0 otherwise
  
  #detete from the graph the nodes through which we cannot pass
  if (length(noNodes>0)){
    noNodes<-intersect(noNodes,graph::nodes(graph))
    gg<-graph::removeNode(noNodes,graph)
  }else{
	  gg<-graph
  }
  
  nod<-graph::nodes(gg)
  edg<-graph::edges(gg)
  
  if (any(which(nod==node1))){
    queue<-edg[[which(nod==node1)]]
	gg<-graph::removeNode(node1,gg)
	nod<-graph::nodes(gg)
	edg<-graph::edges(gg)
  }else{
    queue<-vector()
  }
  
  ck<-0
  while (length(queue)>0 && ck!=1){
	if (any(which(nod==queue[1]))){
      queue<-c(queue, edg[[which(nod==queue[1])]])
	}
    if(any(queue==node2)){
      ck<-1
    }
	if (any(which(nod==queue[1]))){
	  gg<-graph::removeNode(queue[1],gg)
	  nod<-graph::nodes(gg)
	  edg<-graph::edges(gg)
	}
	queue<-queue[-1]
  }
  return(ck)
}
