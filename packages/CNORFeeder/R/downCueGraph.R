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
downCueGraph <-
function(cue,graph,stopNodes){
  
  gg<-graph
  
  nod<-graph::nodes(gg)
  edg<-graph::edges(gg)
  
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
