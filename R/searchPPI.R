#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv2 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-2.0.html
#
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/software.html
#
##############################################################################
# $Id$
searchPPI <- function(node1, node2, UniprotID, graph){
  
  ck<-0
  #node1 is always a 1 value vector
  #node2 can have more than 1 element
  node1ID<-UniprotID[[which(names(UniprotID)==node1)]]
  for (i1 in 1:length(node1ID)){
    for (i11 in 1:length(PINgraphs$g1$connComp)){
      ix1<-grep(pattern=node1ID[i1],graph$connComp[[i11]],perl=TRUE,ignore.case=TRUE)
      if (length(ix1)>0){
        ckVec<-rep(0,length(node2))
        for (j in 1:length(node2)){
          node2ID<-UniprotID[[which(names(UniprotID)==node2[j])]]
          for (j1 in 1:length(node2ID)){
            ix2<-grep(pattern=node2ID[j1],graph$connComp[[i11]],perl=TRUE,ignore.case=TRUE)
            if (length(ix2)>0){
              ckVec[j]<-1
            }
          } 
        }
        if (sum(ckVec)==length(ckVec)){
          ck<-1
        }
      }
    }
  }
  return(ck)
}
