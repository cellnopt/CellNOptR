#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2013 - EBI
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

# randomize a cnolist by adding noise to the signal matrices at all times except time=0
# The random noise is a gaussian distribution, which standard deviation can be
# changed.
# The random noise can also be a uniform noise between 0 and 1, in which case
# the data is not used anymore.
randomizeCNOlist <- function(cnolist, sd=0.1, minValue=0, maxValue=1, mode="gaussian"){

     if (mode %in% list("gaussian", "uniform")==FALSE){
         stop("mode can be only 'gaussian'")
     }

     cnames = colnames(cnolist@signals$`0`)
     dimensions = dim(cnolist@signals$`0`)
     times = cnolist@timepoints


     for (time in 2:length(times)){
         if (mode=="gaussian"){
             cnolist@signals[[time]] = rnorm(cnolist@signals[[time]], sd=sd)
         }
         if (mode=="uniform"){
             cnolist@signals[[time]] = runif(cnolist@signals[[time]])
         }
         cnolist@signals[[time]][cnolist@signals[[time]]<0] <- 0
         cnolist@signals[[time]][cnolist@signals[[time]]>1] <- 1
         dim(cnolist@signals[[time]])<-dimensions
         colnames(cnolist@signals[[time]])<-cnames
    }

     return(cnolist)
 }

