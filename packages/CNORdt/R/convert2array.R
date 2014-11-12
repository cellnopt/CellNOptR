# This file is part of the CNO software
# 
# Copyright (c) 2011-2013 - EBI
# 
# File author(s): CNO developers (cno-dev@ebi.ac.uk)
# 
# Distributed under the GPLv3 License.  See accompanying file LICENSE.txt or copy at
# http://www.gnu.org/licenses/gpl-3.0.html
# 
# CNO website: http://www.cellnopt.org
# 
# $Id$
    	
convert2array <- function(x, nCond, nSpecies, nBool) {
  v1 = c(x)
  count = 1
  out1 = array(NA, dim = c(nCond, nSpecies, nBool))
  for (d in 1:nBool) {
    for (a in 1:nCond) {
      for (b in 1:nSpecies) {
        out1[a, b, d] = v1[count]
        count = count + 1
      }
    }
  }
  
  return(out1)
}
