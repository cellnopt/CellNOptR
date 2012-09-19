# This file is part of the CNO software
# 
# Copyright (c) 2011-2012 - EBI
# 
# File author(s): CNO developers (cno-dev@ebi.ac.uk)
# 
# Distributed under the GPLv2 License.  See accompanying file LICENSE.txt or copy at
# http://www.gnu.org/licenses/gpl-2.0.html
# 
# CNO website: http://www.ebi.ac.uk/saezrodriguez/software.html
# 
# $Id$
    	
    	
    	convert2array <- function(x, nRow, nCol, nBool) {
    v1 = c(x)
    count = 1
    out1 = array(NA, dim = c(nRow, nCol, nBool))
    for (d in 1:nBool) {
        for (a in 1:nRow) {
            for (b in 1:nCol) {
                out1[a, b, d] = v1[count]
                count = count + 1
            }
        }
    }
    
    return(out1)
	}