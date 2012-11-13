#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GpLv2 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-2.0.html
#
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/software.html
#
##############################################################################
# $Id$

feedbackFinder1 <- function(model) {

v = 1:length(model$namesSpecies)
gamma = vector("list",length(v))
names(gamma) = model$namesSpecies
for(a in 1:length(v)) {
	pickCols = which(model$interMat[v[a],]==-1)
	for(b in pickCols) {
		gamma[[a]] = c(gamma[[a]], which(model$interMat[,b]==1))
	}
	if(!is.null(gamma[[a]])) {
		gamma[[a]] = unique(gamma[[a]])
	}
}

h = matrix(0, nrow=length(model$namesSpecies), ncol=length(model$namesSpecies))
k = 1
p = 0
p[1] = v[1]
flagT = 1
flagUpdate = c()
j = 1
loops = list()

while (flagT == 1) {

	if(!is.null(gamma[[p[k]]])) {
		if(length(setdiff(gamma[[p[k]]],h[p[k],]))) {
			whatDiff = setdiff(gamma[[p[k]]],h[p[k],])
			openEdges = c() 
			for(a in 1:length(whatDiff)) {
				openEdges = c(openEdges, which(gamma[[p[k]]] == whatDiff[a]))
			}
		} else {
			openEdges = 1	
		}
	
		if (gamma[[p[k]]][openEdges[1]] > p[1] && !any(p == gamma[[p[k]]][openEdges[1]]) && !any(h[p[k],] == gamma[[p[k]]][openEdges[1]])) {
			flagUpdate = 1	
		} else {
			flagUpdate = 0	
		}
	} else {
		openEdges = 1
		flagUpdate = 0
	}
	
	while (flagUpdate == 1) {
		whatDiff = setdiff(gamma[[p[k]]],h[p[k],])
		openEdges = c() 
		for(a in 1:length(whatDiff)) {
			openEdges = c(openEdges, which(gamma[[p[k]]] == whatDiff[a]))
		}
		if (gamma[[p[k]]][openEdges[1]] > p[1] && !any(p == gamma[[p[k]]][openEdges[1]]) && !any(h[p[k],] == gamma[[p[k]]][openEdges[1]])) {
			k = k+1
			p[k] = gamma[[p[k-1]]][openEdges[1]]
			if(!is.null(gamma[[p[k]]])) {
				flagUpdate = 1
			} else {
				flagUpdate = 0	
			}	
		} else {
			flagUpdate = 0	
		}	
	}
	
	# check to report circuit here
	if(!is.null(gamma[[p[k]]])) {
		if(p[1] == gamma[[p[k]]][openEdges[1]]) {
			loops = c(loops,list(p))	
		}
	}
	
	# check for closure here
	if(k==1) {
		if(p[1] == length(model$namesSpecies)) {
			# terminate
			flagT = 0	
		} else {
			# go to next node
			p[1] = p[1] + 1
			k = 1
			h = matrix(0, nrow=length(model$namesSpecies), ncol=length(model$namesSpecies))
		}
		
	} else {
		m = which(h[p[k-1],] == 0)[1]
		h[p[k],] = 0
		h[p[k-1],m] = p[k]
		p[k] = 0			
		k = k-1
	}
}

return(loops)
}



