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

cSimulator <- function(CNOlist, model, simList, indices) {
	
	# CONSTANTS
	# cnolist
	nStimuli <- as.integer(length(indices$stimulated))
	nInhibitors <- as.integer(length(indices$inhibited))
	nSignals <- as.integer(length(indices$signals))
	nCond <- as.integer(dim(CNOlist$valueSignals[[1]])[1])
	nTimes <- as.integer(length(CNOlist$timeSignals))

	# model
	nReacs <- as.integer(length(model$reacID))
	nSpecies <- as.integer(length(model$namesSpecies))

	# simList
	nMaxInputs <- as.integer(dim(simList$finalCube)[2])
	finalCube = as.integer(as.vector(t(simList$finalCube))-1)
	ixNeg = as.integer(as.vector(t(simList$ixNeg)))
	ignoreCube = as.integer(as.vector(t(simList$ignoreCube)))
	maxIx = as.integer(simList$maxIx-1)
	
	# index
	indexSignals <- as.integer(as.vector(indices$signals)-1)
	indexStimuli <- as.integer(as.vector(indices$stimulated)-1)
	indexInhibitors <- as.integer(as.vector(indices$inhibited)-1)

	# STRUCTURES
	# cnolist
	valueInhibitors <- as.integer(t(CNOlist$valueInhibitors))
	valueStimuli <- as.integer(t(CNOlist$valueStimuli))

	res = .Call("CellNOptR",
		# variables	
		nStimuli,
		nInhibitors,
		nSignals,
		nCond,
		nTimes,
		nReacs,
		nSpecies,
		nMaxInputs,
		# vectors
		maxIx,		
		indexSignals, 
		indexStimuli, 
		indexInhibitors,
		# matrices
		finalCube,
		ixNeg,
		ignoreCube,
		valueInhibitors,
		valueStimuli
	)
	
	return(res)
}
