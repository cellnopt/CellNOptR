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
# $Id: cSimulator.R 3942 2013-08-28 14:16:59Z cokelaer $

cSimulator <- function(CNOlist, model, simList) {

#    if ((class(CNOlist)=="CNOlist")==FALSE){
#        CNOlist = CellNOptR::CNOlist(CNOlist)
#    }
#	
#	# check the structures
#	if(is.null(CNOlist@stimuli) || is.null(CNOlist@inhibitors)) {
#		stop("This function needs 'valueStimuli' and 'valueInhibitors' in CNOlist")
#	}
#	
#	if(is.null(model$reacID) || is.null(model$namesSpecies)) {
#		stop("This function needs 'reacID' and 'namesSpecies' in model")
#	}
#
        indexList<-indexFinder(CNOlist=CNOlist,model=model, verbose=FALSE)
	# variables
	nCond <- as.integer(dim(CNOlist@stimuli)[1])
    nReacs <- as.integer(length(model$reacID))
	nSpecies <- as.integer(length(model$namesSpecies))

    if (is.null(dim(simList$finalCube))){
    	nSpecies <- length(model$interMat)
        maxIx <- as.integer(dim(simList$finalCube)[2])
    }

	# simList
	finalCube = as.integer(simList$finalCube-1)
	ixNeg = as.integer(simList$ixNeg)
	#ignoreCube = as.integer(as.vector(t(simList$ignoreCube)))
	ignoreCube = as.integer(simList$ignoreCube)
	maxIx = as.integer(simList$maxIx-1)
	
	# index
	indexSignals <- as.integer((indexList$signals)-1)
	indexStimuli <- as.integer((indexList$stimulated)-1)
	indexInhibitors <- as.integer((indexList$inhibited)-1)
    nSignals <- length(indexSignals)

	# cnolist changed to float in bug report 222
	valueInhibitors <- as.numeric(CNOlist@inhibitors)
	valueStimuli <- as.numeric(CNOlist@stimuli)


	res = .Call("simulatorT1",
		# variables	
	    as.integer(length(indexList$stimulated)),
	    as.integer(length(indexList$inhibited)),
		nCond,
		nReacs,
		nSpecies,
        nSignals,
		as.integer(simList$maxInput),
		# simList
		finalCube,
		ixNeg,
		ignoreCube,
		maxIx,	
        as.double(as.vector(simList$gCube)),
        as.double(as.vector(simList$kCube)),
        as.double(as.vector(simList$nCube)),
		# index
		indexSignals, 
		indexStimuli, 
		indexInhibitors,
		# cnolist
		valueInhibitors,
		valueStimuli
	)

	return(res)
}
