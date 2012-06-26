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
# $Id: $

cutModel <- function(model, simList, bitString) {
		
	# cut the model according to bitstring	
	modelCut <- model
	modelCut$interMat <- modelCut$interMat[,as.logical(bitString)]
	modelCut$notMat <- modelCut$notMat[,as.logical(bitString)]
	modelCut$reacID <- modelCut$reacID[as.logical(bitString)]
	simListCut <- simList
	simListCut$finalCube <- simListCut$finalCube[as.logical(bitString),]
	simListCut$ixNeg <- simListCut$ixNeg[as.logical(bitString),]
	simListCut$ignoreCube <- simListCut$ignoreCube[as.logical(bitString),]
	simListCut$maxIx <- simListCut$maxIx[as.logical(bitString)]
	dataOut = list(modelCut, simListCut)
}