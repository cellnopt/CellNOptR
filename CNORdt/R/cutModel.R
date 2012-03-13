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

cutModel <- function(Model, SimList, bitString) {
		
	# cut the model according to bitstring	
	ModelCut <- Model
	ModelCut$interMat <- ModelCut$interMat[,as.logical(bitString)]
	ModelCut$notMat <- ModelCut$notMat[,as.logical(bitString)]
	ModelCut$reacID <- ModelCut$reacID[as.logical(bitString)]
	SimListCut <- SimList
	SimListCut$finalCube <- SimListCut$finalCube[as.logical(bitString),]
	SimListCut$ixNeg <- SimListCut$ixNeg[as.logical(bitString),]
	SimListCut$ignoreCube <- SimListCut$ignoreCube[as.logical(bitString),]
	SimListCut$maxIx <- SimListCut$maxIx[as.logical(bitString)]
	dataOut = list(ModelCut, SimListCut)
}