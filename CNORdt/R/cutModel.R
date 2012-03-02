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