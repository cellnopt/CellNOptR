#
#  Copyright 2011 EBI
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
# $Id$
cutSimList <- function(SimList, bitString)
{
    SimListCut<-SimList
    finalCube<-SimListCut$finalCube[as.logical(bitString),]
    ixNeg<-SimListCut$ixNeg[as.logical(bitString),]
    ignoreCube<-SimListCut$ignoreCube[as.logical(bitString),]
    maxIx<-SimListCut$maxIx[as.logical(bitString)]
    # in some cases the finalcube is a matrix but list of integer, so we
    # need to convert back to a matrix. Happens for simple models only hence
    # the warning.
    if (is.matrix(finalCube) == FALSE){
        warning("converting back to matrix in prep4sim")
        SimListCut$finalCube<-matrix(finalCube,
            dimnames=list(names(finalCube), 1))
        SimListCut$ixNeg<-matrix(ixNeg, dimnames=list(names(ixNeg), 1))
        SimListCut$ignoreCub<-matrix(ignoreCube,dimnames=list(names(ignoreCube), 1))
        SimListCut$maxIx<-matrix(maxIx,dimnames=list(names(maxIx), 1))
    }
    else{
        SimListCut$finalCube = finalCube
        SimListCut$ixNeg<-ixNeg
        SimListCut$ignoreCube<-ignoreCube
        SimListCut$maxIx<-maxIx
    }
    return(SimListCut)
}
