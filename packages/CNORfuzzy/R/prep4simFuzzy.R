#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI - Massachusetts Institute of Technology
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv2 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-2.0.html
#
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/cno
#
##############################################################################
# $id: prep4simFuzzy.R 1506 2012-06-19 16:09:59Z cokelaer $
prep4simFuzzy <-function(model, paramsList, verbose=TRUE){

    # Get the fields from the CellNOptR package

    fields4sim <- prep4sim(model)

    # now, we can fill fields related to the fuzzy package

    #set default for parameters that aren't likely to be there
    if (!exists('paramsList$defaultg')) {
        paramsList$defaultg = 1
    }
    if (!exists('paramsList$defaultn')) {
        paramsList$defaultn = 3
    }
    if (!exists('paramsList$defaultk')) {
        paramsList$defaultk = 0.5503
    }
    if (!exists('paramsList$type2Def')) {
        paramsList$type2Def = 'Stim'
    }

    #compute indices of stimuli and inhibitors to be used to determine type (this should probably actually be an input)
    Indices<-indexFinder(CNOlist=paramsList$data, model=model,verbose=verbose)

    # reference types according to indices
    typeIDs = array(1,dim(model$interMat)[1]) # same as length(nameSpecies)
    if (grepl('Stim',paramsList$type2Def,ignore.case = TRUE)) {typeIDs[Indices$stimulated] = 2}
    if (grepl('Inhib',paramsList$type2Def,ignore.case = TRUE)) {typeIDs[Indices$inhibited] = 2}

    #Compute the max number of inputs observed in the Model for a single reaction

    # no need to compute it again, it is stored in fields4sim$maxInput
    # maxInput<-colSums(model$interMat)
    # maxInput<-abs(min(maxInput))+1
    # maxInput <- fields4sim$maxInput

    #Make the empty matrices
    gCube<-matrix(1, nrow=length(model$reacID),ncol=fields4sim$maxInput)
    nCube<-matrix(1, nrow=length(model$reacID),ncol=fields4sim$maxInput)
    kCube<-matrix(0, nrow=length(model$reacID),ncol=fields4sim$maxInput)
    typeCube<-matrix(NA, nrow=length(model$reacID),ncol=fields4sim$maxInput)

    #Fill the matrices finalCube, ignoreCube and ixNeg, and maxIx

    for(r in 1:length(model$reacID)){
        input<-which(model$interMat[,r] == -1)
        gCube[r,1:length(input)]<- paramsList$defaultg
        nCube[r,1:length(input)]<- paramsList$defaultn
        kCube[r,1:length(input)]<- paramsList$defaultk
        typeCube[r,1:length(input)]<- typeIDs[input]
    }
    reshapeIx <- which(!fields4sim$ignoreCube)
    reshapeType1 <- which(typeCube==1)
    reshapeType2 <- which(typeCube==2)
    numType1 <- length(reshapeType1)
    numType2 <- length(reshapeType2)

    fields4sim$gCube=gCube
    fields4sim$nCube=nCube
    fields4sim$kCube=kCube
    fields4sim$typeCube=typeCube
    fields4sim$reshapeIx=reshapeIx
    fields4sim$reshapeType1=reshapeType1
    fields4sim$reshapeType2=reshapeType2
    fields4sim$numType1=numType1
    fields4sim$numType2=numType2

    return(fields4sim)
}

