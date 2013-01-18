#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EMBL - European Bioinformatics Institute
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/cno
#
##############################################################################
# $Id$

## given a MIDAS file, returns a CNOlist object
## You may want to make this (and the existing readMIDAS)
## into S4 generics in order to avoid a name conflict.





## Class definition
setClass("CNOlist",
    representation(
        cues="matrix",
        inhibitors="matrix",
        stimuli="matrix",
        signals="list",
        variances="list",
        timepoints="vector"),
    ## validity method
    validity=function(object) {
        msg <- NULL
        nrow <- nrow(cues(object))
        signalNrows <- unique(sapply(signals(object), nrow))
        if (nrow != nrow(inhibitors(object)) ||
            nrow != nrow(stimuli(object)) ||
            length(signalNrows) != 1 ||
            nrow != signalNrows)
        {
            msg <- "'nrow' differs between elements"
        }
        if (is.null(msg)) TRUE else msg
    })


## constructor
CNOlist <-function(data, subfield=FALSE, verbose=FALSE){

    # input can a filename or the old (still used) CNOlist returned by
    # makeCNOlist function. subfield and verbose used only if MIDASfile is a
    # string.
    if (is.character(data)== TRUE){
        res = internal_CNOlist_from_file(data, subfield, verbose)
    }else {
        if (is.list(data)==TRUE){
            if ("namesCues" %in% names(data) == TRUE){
                res = internal_CNOlist_from_makeCNOlist(data) 
            }else{
                stop("Not a valid list. Does not seem to be returned by CellNOptR::makeCNOlist")
            }
        }else{
        stop("Input data must be a filename or the output of CellNOptR::makeCNOlist function")
        }
    }

    new("CNOlist", cues=res$cues, inhibitors=res$inhibitors,
        stimuli=res$stimuli, signals=res$signals, variances=res$variances, timepoints=res$timepoints)
}



## accessors
cues <- function(cnoList, ...) cnoList@cues
inhibitors <- function(cnoList, ...) cnoList@inhibitors
stimuli <- function(cnoList, ...) cnoList@stimuli
signals <- function(cnoList, ...) cnoList@signals
variances <- function(cnoList, ...) cnoList@variances
timepoints <- function(cnoList, ...) cnoList@timepoints
#plot <- function(cnoList, ...) plotCNOlist(cnoList)


## show method
setMethod(show, "CNOlist", function(object) {
    cat("class:", class(object), "\n")
    cat("cues:", colnames(cues(object)), "\n")
    cat("inhibitors:", colnames(inhibitors(object)), "\n")
    cat("stimuli:", colnames(stimuli(object)), "\n")
    cat("timepoints:", names(signals(object)), "\n")
    cat("signals:", colnames(signals(object)[[1]]), "\n")
    cat("variances:", colnames(signals(object)[[1]]), "\n")
    cat("--\nTo see the values of any data contained in this instance, just use the
method (e.g., cues(cnolist), signals(cnolist), variances(cnolist), ...\n\n")
})

#setMethod("plot", signature(x="CNOlist", y="missing"), function(x, y, ...){
#    plotCNOlist(x)
#})

setMethod("plot", "CNOlist", function(x, y, ... ){
    plotCNOlist(x)
})
setMethod("plot", signature(x="CNOlist", y="CNOlist"), function(x, y, ... ){
    plotCNOlist2(x,y)
})

setMethod("length", "CNOlist", function(x) length(x@timepoints))

#if (isGeneric("randomize")==FALSE){
    setGeneric(
        name="randomize",
        def=function(object,sd=0.1, minValue=0,maxValue=1,mode="gaussian"){standardGeneric("randomize")}
    )
#}
#lockBinding("randomize", .GlobalEnv)


setMethod("randomize", "CNOlist", 
    definition=function(object, sd=0.1, minValue=0, maxValue=1,mode="uniform"){
        res = randomizeCNOlist(object, sd=sd, mode=mode)
        return(res)
    }
)


# used by the constructor not for export.
# open a MIDAS file (given a  filename) and create the instance of CNOlist
internal_CNOlist_from_file <- function(MIDASfile, subfield=FALSE, verbose=FALSE)
{
    x <- readMIDAS(MIDASfile, verbose=verbose)
    cnolist <- makeCNOlist(x, subfield=subfield, verbose=verbose)
    res <- internal_CNOlist_from_makeCNOlist(cnolist)
    return(res)
}


# used by the constructor not for export.
# open a MIDAS file (given a  filename) and create the instance of CNOlist
internal_CNOlist_from_makeCNOlist <- function(cnolist)
{

    myCues <- cnolist$valueCues
    colnames(myCues) <- cnolist$namesCues

    myInhibitors <- cnolist$valueInhibitors
    colnames(myInhibitors) <- cnolist$namesInhibitors

    myStimuli <- cnolist$valueStimuli
    colnames(myStimuli) <- cnolist$namesStimuli

    mySignals <- cnolist$valueSignals
    names(mySignals) <- cnolist$timeSignals
    mySignals <- lapply(mySignals, "colnames<-", cnolist$namesSignals)

    if ("valueVariances" %in% names(cnolist)){
        myVars <- cnolist$valueVariances
        names(myVars) <- cnolist$timeSignals
        myVars <- lapply(myVars, "colnames<-", cnolist$namesSignals)
    } else{
        myVars = mySignals
        for (time in 1:length(mySignals)){
            myVars[[time]] = myVars[[time]] * NA
        }
    }

    myTimePoints <- cnolist$timeSignals

    #CNOlist(myCues, myInhibitors, myStimuli, mySignals)
    return( list(cues=myCues, inhibitors=myInhibitors, stimuli=myStimuli,
        signals=mySignals, variances=myVars, timepoints=myTimePoints))
}


