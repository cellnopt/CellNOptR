readSif <- function(sifFile){
    warning("readSif is a deprecated function. Use readSIF instead. Calling readSIF for you")
    return(readSIF(sifFile))
}

prep4Sim <- function(model, params){
    warning("prep4Sim is a deprecated function. Use prep4sim instead. ")
    return(prep4sim(model))
}

simulateT1 <- function(CNOlist, model, bStringT1,simList, indexList){
    warning("simulateT1 is a deprecated function. Use simulate instead. ")
    return(simulate(CNOlist, model, bStrings=list(bStringT1)))
}

computeScoreT2 <- function(CNOlist, model, simList=NULL, indexList=NULL, 
    simResT1, bStringT1, bStringT2, sizeFac=0.0001, NAFac=1){

    warning("computeScoreT2 is a deprecated function. Use computeScoreTN instead.")
    return(computeScoreTN(CNOlist, model, simList=simList, indexList=indexList, 
        simResT1, bStringT1,  bStringT2, timeIndex=3, sizeFac=sizeFac,
        NAFac=NAFac))
 
}


