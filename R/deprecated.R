readSif <- function(sifFile){
    warning("readSif is a deprecated function. Use readSIF instead. Calling readSIF for you")
    return(readSIF(sifFile))
}

prep4Sim <- function(model, params){
    warning("prep4Sim is a deprecated function. Use prep4sim instead. ")
    return(prep4sim(model))
}
