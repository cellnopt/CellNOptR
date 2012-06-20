readSif <- function(sifFile){
    warning("readSif is a deprecated function. Use readSIF instead. Calling readSIF for you")
    return(readSIF(sifFile))
}

prep4Sim <- function(Model, params){
    warning("readSif is a deprecated function. Use readSIF instead. ")
    return(prep4sim(Model, params))
}
