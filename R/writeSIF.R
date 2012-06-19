writeSIF <- function(model, filename){

    # convert internal model structure to a SIF matrix
    sif = model2sif(model)

    # and save it into a file
    if (file.exists(filename)==FALSE){
        write.table(sif[,1:3],file=filename,
            row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    }
    else{
       stop(paste("File ", filename, "already exists.",  sep=""))
    }


}

