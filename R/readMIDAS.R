readMIDAS<-function(MIDASfile, verbose=TRUE){

	#Read the data
	data<-read.csv(file=MIDASfile,header=TRUE,sep=',',fill=TRUE,as.is=TRUE,check.names=FALSE)

#Determine which are the informative columns (i.e. columns with useful data info and values)
	TRcol<-grep(pattern="TR",x=colnames(data),ignore.case=FALSE)
	DAcol<-grep(pattern="DA",x=colnames(data),ignore.case=FALSE)
	DVcol<-grep(pattern="DV",x=colnames(data),ignore.case=FALSE)

#Print information about the data set 
	if (verbose){
        print(paste(
		    "Your data set comprises ", nrow(data),
    		"conditions (i.e. combinations of time point and treatment)"))
    }

#check that the right number of columns are present			
	if(length(DAcol) != length(DVcol)){

		if(length(grep("DA:ALL",colnames(data))) != 0 && length(DAcol) == 1){

			if (verbose){
                print("Your data comprises only a DA:ALL column;  all readouts are assumed to have been acquired at the same time.")
            }

			#make the additional columns
			newDA<-matrix(rep(data[,DAcol],length(DVcol)),byrow=FALSE,ncol=length(DVcol))
			newDAnames<-colnames(data)[DVcol]
			colnames(newDA)<-sub(pattern="DV:",x=newDAnames,replacement="DA:",ignore.case=FALSE)
			data<-cbind(data,newDA)
			data<-data[,-DAcol]

			}else{
				warning("You have more data values columns (DV columns) than data points columns (DV columns)")
				}
		}

	if (verbose){
		print(paste("Your data set comprises measurements on ", length(DVcol)," different species"))
	}

	CellLine<-grep(pattern="(TR:\\w*:CellLine)",
		x=colnames(data),ignore.case=TRUE,perl=TRUE,value=TRUE)

	if(length(CellLine) != 0){

		CellLine<-sub(pattern="TR:",x=CellLine,replacement="",ignore.case=FALSE)
		CellLine<-sub(pattern=":CellLine",x=CellLine,replacement="",ignore.case=TRUE)

		if (verbose){
            print(paste(
			    "Your data set comprises ", (length(TRcol)-length(CellLine)),
                "stimuli/inhibitors and", length(CellLine),"cell line(s) (",CellLine,")" ))
        }

		data<-data[,-grep(pattern="(TR:\\w*:CellLine)",x=colnames(data),ignore.case=TRUE,perl=TRUE,value=FALSE)]

		}else{

		if (verbose){
			print(paste("Your data set comprises ", length(TRcol),"stimuli and inhibitors"))
        }
			warning("There is no cell line information. If some of your TR columns represents the cell lines, please indicate it in your file by naming them 'TR:name:CellLine'")

		}
	if (verbose){	
        print("Please be aware that CNO only handles measurements on one cell line at this time.")
    }

	TRcol<-grep(pattern="TR",x=colnames(data),ignore.case=FALSE)
	DAcol<-grep(pattern="DA",x=colnames(data),ignore.case=FALSE)
	DVcol<-grep(pattern="DV",x=colnames(data),ignore.case=FALSE)
	data<-data[,c(TRcol,DAcol,DVcol)]


    # replace NaN character by NA. as.matrix is required to scan all columns AND rows
    # Note that on some older R version, the is.nan does not seem to work well, hence the
    # try catch (see Changelog 0.99.6
    conversion <- tryCatch({
        if (any(is.nan(as.matrix(data)))){
            data[is.nan(as.matrix(data))] <- NA
        }
        conversion = TRUE},
        error=function(e) {return(FALSE)},
        finally={
            if (verbose){
                print("Your data file contained 'NaN'. We have assumed that these were missing values and replaced them by NAs.")
            }}
    )

    if ( conversion == FALSE){
        if(any(as.matrix(data == "NaN"))){
			for(c in 1:dim(data)[2]){
	            for(r in 1:dim(data)[1]){
                     if(data[r,c] == NA) {
                         data[r,c]<-NA
                     }
                }
		    }
            if (verbose){
                print("Your data file contained 'NaN'. We have assumed that these were missing values and replaced them by NAs.")
            }
	    }
    }


	return(list(
		dataMatrix=data,
		TRcol=grep(pattern="TR",x=colnames(data),ignore.case=FALSE),
		DAcol=grep(pattern="DA",x=colnames(data),ignore.case=FALSE),
		DVcol=grep(pattern="DV",x=colnames(data),ignore.case=FALSE)))

	}

