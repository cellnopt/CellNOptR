readMIDAS <-
function(MIDASfile){
	#Read the data
	data<-read.csv(file=MIDASfile,header=TRUE,sep=',',fill=TRUE,as.is=TRUE,check.names=FALSE)
#Determine which are the informative columns (i.e. columns with useful data info and values)
	TRcol<-grep(pattern="TR",x=colnames(data),ignore.case=FALSE)
	DAcol<-grep(pattern="DA",x=colnames(data),ignore.case=FALSE)
	DVcol<-grep(pattern="DV",x=colnames(data),ignore.case=FALSE)
#Print information about the data set and check that the right number of columns are present
	print(paste("Your data set comprises ", nrow(data),"conditions (i.e. combinations of time point and treatment)"))
	if(length(DAcol) != length(DVcol)){
		warning("You have more data values columns (DV columns) than data points columns (DV columns)")
		}
	print(paste("Your data set comprises measurements on ", length(DVcol)," different species"))	
	CellLine<-grep(pattern="(TR:\\w*:CellLine)",x=colnames(data),ignore.case=TRUE,perl=TRUE,value=TRUE)
	if(length(CellLine) != 0){
		CellLine<-sub(pattern="TR:",x=CellLine,replacement="",ignore.case=FALSE)
		CellLine<-sub(pattern=":CellLine",x=CellLine,replacement="",ignore.case=TRUE)
		print(paste("Your data set comprises ", (length(TRcol)-length(CellLine)),"stimuli/inhibitors and", length(CellLine),"cell line(s) (",CellLine,")" ))
		TRcol<-TRcol[-match(grep(pattern="(TR:\\w*:CellLine)",x=colnames(data),ignore.case=TRUE,perl=TRUE,value=FALSE),TRcol)]
		}else{
			print(paste("Your data set comprises ", length(TRcol),"stimuli and inhibitors"))
			warning("There is no cell line information. If some of your TR columns represents the cell lines, please indicate it in your file by naming them 'TR:name:CellLine'")
			}
	print("Please be aware that CNO only handles measurements on one cell line at this time.")
	relCol<-sort(c(TRcol,DAcol,DVcol))
	data<-data[,relCol]
	if(any(is.nan(data))){
		for(c in 1:dim(data)[2]){for(r in 1:dim(data)[1]){if(data[r,c] == "NaN") data[r,c]<-NA}}
		print("Your data file contained 'NaN'. We have assumed that these were missing values and replaced them by NAs.")
		}
	return(list(dataMatrix=data,TRcol=(TRcol-min(relCol)+1),DAcol=(DAcol-min(relCol)+1),DVcol=(DVcol-min(relCol)+1)))
	}

