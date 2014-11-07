test_inference <- function() {
  # load the data already formatted as CNOlist
  data(CNOlistDREAM,package="CellNOptR")
  # load the model (PKN) already in the CNO format
  data(DreamModel,package="CellNOptR")
  # see CellNOptR documentation to import other data/PKNs)
  
  if ((class(CNOlistDREAM)=="CNOlist")==FALSE){
    CNOlistDREAM = CellNOptR::CNOlist(CNOlistDREAM)
  }	
  
  
  BTable <- makeBTables(CNOlist=CNOlistDREAM, k=2, measErr=c(0.1, 0))
  
  checkEquals(length(BTable$namesSignals), length(colnames(CNOlistDREAM@signals[[2]])))
  checkEquals(length(BTable$tables), length(colnames(CNOlistDREAM@signals[[2]])))
}