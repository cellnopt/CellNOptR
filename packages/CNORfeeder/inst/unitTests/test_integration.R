test_integration <- function() {
  # load the data already formatted as CNOlist
  data(CNOlistDREAM,package="CellNOptR")
  # load the model (PKN) already in the CNO format
  data(DreamModel,package="CellNOptR")
  # see CellNOptR documentation to import other data/PKNs)
  
  if ((class(CNOlistDREAM)=="CNOlist")==FALSE){
    CNOlistDREAM = CellNOptR::CNOlist(CNOlistDREAM)
  }	
  
  
  BTable <- makeBTables(CNOlist=CNOlistDREAM, k=2, measErr=c(0.1, 0))
  
  model<-preprocessing(data=CNOlistDREAM, model=DreamModel)
  
  modelIntegr <- mapBTables2model(BTable=BTable,model=model,allInter=TRUE)
  
  checkTrue(is.vector(modelIntegr$indexIntegr))

}