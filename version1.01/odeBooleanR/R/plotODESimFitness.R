plotODESimFitness <-
function(SimResults,CNOList)
{ 
  times=CNOlist$timeSignals;
  expResults=CNOlist$valueSignals;  
  SimResults=convertTimeSeriesToList(SimResults,CNOlist);
  
  namesSignals=CNOlist$namesSignals;
  
  namesCues=c(CNOlist$namesStimuli,CNOlist$namesInhibitors);
  
  valueCues=cbind(CNOlist$valueStimuli,CNOlist$valueInhibitors)
  valueCues=as.matrix(valueCues);
  valueCues[which(valueCues>0)]=1;
  names(valueCues)=namesCues;
  
  plotOptimResults(SimResults=SimResults,expResults=expResults,times=times,namesCues=namesCues,namesSignals=namesSignals,valueCues=valueCues); 
}

