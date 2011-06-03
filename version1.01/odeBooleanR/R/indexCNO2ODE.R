indexCNO2ODE <-
function(CNOlist,odePars)
{ 
  CNOindexSignals=c();
  ODEindexSignals=c();
  for(i in 1:length(odePars$outputs))
  {
    tempIndex=which(CNOlist$namesSignals==odePars$outputs[i]);
    CNOindexSignals=c(CNOindexSignals,tempIndex);
    if(length(tempIndex)>0)ODEindexSignals=c(ODEindexSignals,i);
  }
  
  CNOindexStimuli=c();
  ODEindexStimuli=c();
  for(i in 1:length(odePars$inputs))
  {
    tempIndex=which(CNOlist$namesStimuli==odePars$inputs[i]);
    CNOindexStimuli=c(CNOindexStimuli,tempIndex);
    if(length(tempIndex)>0) ODEindexStimuli=c(ODEindexStimuli,i);
  }
  
  CNOindexInhibitors=c();
  ODEindexInhibitors=c();
  for(i in 1:length(odePars$outputs))
  {
    tempIndex=which(CNOlist$namesInhibitors==odePars$outputs[i]);
    CNOindexInhibitors=c(CNOindexInhibitors,tempIndex);
    if(length(tempIndex)>0)ODEindexInhibitors=c(ODEindexInhibitors,i);
  }
  res=list(CNOindexSignals=CNOindexSignals,ODEindexSignals=ODEindexSignals,
    CNOindexStimuli=CNOindexStimuli,ODEindexStimuli=ODEindexStimuli,
      CNOindexInhibitors=CNOindexInhibitors,ODEindexInhibitors=ODEindexInhibitors);
  
  return(res);
}

