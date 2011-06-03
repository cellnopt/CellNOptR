simulateData <-
function(CNOlist,indexes,odePars)
{
  numExperiments=length(CNOlist$valueStimuli[,1]);
  times=CNOlist$timeSignals;
  res=list();
  for(i in 1:numExperiments)
  {
    odePars$paramValues[odePars$index_inputs[indexes$ODEindexStimuli]]=
      CNOlist$valueStimuli[i,indexes$CNOindexStimuli];
    odePars$paramValues[odePars$index_inh[indexes$ODEindexInhibitors]]=
      CNOlist$valueInhibitors[i,indexes$CNOindexInhibitors];
      y0=matrix(0,1,length(odePars$outputs));
    y0[indexes$ODEindexSignals]=CNOlist$valueSignals[[1]][i,indexes$CNOindexSignals]
      CNOlist$valueInhibitors[i,indexes$CNOindexInhibitors];
    tempRes=createInSilicoData(odePars,y0,times,"myODE");
    res[[i]]=tempRes[,indexes$ODEindexSignals];    
  }
  
  return(res);
}

