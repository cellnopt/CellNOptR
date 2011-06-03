computeSensitivities <-
function(mathModel,pDerivsInfo,CNOlist,indexes,times)
{
  pderivs=pDerivsInfo$listPDeriv;
  numParameters=length(mathModel$index_k)+length(mathModel$index_n)+
      length(mathModel$index_k);
  numStates=length(mathModel$expressions);
  numExperiments=length(CNOlist$valueStimuli[,1]);
  simData=simulateModel(CNOlist,indexes,mathModel,times);
  numInputs=length(mathModel$inputs);
  numPderivs=numParameters*numStates;
  res=list();
  
  for(i in 1:length(mathModel$index_inh))assign(mathModel$parNames[mathModel$index_inh[i]],0);
  for(i in 1:numParameters)assign(mathModel$parNames[i],mathModel$paramValues[i]);  
  for(i in 1:numExperiments)
  { 
    expRes=vector('list',numPderivs); 
    for(j in 1:numInputs)
    {
      assign(mathModel$inputs[j],CNOlist$valueStimuli[i,j]); 
    }
    
    for(j in indexes$CNOindexInhibitors)
    {
      print(j) 
    }
    index[i]=1;
    
    for(j in 1:length(times))
    {
      for(k in 1:numStates)
      {
        assign(mathModel$outputs[k],simData[[i]][j,k]); 
      }
      count=0;
      for(k in 1:numStates)
      {
        for(r in 1:numParameters)
        {
          count=count+1;
          if(pderivs[[count]]!=0)
          {
            expRes[[count]]=c(expRes[[count]],eval(pderivs[[count]]))
          }
          else{expRes[[count]]=0;}
        }
      }
    }
    res[[i]]=expRes
  }
  return(res)
}

