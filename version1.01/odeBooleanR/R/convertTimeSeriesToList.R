convertTimeSeriesToList <-
function(sim,CNOlist)
{
  nTimePoints=length(data.frame(sim[1])[,1]);
  nSpecies=length(data.frame(sim[1])[1,])
  nExperiments=length(sim);
  res=list();
  for(i in 1:nTimePoints)res[[i]]=data.frame();
  for(i in 1:nExperiments)
  {
    temp=data.frame(sim[i]);
    for(j in 1:nTimePoints)
    {
      res[[j]]=rbind(res[[j]],temp[j,]);
    } 
  }
  for(i in 1:nTimePoints)colnames(res[[i]])=CNOlist$namesSignals;
  return(res)  
}

