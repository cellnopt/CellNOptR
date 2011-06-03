getPartialDerivs <-
function(mathModel)
{
  numParameters=length(mathModel$index_k)+length(mathModel$index_n)+
      length(mathModel$index_k);
  numStates=length(mathModel$expressions);
  listPDeriv=list();
  count=0;
  indexNonZeroPDeriv=c();
  tag=c();
  for(i in 1:numStates)
  {
    for(j in 1:numParameters)
    {
      count=count+1;
      tag=c(tag,paste('d',mathModel$outputs[i],'/d',mathModel$parNames[j],sep=""));
      listPDeriv[[count]]=D(mathModel$expressions[i],mathModel$parNames[j]);
      if(listPDeriv[[count]]!=0)indexNonZeroPDeriv=c(indexNonZeroPDeriv,count);
    }
  }
  res=list(listPDeriv=listPDeriv,tag=tag,indexNonZeroPDeriv=indexNonZeroPDeriv);
  return(res)
}

