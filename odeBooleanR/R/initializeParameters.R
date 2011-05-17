initializeParameters <-
function(paramInfo,typeInit)
{
  if(typeInit=="randomNK")
  {
    indexNK=c(paramInfo$index_k,paramInfo$index_n)
    paramInfo$paramValues[indexNK]=runif(length(indexNK));
  }
  else if(typeInit=="randomBoundedNK")
  {
    indexK=paramInfo$index_k;
    indexN=paramInfo$index_n;
    paramInfo$paramValues[indexK]=runif(length(indexK),0.5,5);
    paramInfo$paramValues[indexN]=runif(length(indexN),0.5,7);
    paramInfo$UB[indexK]=5;
    paramInfo$LB[indexK]=0.5;
    paramInfo$UB[indexN]=7;
    paramInfo$LB[indexN]=0.5;
  }
  else
  {
    print("Check you have selected any valid Option");
  }
  return(paramInfo);
}

