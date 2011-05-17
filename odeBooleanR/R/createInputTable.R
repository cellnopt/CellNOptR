createInputTable <-
function(numInputs)
{ 
  mat=c();
  for(i in 1:2^numInputs)mat=rbind(mat,as.logical(dec2bin(i-1,numInputs)));
  return(mat);
}

