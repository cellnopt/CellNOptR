genNormHillCubes <-
function(input,output)
{
  k=paste(input,"_k_",output,sep="");
  n=paste(input,"_n_",output,sep="");
  strODE=paste("pow(",input,",",n,")/(pow(",input,",",n,")+pow(",k,",",n,"))*(1+pow(",k,",",n,"))",sep="");
  res=list(ODE=str)
  return(strODE);
}

