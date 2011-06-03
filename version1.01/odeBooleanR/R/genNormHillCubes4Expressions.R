genNormHillCubes4Expressions <-
function(input,output)
{
  k=paste(input,"_k_",output,sep="");
  n=paste(input,"_n_",output,sep="");
  strODE=paste(input,"^",n,"/(",input,"^",n,"+",k,"^",n,")*(1+",k,"^",n,")",sep="");
  res=list(ODE=str)
  return(strODE);
}

