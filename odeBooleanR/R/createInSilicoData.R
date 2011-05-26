createInSilicoData <-
function(paramInfo,y0,times,cFileName)
{
   if(cFileName=="mymod") 
   {
      out<-lsode(y0, times, func = "derivs", parms=paramInfo$paramValues, jacfunc=NULL, dllname="mymod",initfunc ="initmod",nout = 0);
   }
   else if(cFileName=="myODE")
   {
     out<-cvodes(y0,times,"myODE","rhs",fndata=paramInfo$paramValues,verbose =FALSE,maxstep=0.1);
   }
   else
   {       
      out=NA;
      print("A valid File Name was not provided");
   }
   return(out);
}

