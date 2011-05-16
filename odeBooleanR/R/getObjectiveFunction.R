getObjectiveFunction <-
function(paramInfo1,y01,expData1,times1,obsIndex1,fullpars1,indexPars1,cFileName)
{
  if(cFileName=="mymod")
  {
    objectiveFunction<-function(x, y0=y01,paramInfo=paramInfo1,expData=expData1,times=times1,obsIndex=obsIndex1,fullpars=fullpars1,indexPars=indexPars1)
    {    
       fullpars[indexPars]=x;
       simData<-lsode(y0, times, func = "derivs", parms=fullpars, jacfunc=NULL, dllname="mymod",initfunc ="initmod",nout = 0);
       res=(expData[,obsIndex]-simData[,obsIndex])^2;
      
       return(res);
    }             
  }
  else if(cFileName=="myODE")
  {   
    objectiveFunction<-function(x, y0=y01,paramInfo=paramInfo1,expData=expData1,times=times1,obsIndex=obsIndex1,fullpars=fullpars1,indexPars=indexPars1)
    {    
       fullpars[indexPars]=x;
       simData<-cvodes(y0,times,"myODE","rhs",fndata =fullpars,verbose =FALSE);
       res=sum((expData[,obsIndex]-simData[,obsIndex])^2); 
       if(is.nan(res))res=10000000000;
       print(res);
       return(res);
    }
  }
  else 
  {  
    objectiveFunction<-function(x)
    {
      print("NO VALID OBJECTIVE FUNCTION WAS PROVIDED")
      return(NA);
    }
    print("A valid File Name was not provided");
  }
  
 return(objectiveFunction);   
}

