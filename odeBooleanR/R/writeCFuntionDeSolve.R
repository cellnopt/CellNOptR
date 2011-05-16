writeCFuntionDeSolve <-
function(inputs,outputs,strODE,numOutputs,parVecForPrint)
{
   outputsStr=paste(outputs,sep=",",collapse=",");
   outputsStr=paste("double ",outputsStr,";",sep="");
   
   outputAttr=c();
    
   for(i in 1:numOutputs)
   {
      outputAttr[i]=paste("double ",outputs[i],"=y[",i-1,"];",sep="");
   }
   
   for(i in 1:numOutputs)
   {
      strODE[i]=paste("ydot[",i-1,"]=",strODE[i],";",sep="");
   }
       
   parVecForPrintInclude=c();    
   
   for(i in 1:length(parVecForPrint))
   {
     parVecForPrintInclude[i]=paste("#define ",parVecForPrint[i]," parms[",i-1,"]",sep="");
   }
   
   write("#include <R.h>",file="mymod.c",sep="",append=FALSE);
   write("#include <math.h>",file="mymod.c",sep="",append=TRUE);
   tempStr=paste("static double parms[",length(parVecForPrintInclude),"];",sep="");
   write(tempStr,file="mymod.c",sep="",append=TRUE)
   write(parVecForPrintInclude,file="mymod.c",sep="",append=TRUE);
   tempStr=paste("void initmod(void (* odeparms)(int *, double *))\n{\n int N=",length(parVecForPrint),";",sep="")
   write(tempStr,file="mymod.c",sep="",append=TRUE);
   tempStr="odeparms(&N, parms);\n}\n/* Derivatives and 1 output variable */\n void derivs (int *neq, double *t, double *y, double *ydot,\n double *yout, int *ip)\n{";
   write(tempStr,file="mymod.c",sep="",append=TRUE);
   write(outputAttr,file="mymod.c",sep="",append=TRUE);
   write(strODE,file="mymod.c",sep="",append=TRUE);
   write("}",file="mymod.c",sep="",append=TRUE);
}

