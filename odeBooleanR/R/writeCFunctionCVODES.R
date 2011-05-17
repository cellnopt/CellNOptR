writeCFunctionCVODES <-
function(inputs,outputs,strODE,numOutputs,parVecForPrint)
{   
   headerDir=paste(.find.package("odeBooleanR"),"/textHeaders/startStringForFile.txt",sep="");
   startString=readChar(headerDir,nchar=10000);
   outputsStr=paste(outputs,sep=",",collapse=",");
   outputsStr=paste("double ",outputsStr,";",sep="");
   
   declareParVecForPrint=paste(parVecForPrint,sep="",collapse=",");
   declareParVecForPrint=paste("double ",parVecForPrint,";",sep="");
   
   outputAttr=c();

   for(i in 1:numOutputs)
   {
      outputAttr[i]=paste(outputs[i],"=Ith(y,",i,");",sep="");  
   }
   
   for(i in 1:numOutputs)
   {
      strODE[i]=paste("Ith(ydot,",i,")=",strODE[i],";",sep="");        
   }
   
   sizeParVec=length(parVecForPrint);
   count=0;
   
   for(i in 1:length(parVecForPrint))
   {
     parVecForPrint[i]=paste(parVecForPrint[i],"=ptr[",i-1,"];",sep=""); 
   }
   
   cFileDir=paste(.find.package("odeBooleanR"),"/models/myODE.c",sep="");  
   write(startString,file=cFileDir,sep="",append=FALSE)
   write(outputsStr,file=cFileDir,sep="",append=TRUE);
   write(declareParVecForPrint,file=cFileDir,sep="",append=TRUE);
   write(parVecForPrint,file=cFileDir,sep="",append=TRUE);
   write(outputAttr,file=cFileDir,sep="",append=TRUE);
   write(strODE, file=cFileDir, sep="", append = TRUE);  
   write("return(0);\n}",file=cFileDir, sep="", append = TRUE); 
}

