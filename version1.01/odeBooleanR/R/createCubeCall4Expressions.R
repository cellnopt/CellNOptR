createCubeCall4Expressions <-
function(inputNames,output,truthTable)
{
    numInputs=length(inputNames)
    numStates=2^numInputs;
    firstFlag=TRUE;
    strODE="(";
    if(length(truthTable)>0)
    {
      for(i in 1:numStates)
      {
        if(truthTable[i])
        {
          if(firstFlag){firstFlag=FALSE;}
          else strODE=paste(strODE,"+",sep="");
          tempBin=dec2bin(i-1,numInputs)
          for(j in 1:numInputs)
          {  
            ij=tempBin[j];
            if(!as.logical(ij))
            {
              tempStr=genNormHillCubes4Expressions(inputNames[j],output);
              strODE=paste(strODE,"(1-",tempStr,")",sep="");
            }
            else if(as.logical(ij))
            {
              tempStr=genNormHillCubes4Expressions(inputNames[j],output);
              strODE=paste(strODE,tempStr,sep="");
            }
            if(j<numInputs) strODE=paste(strODE,"*",sep="");
          }
        }      
      }
    }
    else{strODE='0'}; 
    parTau=paste(output,"_Tau",sep="");
    parInh=paste(output,"_Inh",sep=""); 
    strODE=paste("(",strODE,"-",output,")*",parTau,")","*(1-",parInh,")",sep="");
    return(strODE);
}

