writeCfunction <-
function(model)
{
   parKvec=c();
   parNvec=c();
   parTauVec=c();
   parInhVec=c();
   inputs=c();
   outputs=c();
   numOutputs=length(model);
   parStr="";
   countPars=0;
   countInputs=0;
   
   for(i in 1:numOutputs)
   {
      outputs[i]=model[[i]]$output;
   }
   
   for(i in 1:length(outputs))
   {
      inputVec=model[[i]]$inputs;
      for(j in 1:length(inputVec))
      { 
        countPars=countPars+1;
        if(is.na(which(inputs==inputVec[j])[1]) && is.na(which(outputs==inputVec[j])[1]))
        {
          countInputs=countInputs+1;
          inputs[countInputs]=inputVec[j];
        }
        parKvec[countPars]=paste(inputVec[j],"_k_",outputs[i],sep="");
        parNvec[countPars]=paste(inputVec[j],"_n_",outputs[i],sep="");
      }
      parTauVec[i]=paste(outputs[i],"_Tau",sep="");
      parInhVec[i]=paste(outputs[i],"_Inh",sep=""); 
   }
   

   strODE=odefy(model);
   
   parVecForPrint=append(parKvec,append(parNvec,append(parTauVec,parInhVec)));
     
   sizeParVec=length(parVecForPrint);
   count=0;
   
   for(i in (sizeParVec+1):(sizeParVec+length(inputs)))
   {
      count=count+1;
      parVecForPrint[i]=inputs[count];
   }

   #Parameters Index
   parNames=parVecForPrint;
   size=1;
   index_k=seq(size,length(parKvec));
   size=length(index_k);
   index_n=seq(size+1,size+length(parNvec));
   size=size+length(index_n);
   index_tau=seq(size+1,size+length(parTauVec));
   size=size+length(parTauVec);
   index_inh=seq(size+1,size+length(parInhVec));
   size=size+length(parInhVec);
   index_inputs=seq(size+1,size+countInputs);
   size=size+countInputs;
   paramValues=c();
   paramValues[index_k]=0.5;
   paramValues[index_n]=3;
   paramValues[index_tau]=1;
   paramValues[index_inh]=0;
   paramValues[index_inputs]=1;
   paramInfo=list(parNames=parNames,paramValues=paramValues,index_k=index_k,index_n=index_n,index_tau=index_tau,index_inh=index_inh,index_inputs=index_inputs,
   numOutputs=numOutputs,outputs=outputs,numInputs=countInputs,inputs=inputs);
        
   #writeCFuntionDeSolve(inputs,outputs,strODE,numOutputs,parVecForPrint)
   writeCFunctionCVODES(inputs,outputs,strODE,numOutputs,parVecForPrint);
   
   return(paramInfo)
}

