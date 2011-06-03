odefy <-
function(model)
{
   fileStr="";
   numOutputs=length(model);
   fileStr=c();
   for(i in 1:numOutputs)
   {
      inputVec=model[[i]]$inputs;
      output=model[[i]]$output;
      truthTable=model[[i]]$truthTable;
      fileStr[i]=createCubeCall(inputVec,output,truthTable);
   }
   return(fileStr);
}

