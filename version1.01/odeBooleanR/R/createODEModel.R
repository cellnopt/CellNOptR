createODEModel <-
function(CNOlist,CNOModel)
{
  TTModel=createOdefyModel(CNOModel$interMat,CNOModel$notMat,CNOModel$namesSpecies);
                                    
  paramInfo=writeCfunction(TTModel);
  
  compileResult=compileAndLoad("myODE");
  
  return(paramInfo);
}

