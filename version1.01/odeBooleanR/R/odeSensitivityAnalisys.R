odeSensitivityAnalisys <-
function(CNOlist,CNOModel,odePars,times)
{ 
  compileResult=compileAndLoad("myODE");
  
  indexes=indexCNO2ODE(CNOlist,odePars);
  
  TTModel=createOdefyModel(CNOModel$interMat,CNOModel$notMat,CNOModel$namesSpecies);
  
  expressionsInfo=getODEModelExpressions(TTModel)
  
  pderivs=getPartialDerivs(expressionsInfo)
  
  res=computeSensitivities(expressionsInfo,pderivs,CNOlist,indexes,times)
}

