compileAndLoad <-
function(cFileName)
{
  try(dyn.unload(paste(cFileName, .Platform$dynlib.ext, sep = "")));
  compileResult=system(paste('R CMD SHLIB ',cFileName,".c",sep=""));
  if(as.logical(compileResult))print("Compilation Failed");
  dyn.load(paste(cFileName, .Platform$dynlib.ext, sep = ""));   
  if(cFileName=="mymod") {library(deSolve);}
  else if(cFileName=="myODE"){library(Rsundials);}
  else {print("A valid File Name was not provided");}
  return(compileResult);
}

