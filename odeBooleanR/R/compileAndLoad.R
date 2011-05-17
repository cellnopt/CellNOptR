compileAndLoad <-
function(cFileName)
{  
  actualDir=getwd();
  cFileDir=paste(.find.package("odeBooleanR"),"/models/",sep="")
  setwd(cFileDir);
  try(dyn.unload(paste(cFileName, .Platform$dynlib.ext, sep = "")));
  compileResult=system(paste('R CMD SHLIB ',cFileName,".c",sep=""));
  if(as.logical(compileResult))print("Compilation Failed");
  dyn.load(paste(cFileDir,cFileName, .Platform$dynlib.ext, sep = ""));   
  if(cFileName=="mymod") {library(deSolve);}
  else if(cFileName=="myODE"){library(Rsundials);}
  else {print("A valid File Name was not provided");}
  setwd(actualDir);
  return(compileResult);
}

