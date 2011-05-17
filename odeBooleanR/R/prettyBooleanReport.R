prettyBooleanReport <-
function(boolModel)
{ 
  library(xtable);
  numOutputs=length(boolModel) 
  startString=readChar("latexHeader.txt",nchar=10000);
  repFile="boolReport/report.tex";
  write(startString,file=repFile,sep="",append=FALSE)
  for(i in 1:numOutputs)
  {
    numInputs=length(boolModel[[i]]$inputs);
    ttable=createInputTable(numInputs);
    ttable=cbind(ttable,boolModel[[i]]$truthTable)
    namesTruthTable=c(boolModel[[i]]$inputs,paste(boolModel[[i]]$output,"output",sep="-"));
    colnames(ttable)=namesTruthTable;
    ttable[which(ttable)]=as.integer(1);
    ttable[which(!ttable)]=as.integer(0);
    ttable=xtable(ttable);
    print(ttable, type="latex", file=repFile, append=TRUE,tabular.environment="longtable",floating=FALSE);
    write("\\clearpage",file=repFile,sep="\n",append=TRUE)
  }
  write("\\end{document}",file=repFile,sep="\n",append=TRUE)   
}

