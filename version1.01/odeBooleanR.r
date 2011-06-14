dec2bin<-function(x, ndigits)
{
        b=2;
        xi <- as.integer(x)
        if(any(is.na(xi) | ((x-xi)!=0)))
                print(list(ERROR="x not integer", x=x))
        N <- length(x)
        xMax <- max(x)
        #ndigits <- (floor(logb(xMax, base=2))+1)
        Base.b <- array(NA, dim=c(N, ndigits))
        for(i in 1:ndigits){#i <- 1
                Base.b[, ndigits-i+1] <- (x %% b)
                x <- (x %/% b)
        }
        if(N ==1) Base.b[1, ] else Base.b
}

bin2dec<-function(binaryvector)
{
  sum(2^(which(rev(binaryvector)==TRUE)-1))
}

genNormHillCubes<-function(input,output)
{
  k=paste(input,"_k_",output,sep="");
  n=paste(input,"_n_",output,sep="");
  strODE=paste("pow(",input,",",n,")/(pow(",input,",",n,")+pow(",k,",",n,"))*(1+pow(",k,",",n,"))",sep="");
  res=list(ODE=str)
  return(strODE);
}

getObjectiveFunction <-function(paramInfo1,y01,expData1,times1,obsIndex1,fullpars1,indexPars1,cFileName)
{
  if(cFileName=="mymod")
  {
    objectiveFunction<-function(x, y0=y01,paramInfo=paramInfo1,expData=expData1,times=times1,obsIndex=obsIndex1,fullpars=fullpars1,indexPars=indexPars1)
    {    
       fullpars[indexPars]=x;
       simData<-lsode(y0, times, func = "derivs", parms=fullpars, jacfunc=NULL, dllname="mymod",initfunc ="initmod",nout = 0);
       res=(expData[,obsIndex]-simData[,obsIndex])^2;
      
       return(res);
    }             
  }
  else if(cFileName=="myODE")
  {   
    objectiveFunction<-function(x, y0=y01,paramInfo=paramInfo1,expData=expData1,times=times1,obsIndex=obsIndex1,fullpars=fullpars1,indexPars=indexPars1)
    {    
       fullpars[indexPars]=x;
       simData<-cvodes(y0,times,"myODE","rhs",fndata =fullpars,verbose =FALSE);
       res=sum((expData[,obsIndex]-simData[,obsIndex])^2); 
       if(is.nan(res))res=10000000000;
       print(res);
       return(res);
    }
  }
  else 
  {  
    objectiveFunction<-function(x)
    {
      print("NO VALID OBJECTIVE FUNCTION WAS PROVIDED")
      return(NA);
    }
    print("A valid File Name was not provided");
  }
  
 return(objectiveFunction);   
}

initializeParameters <-function(paramInfo,typeInit)
{
  if(typeInit=="randomNK")
  {
    indexNK=c(paramInfo$index_k,paramInfo$index_n)
    paramInfo$paramValues[indexNK]=runif(length(indexNK));
  }
  else if(typeInit=="randomBoundedNK")
  {
    indexK=paramInfo$index_k;
    indexN=paramInfo$index_n;
    paramInfo$paramValues[indexK]=runif(length(indexK),0.5,5);
    paramInfo$paramValues[indexN]=runif(length(indexN),0.5,7);
    paramInfo$UB[indexK]=5;
    paramInfo$LB[indexK]=0.5;
    paramInfo$UB[indexN]=7;
    paramInfo$LB[indexN]=0.5;
  }
  else
  {
    print("Check you have selected any valid Option");
  }
  return(paramInfo);
}

compileAndLoad <-function(cFileName)
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

createCubeCall <-function(inputNames,output,truthTable)
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
              tempStr=genNormHillCubes(inputNames[j],output);
              strODE=paste(strODE,"(1-",tempStr,")",sep="");
            }
            else if(as.logical(ij))
            {
              tempStr=genNormHillCubes(inputNames[j],output);
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

createInputTable <-function(numInputs)
{ 
  mat=c();
  for(i in 1:2^numInputs)mat=rbind(mat,as.logical(dec2bin(i-1,numInputs)));
  return(mat);
}

createInSilicoData<-function(paramInfo,y0,times,cFileName)
{
   if(cFileName=="mymod") 
   {
      out<-lsode(y0, times, func = "derivs", parms=paramInfo$paramValues, jacfunc=NULL, dllname="mymod",initfunc ="initmod",nout = 0);
   }
   else if(cFileName=="myODE")
   {
     out<-cvodes(y0,times,"myODE","rhs",fndata=paramInfo$paramValues,verbose=FALSE,maxnumsteps=1000);
   }
   else
   {       
      out=NA;
      print("A valid File Name was not provided");
   }
   return(out);
}

createOdefyModel <-
function(intermat,notmat,specID)
{
  res=list();
  rules=transitionRules(intermat,notmat);
  count=0;
  for(i in 1:length(specID))
  { 
    tempRules=rules[[i]]$mat; 
     
    numRules=length(tempRules[1,])
    indexInputs=as.matrix(which(tempRules!=0,arr.ind=TRUE));
    indexInputs=unique(indexInputs[,1]);
    numInputs=length(indexInputs);
    nameInputs=specID[indexInputs];
    if(numRules>0)
    { 
      inputTable=createInputTable(numInputs);
      count=count+1
      res[count]=list(inputs=NULL,output=NULL,truthTable=NULL)
      truthTable=logical();
      truthTable[seq(1,2^numInputs)]=FALSE; 
      for(j in 1:numRules)
      {  
        rule=tempRules[,j];
        truth=which(rule!=0)
        rule[which(rule==-1)]=0;
        rule=as.logical(rule);
        stateVecIndex=c();
        
        for(k in 1:length(truth))stateVecIndex[k]=which(indexInputs==truth[k]);
        for(k in 1:2^numInputs)
        {
          if(!any(which((inputTable[k,stateVecIndex]==rule[truth])==FALSE)))
          {
            truthTable[k]=TRUE;
          }
        } 
      }
      res[[count]]$inputs=nameInputs;
      res[[count]]$output=specID[i];
      res[[count]]$truthTable=truthTable;
    }
  }
    return(res)  
}
  
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

optParamatersDEoptim <-
function(pars,ub,lb,f)
{  
   library(DEoptim);
   DEoptim(fn =f, lower =lb , upper =ub, control = list(NP = 10, itermax = 10, trace = TRUE));
}

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

writeCfunction<-function(model)
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

writeCFunctionCVODES<-function(inputs,outputs,strODE,numOutputs,parVecForPrint)
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

writeCFuntionDeSolve<-function(inputs,outputs,strODE,numOutputs,parVecForPrint)
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

indexCNO2ODE<-function(CNOlist,odePars)
{ 
  CNOindexSignals=c();
  ODEindexSignals=c();
  for(i in 1:length(odePars$outputs))
  {
    tempIndex=which(CNOlist$namesSignals==odePars$outputs[i]);
    CNOindexSignals=c(CNOindexSignals,tempIndex);
    if(length(tempIndex)>0)ODEindexSignals=c(ODEindexSignals,i);
  }
  
  CNOindexStimuli=c();
  ODEindexStimuli=c();
  for(i in 1:length(odePars$inputs))
  {
    tempIndex=which(CNOlist$namesStimuli==odePars$inputs[i]);
    CNOindexStimuli=c(CNOindexStimuli,tempIndex);
    if(length(tempIndex)>0) ODEindexStimuli=c(ODEindexStimuli,i);
  }
  
  CNOindexInhibitors=c();
  ODEindexInhibitors=c();
  for(i in 1:length(odePars$outputs))
  {
    tempIndex=which(CNOlist$namesInhibitors==odePars$outputs[i]);
    CNOindexInhibitors=c(CNOindexInhibitors,tempIndex);
    if(length(tempIndex)>0)ODEindexInhibitors=c(ODEindexInhibitors,i);
  }
  res=list(CNOindexSignals=CNOindexSignals,ODEindexSignals=ODEindexSignals,
    CNOindexStimuli=CNOindexStimuli,ODEindexStimuli=ODEindexStimuli,
      CNOindexInhibitors=CNOindexInhibitors,ODEindexInhibitors=ODEindexInhibitors);
  
  return(res);
}

simulateData<-function(CNOlist,indexes,odePars)
{
  numExperiments=length(CNOlist$valueStimuli[,1]);
  times=CNOlist$timeSignals;
  res=list();
  for(i in 1:numExperiments)
  {
    odePars$paramValues[odePars$index_inputs[indexes$ODEindexStimuli]]=
      CNOlist$valueStimuli[i,indexes$CNOindexStimuli];
    odePars$paramValues[odePars$index_inh[indexes$ODEindexInhibitors]]=
      CNOlist$valueInhibitors[i,indexes$CNOindexInhibitors];
    y0=matrix(0,1,length(odePars$outputs));
    y0[indexes$ODEindexSignals]=CNOlist$valueSignals[[1]][i,indexes$CNOindexSignals]
      CNOlist$valueInhibitors[i,indexes$CNOindexInhibitors];
    tempRes=createInSilicoData(odePars,y0,times,"myODE");
    res[[i]]=tempRes[,indexes$ODEindexSignals];    
  }
  
  return(res);
}

simulateModel<-function(CNOlist,indexes,odePars,times)
{
  numExperiments=length(CNOlist$valueStimuli[,1]);
  res=list();
  for(i in 1:numExperiments)
  {
    odePars$paramValues[odePars$index_inputs[indexes$ODEindexStimuli]]=
      CNOlist$valueStimuli[i,indexes$CNOindexStimuli];
    odePars$paramValues[odePars$index_inh[indexes$ODEindexInhibitors]]=
      CNOlist$valueInhibitors[i,indexes$CNOindexInhibitors];
      y0=matrix(0,1,length(odePars$outputs));
    y0[indexes$ODEindexSignals]=CNOlist$valueSignals[[1]][i,indexes$CNOindexSignals]
      CNOlist$valueInhibitors[i,indexes$CNOindexInhibitors];
    tempRes=createInSilicoData(odePars,y0,times,"myODE");
    res[[i]]=tempRes;    
  }
  
  return(res);
}

convertTimeSeriesToList<-function(sim,CNOlist)
{
  nTimePoints=length(data.frame(sim[1])[,1]);
  nSpecies=length(data.frame(sim[1])[1,])
  nExperiments=length(sim);
  res=list();
  for(i in 1:nTimePoints)res[[i]]=data.frame();
  for(i in 1:nExperiments)
  {
    temp=data.frame(sim[i]);
    for(j in 1:nTimePoints)
    {
      res[[j]]=rbind(res[[j]],temp[j,]);
    } 
  }
  for(i in 1:nTimePoints)colnames(res[[i]])=CNOlist$namesSignals;
  return(res)  
}

plotODESimFitness<-function(SimResults,CNOList)
{ 
  times=CNOlist$timeSignals;
  expResults=CNOlist$valueSignals;  
  SimResults=convertTimeSeriesToList(SimResults,CNOlist);
  
  namesSignals=CNOlist$namesSignals;
  
  namesCues=c(CNOlist$namesStimuli,CNOlist$namesInhibitors);
  
  valueCues=cbind(CNOlist$valueStimuli,CNOlist$valueInhibitors)
  valueCues=as.matrix(valueCues);
  valueCues[which(valueCues>0)]=1;
  names(valueCues)=namesCues;
  
  plotOptimResults(SimResults=SimResults,expResults=expResults,times=times,namesCues=namesCues,namesSignals=namesSignals,valueCues=valueCues); 
}

getODEModelExpressions<-function(model)
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
   

   strODE=odefy4Expressions(model);
   
   expressions=c();
   for(i in 1:length(strODE))   
   {
    expressions[i]=parse(text=strODE[i]);
   }
   
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
   numOutputs=numOutputs,outputs=outputs,numInputs=countInputs,inputs=inputs,expressions=expressions);
   
   return(paramInfo)
}

odefy4Expressions<-function(model)
{
   fileStr="";
   numOutputs=length(model);
   fileStr=c();
   for(i in 1:numOutputs)
   {
      inputVec=model[[i]]$inputs;
      output=model[[i]]$output;
      truthTable=model[[i]]$truthTable;
      fileStr[i]=createCubeCall4Expressions(inputVec,output,truthTable);
   }
   return(fileStr);
}

genNormHillCubes4Expressions<-function(input,output)
{
  k=paste(input,"_k_",output,sep="");
  n=paste(input,"_n_",output,sep="");
  strODE=paste(input,"^",n,"/(",input,"^",n,"+",k,"^",n,")*(1+",k,"^",n,")",sep="");
  res=list(ODE=str)
  return(strODE);
}

createCubeCall4Expressions<-function(inputNames,output,truthTable)
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

getPartialDerivs<-function(mathModel)
{
  numParameters=length(mathModel$index_k)+length(mathModel$index_n)+
      length(mathModel$index_k);
  numStates=length(mathModel$expressions);
  listPDeriv=list();
  count=0;
  indexNonZeroPDeriv=c();
  tag=c(); 
  
  for(i in 1:numStates)
  { 
    for(j in 1:numParameters)
    {
      count=count+1;
      tag=c(tag,paste('d',mathModel$outputs[i],'/d',mathModel$parNames[j],sep=""));
      listPDeriv[[count]]=D(mathModel$expressions[i],mathModel$parNames[j]);
      if(listPDeriv[[count]]!=0)indexNonZeroPDeriv=c(indexNonZeroPDeriv,count);
    }
  }
  res=list(listPDeriv=listPDeriv,tag=tag,indexNonZeroPDeriv=indexNonZeroPDeriv);
  return(res)
}

odeSensitivityAnalisys<-function(CNOlist,CNOModel,odePars,times)
{ 
  #compileResult=compileAndLoad("myODE");
  
  indexes=indexCNO2ODE(CNOlist,odePars);
  
  TTModel=createOdefyModel(CNOModel$interMat,CNOModel$notMat,CNOModel$namesSpecies);
  
  expressionsInfo=getODEModelExpressions(TTModel)
  
  pderivs=getPartialDerivs(expressionsInfo)
  
  res=writeSensitivityCode(expressionsInfo, pderivs)
  
  return(res)
  #res=computeSensitivities(expressionsInfo,pderivs,CNOlist,indexes,times)
}

writeSensitivityCode<-function(expressionsInfo, pderivs)
{
  exps=expressionsInfo$expressions;
  numExps=length(exps);
  res=c();                         
  ode=c();
  for(i in 1:numExps){ode[i]=parseExpression(exps[i]);}
  loc=regexpr.simple2('\\w+\\^\\w+',str)
  expr=ode[1];
  while(loc!=-1)
  {   
     tempSub=substring(expr,loc$from,loc$to);
     tempSub=strsplit(tempSub,"\\^");
     insertSub=paste("pow(",tempSub[[1]][1],",",tempSub[[1]][2],")",sep="",collapse="")
     print(insertSub)
     print(tempSub)
     exp1=substr(expr,1,loc$from-1)
     exp2=substr(expr,loc$to+1,nchar(expr));
     expr=paste(exp1,exp2,sep="");
     loc=regexpr.simple2('\\w+\\^\\w+',expr);
     
     
  }
  return(res);
}

parseExpression<-function(expr)
{
    expr=paste(deparse(expr,width.cutoff=100000),collapse="",sep="");
    expr=sub('expression\\(\\(', '', expr, perl = TRUE); ## Perl-style white space
    expr=sub('\\)$', '',expr, perl = TRUE);
    expr=gsub(' ','',expr,perl=TRUE);  
}

simulateODEAndPlotFitness<-function(CNOlist,CNOModel,odePars)
{
  indexes=indexCNO2ODE(CNOlist,odePars);
  
  simData=simulateData(CNOlist,indexes,odePars)
   
  plotODESimFitness(simData,CNOlist); 
}

createODEModel<-function(CNOlist,CNOModel)
{
  TTModel=createOdefyModel(CNOModel$interMat,CNOModel$notMat,CNOModel$namesSpecies);
                                    
  paramInfo=writeCfunction(TTModel);
  
  #compileResult=compileAndLoad("myODE");
  
  
  return(paramInfo);
}

regexpr.simple2<-function(a,x)
{
	# Use as follows: >regexpr.simple2("[:][0-9.]+[,]",x)
	# Returns the first and last positions in string "x" of the first match with string "a"
	result<-list()
	w<-regexpr(a,x)
	result$from<-w[1]
	if(w[1]>-1)
  {
	   result$to<-w[1]+attr(w, "match.length")-1
	}
	return(result)
}




library(CellNOptR);
library(odeBooleanR)
setwd("C:/Users/David/Desktop/testR/");

load("CNOlistSilico")

CNOlist=CNOlistSilico;

model=readSif(sifFile = "modNet2.sif")

odePars=createODEModel(CNOlist,model);

indexes=indexCNO2ODE(CNOlist,odePars)

#simulateODEAndPlotFitness(CNOlist,CNOModel,odePars)

res=odeSensitivityAnalisys(CNOlist,model,odePars,times);
#es=sensitivityMatrix(odePars,CNOlist)


