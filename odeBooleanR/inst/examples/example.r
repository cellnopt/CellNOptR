library(CNelephanteR);
library(odeBooleanR);
#setwd('testBoolean');
boolModel=readSif('Network.sif');
odeModel=createOdefyModel(boolModel$interMat,boolModel$notMat,boolModel$namesSpecies);
paramInfo=writeCfunction(odeModel)
y0=c();
y0[seq(1:13)]=0.2;
times=seq(0,10,0.1);
cFileName="myODE";
compileResult=compileAndLoad(cFileName);
if(compileResult==1)stop();
inSilico=createInSilicoData(paramInfo,y0,times,cFileName)
inSilicoPlot=ts(inSilico[seq(1,101),seq(2,10)], frequency = 1)
plot(inSilicoPlot,nc=3,plot.type ="multi",col=c(seq(2,16)))