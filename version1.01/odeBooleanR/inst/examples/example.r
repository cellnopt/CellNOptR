library(CNelephanteR);
library(odeBooleanR)

setwd(paste(.find.package("odeBooleanR"),"/examples/",sep=""));

load("CNOlistSilico")

CNOlist=CNOlistSilico;

model=readSif(sifFile = "modNet2.sif")

odePars=createODEModel(CNOlist,model);

#indexes=indexCNO2ODE(CNOlist,odePars)

simulateODEAndPlotFitness(CNOlist,CNOModel,odePars)

#res=odeSensitivityAnalisys(CNOlist,model,odePars,times);