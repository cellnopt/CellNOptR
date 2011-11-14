# authors: TC
# This example provide an example and function to call a C function
# from CNORinterface shared library.


# You can test using :
# R --no-restore --no-save  < example.R
#
# or within a R shell by copy/paste the code that follows


# First read the MIDAS and SIF file to get some data


library(CellNOptR)
#install.packages("CNORode_1.0.zip",repos=NULL);
#setwd("tests/test1");
library("CNORode");
#m = readMIDAS('initialData.csv')
#cnolist = makeCNOlist(m, subfield=FALSE)

load("CNOlistToyFB.RData");
cnolist=CNOlistToyFB;


model = readSif('ToyModelFeedbackDataGenerator.sif')
indices <- indexFinder(cnolist, model)

#indices<- findNONC(model, indices, verbose = TRUE)
#model <- cutNONC(model, indices)
#indices <- indexFinder(cnolist, model)
#model <- compressModel(model, indices)
#indices<- indexFinder(cnolist, model)


#results=logic_based_ode_parameters_estimation_eSSm(cnolist,s)
#logic_based_ode_MINLP_SSm(cnolist,s,ndiverse=10,dim_refset=6)

adjMat=incidence2Adjacency(model);
ode_parameters=makeParameterList(adjMat,model$namesSpecies,default_n=3,random=TRUE);
simulate_and_plot_ode_fitness(cnolist,model,ode_parameters,indices,transfer_function=3);
res=logic_based_ode_parameters_estimation_SSm(cnolist,model,ode_parameters,dim_refset=20,maxtime=1000);
ode_parameters$parValues[ode_parameters$index_opt_pars]=res$xbest;
dev.new();
simulate_and_plot_ode_fitness(cnolist,model,ode_parameters,indices,transfer_function=3);
                                           