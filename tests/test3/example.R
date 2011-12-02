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
setwd("C:/Users/David/Desktop/stuff/CNOR_ode/tests/test3");
library("CNORode")

model = readSif('model.sif')
data = readMIDAS('initialData.csv')
cnolist = makeCNOlist(data, subfield=FALSE)

indices <- indexFinder(cnolist, model, verbose = TRUE)

#results=logic_based_ode_parameters_estimation_SSm(cnolist,s)
#logic_based_ode_MINLP_SSm(cnolist,s,ndiverse=10,dim_refset=6)

adjMat=incidence2Adjacency(s);
ode_parameters=makeParameterList(adjMat,s$namesSpecies,default_n=3);
simulate_and_plot_ode_fitness(cnolist,s,ode_parameters,indices,transfer_function=3);
