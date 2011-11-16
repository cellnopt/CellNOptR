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
library("CNORode")
s = readSif('model.sif')
m = readMIDAS('initialData.csv')
cnolist = makeCNOlist(m, subfield=FALSE)

indices <- indexFinder(cnolist, s, verbose = TRUE)

#results=logic_based_ode_parameters_estimation_SSm(cnolist,s)
#logic_based_ode_MINLP_SSm(cnolist,s,ndiverse=10,dim_refset=6)

adjMat=incidence2Adjacency(s);
ode_parameters=makeParameterList(adjMat,s$namesSpecies,LB_k=0.15,UB_k=0.85,LB_n=1,UB_n=5,
LB_tau=0.001,UB_tau=0.03,random=TRUE);
#print(system.time(get_logic_based_ode_data_simulation(cnolist,s,ode_parameters,indices,transfer_function=3)));
results=logic_based_ode_parameters_estimation_SSm(cnolist,s,ode_parameters=ode_parameters,ndiverse=10,dim_refset=6);

