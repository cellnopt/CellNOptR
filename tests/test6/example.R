# authors: TC
# This example provide an example and function to call a C function
# from CNORinterface shared library.


# You can test using :
# R --no-restore --no-save  < example.R
#
# or within a R shell by copy/paste the code that follows


# First read the MIDAS and SIF file to get some data


library(CellNOptR)
#setwd("/net/nas20/rcloud_data/storage/wdir/davidh/temp/test6")
#install.packages("CNORode_1.0.zip",repos=NULL);
library("CNORode")
#source("logic_based_ode_parameters_estimation_SSm_cluster.R");
#source("simulate_and_plot_ode_fitness.R")
s = readSif('model_espelin.sif')
m = readMIDAS('initialData_espelin.csv')

#s = readSif('model.sif')
#m = readMIDAS('initialData.csv')
cnolist = makeCNOlist(m, subfield=FALSE)
indices <- indexFinder(cnolist, s, verbose = TRUE)

#results=logic_based_ode_parameters_estimation_SSm(cnolist,s)
#logic_based_ode_MINLP_SSm(cnolist,s,ndiverse=10,dims_refset=6)

adjMat=incidence2Adjacency(s);
ode_parameters=makeParameterList(adjMat,s$namesSpecies,LB_k=0.05,UB_k=0.95,LB_n=1,UB_n=5,
LB_tau=0.0005,UB_tau=0.3,random=TRUE,default_tau=0.001,default_n=3);
#ode_parameters=makeParameterList(adjMat,s$namesSpecies,LB_k=0.05,UB_k=0.95,LB_n=1,UB_n=5,
#    LB_tau=10,UB_tau=0.1,random=TRUE,default_tau=0.001,default_n=3);
#print(system.time(get_logic_based_ode_data_simulation(cnolist,s,ode_parameters,indices,transfer_function=3)));
#results=list();
#results=logic_based_ode_parameters_estimation_SSm_cluster(cnolist,s,ode_parameters=ode_parameters,
#ndiverse=15,dim_refset=6,n_nodes=3,maxtime=1,maxeval=Inf,n_iter=1,maxNumSteps=500, transfer_function=2);
#save("results",file="results")
#load("results");
#xb=Inf;
#indexb=-1;
#for(i in 1:length(results))
#{
#    if(results[[i]]$fbest<xb)xb=results[[i]]$fbest;
#    index_b=i;
#}

#ode_parameters=makeParameterList(adjMat,s$namesSpecies,LB_k=0.05,UB_k=0.95,LB_n=1,UB_n=5,
#LB_tau=0.1,UB_tau=10,random=FALSE,default_tau=1,default_n=3);
#ode_parameters$parValues[ode_parameters$index_opt_pars]=results$xbest[[1]];
#simulate_and_plot_ode_fitness(cnolist,s,ode_parameters,transfer_function=3)#,plot_index_signals=c(1,2,3,4,5,6,7,8),
#plot_index_experiments=1:20)
obj_func=get_logic_ode_continuous_objective_function(cnolist,s,ode_parameters,indices,maxNumSteps=500);
#print(obj_func(ode_parameters$parValues))
for(i in 1:10)
{
    ode_parameters=makeParameterList(adjMat,s$namesSpecies,LB_k=0.05,UB_k=0.95,LB_n=1,UB_n=5,
    LB_tau=0.0005,UB_tau=0.3,random=TRUE,default_tau=0.001,default_n=3);
    print(obj_func(ode_parameters$parValues))
}
