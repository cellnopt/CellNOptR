# authors: TC
# This example provide an example and function to call a C function
# from CNORinterface shared library.


# You can test using :
# R --no-restore --no-save  < example.R
#
# or within a R shell by copy/paste the code that follows


# First read the MIDAS and SIF file to get some data
# authors: TC
# This example provide an example and function to call a C function
# from CNORinterface shared library.


# You can test using :
# R --no-restore --no-save  < example.R
#
# or within a R shell by copy/paste the code that follows


# First read the MIDAS and SIF file to get some data


logic_based_ode_continous_PSO <-function
(
		cnolist,				model,			ode_parameters=NULL,
		maxeval=Inf,			maxtime=100,	swarm_size=4,	
		exploration_w=c(1,0),	hybrid=FALSE
)
{	
	if(is.null(ode_parameters))
	{
		adjMat=incidence2Adjacency(model);
		ode_parameters=makeParameterList(adjMat,model$namesSpecies,random=TRUE);
	}
	
	f_obj<-get_logic_ode_continuous_objective_function(cnolist,
			model,ode_parameters,indices);
	x_L <- ode_parameters$LB[ode_parameters$index_opt_pars];
	x_U <- ode_parameters$UB[ode_parameters$index_opt_pars];
	x_0<- ode_parameters$parValues[ode_parameters$index_opt_pars];
	control=list();
	control$maxf=maxeval;
	control$abstol=0;
	control$s=swarm_size;
	control$trace=1;
	res=psoptim(x_0,f_obj,lower=x_L,upper=x_U,control=control)
	return(res);
	
}

library(CellNOptR)
setwd("tests/test4")
#install.packages("CNORode_1.0.tar.gz",type="source");
#source("psoptim.R")
#source("temp.R");
library("CNORode")

s = readSif('model.sif')
m = readMIDAS('initialData.csv')
cnolist = makeCNOlist(m, subfield=FALSE)

indices <- indexFinder(cnolist, s, verbose = TRUE)

#modelNCNOindices <- findNONC(s, indices, verbose = TRUE)
#s <- cutNONC(s, modelNCNOindices);

#results=logic_based_ode_parameters_estimation_SSm(cnolist,s,ndiverse=10,dim_refset=6)
#logic_based_ode_MINLP_SSm(cnolist,s,ndiverse=10,dim_refset=6)

adjMat=incidence2Adjacency(s);
ode_parameters=makeParameterList(adjMat,s$namesSpecies,default_n=3);
#simulate_and_plot_ode_fitness(cnolist,s,ode_parameters,indices,transfer_function=3);
simulate_and_plot_model(cnolist,s,ode_parameters,indices,transfer_function=3);