# authors: TC
# This example provide an example and function to call a C function
# from CNORinterface shared library.


# You can test using :
# R --no-restore --no-save  < example.R
#
# or within a R shell by copy/paste the code that follows


# First read the MIDAS and SIF file to get some data
get_logic_based_ode_model_simulation<-function
(
		cnolist,				model,					ode_parameters=NULL,	
		indices=NULL,			time=1,					verbose=0, 
		transfer_function=3,	reltol=1e-4,			atol=1e-3,	
		maxStepSize=Inf,		maxNumSteps=100000,		maxErrTestsFails=50
)
{
	adjMat=incidence2Adjacency(model);
z
	if(is.null(indices))indices <- indexFinder(cnolist,model,verbose=FALSE);
	if(is.null(ode_parameters))ode_parameters=makeParameterList(adjMat,model$namesSpecies);
	sim_function=get_simulation_function(cnolist,model,adjMat,
			indices1=indices, odeParameters1=ode_parameters1$parValues, time1=time,verbose1=verbose,
			transfer_function1=transfer_function,reltol1=reltol,atol1=atol,maxStepSize1=maxStepSize,
			maxNumSteps1=maxNumSteps,maxErrTestsFails1=maxErrTestsFails);
	return(sim_function(cnolist,model,ode_parameters$parValues));
}

get_logic_based_ode_data_simulation<-function
(
		cnolist,				model,					ode_parameters=NULL,	
		indices=NULL,			time=1,					verbose=0, 
		transfer_function=3,	reltol=1e-4,			atol=1e-3,	
		maxStepSize=Inf,		maxNumSteps=100000,		maxErrTestsFails=50
)
{
	adjMat=incidence2Adjacency(model);
	if(is.null(indices))indices <- indexFinder(cnolist,model,verbose=FALSE);
	if(is.null(ode_parameters))ode_parameters=makeParameterList(adjMat,model$namesSpecies);
	sim_function=get_simulation_function(cnolist,model,adjMat,
			indices1=indices, odeParameters1=ode_parameters1$parValues, time1=time,verbose1=verbose,
			transfer_function1=transfer_function,reltol1=reltol,atol1=atol,maxStepSize1=maxStepSize,
			maxNumSteps1=maxNumSteps,maxErrTestsFails1=maxErrTestsFails);
	sim=sim_function(cnolist,model,ode_parameters$parValues);
	return(lapply(sim,function(x) x[,indices$signals]));
}

plot_fit_ode_simulation<-function
(
		cnolist,				model,					ode_parameters=NULL,	
		indices=NULL,			time=1,					verbose=0, 
		transfer_function=3,	reltol=1e-4,			atol=1e-3,	
		maxStepSize=Inf,		maxNumSteps=100000,		maxErrTestsFails=50
)
{
	sim_data=get_logic_based_ode_data_simulation(cnolist,model,
			ode_parameters,indices,time,verbose,transfer_function,	
			reltol,atol=,maxStepSize=Inf,maxNumSteps,maxErrTestsFails=50);
	
	times=cnolist$timeSignals;
	expResults=cnolist$valueSignals;  
	namesSignals=cnolist$namesSignals;
	namesCues=c(cnolist$namesStimuli,cnolist$namesInhibitors);
	
	valueCues=cbind(cnolist$valueStimuli,cnolist$valueInhibitors)
	valueCues=as.matrix(valueCues);
	valueCues[which(valueCues>0)]=1;
	names(valueCues)=namesCues;
	
	plotOptimResults(SimResults=sim_data,expResults=expResults,
			times=times,namesCues=namesCues,namesSignals=namesSignals,
			valueCues=valueCues); 
	
}

logic_based_ode_continous_PSO <-function
(
		cnolist,				model,			ode_parameters=NULL,
		maxeval=Inf,			maxtime=100,	swarm_size=10,	
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

#library(CellNOptR)
install.packages("CNORode_1.0.zip",repos=NULL);
#source("psoptim.R")
#setwd("tests/test3");
#library("CNORode")

#s = readSif('model.sif')
#m = readMIDAS('initialData.csv')
#cnolist = makeCNOlist(m, subfield=FALSE)

#indices <- indexFinder(cnolist, s, verbose = TRUE)
#modelNCNOindices <- findNONC(s, indices, verbose = TRUE)
#s <- cutNONC(s, modelNCNOindices);

#results=logic_based_ode_parameters_estimation_SSm(cnolist,s,ndiverse=10,dim_refset=6)
#logic_based_ode_MINLP_SSm(cnolist,s,ndiverse=10,dim_refset=6)

#adjMat=incidence2Adjacency(s);
#ode_parameters=makeParameterList(adjMat,s$namesSpecies);
#simulator<-get_simulation_function(cnolist,s,adjMat,indices,ode_parameters,reltol=1e-6,atol=1e-6);

#sim=simulator(cnolist,s,ode_parameters$parValues)
#value_signals<-lapply(sim,function(x) x[,indices$signals]);

#print(unlist(sum((unlist(cnolist$valueSignals)-unlist(value_signals))^2)));
#plotCNOlistLargePDF(cnolist,"original.pdf",nsplit=10)
#plotCNOlist(cnolist);
#cnolist$valueSignals=value_signals;
#windows();
#plotCNOlist(cnolist);

#res=get_logic_based_ode_data_simulation(cnolist,s)
#res=get_logic_based_ode_model_simulation(cnolist,s)
#plot_fit_ode_simulation(cnolist,s);
#plotCNOlistLargePDF(cnolist,"simulated.pdf",nsplit=10)
#print(system.time(logic_based_ode_continous_PSO(cnolist,s,maxeval=500)))
