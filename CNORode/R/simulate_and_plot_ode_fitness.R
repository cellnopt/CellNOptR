simulate_and_plot_ode_fitness<-function
(
		cnolist,				      model,					      ode_parameters=NULL,
		indices=NULL,			    adjMatrix=NULL,			  time=1,
		verbose=0, 				    transfer_function=3,	reltol=1e-4,
		atol=1e-3,				    maxStepSize=Inf,		  maxNumSteps=100000,
		maxErrTestsFails=50
)
{

	if(is.null(indices))indices=indexFinder(cnolist,model);
	if(is.null(adjMatrix))adjMatrix=incidence2Adjacency(model);
	if(is.null(ode_parameters))ode_parameters=makeParameterList(adjMatrix,model$namesSpecies);

	sim_data=get_logic_based_ode_data_simulation(cnolist,model,
			ode_parameters,indices,time,verbose,transfer_function,
			reltol,atol,maxStepSize,maxNumSteps,maxErrTestsFails);

	times=cnolist$timeSignals;
	expResults=cnolist$valueSignals;
	namesSignals=cnolist$namesSignals;
	namesCues=c(cnolist$namesStimuli,cnolist$namesInhibitors);

	valueCues=cbind(cnolist$valueStimuli,cnolist$valueInhibitors);
	valueCues=as.matrix(valueCues);
	valueCues[which(valueCues>0)]=1;
	names(valueCues)=namesCues;

	plotOptimResults(SimResults=sim_data,expResults=expResults,
			times=times,namesCues=namesCues,namesSignals=namesSignals,
			valueCues=valueCues);
			
  return(sim_data);
}