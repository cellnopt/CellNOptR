plotLBodeFitness <-function
(
		cnolist,				    model,					    ode_parameters=NULL,
		indices=NULL,			    adjMatrix=NULL,			 	time=1,
		verbose=0, 				    transfer_function=3,		reltol=1e-4,
		atol=1e-3,				    maxStepSize=Inf,		 	maxNumSteps=100000,
		maxErrTestsFails=50,		plot_index_signals=NULL,	plot_index_experiments=NULL,
		plot_index_cues=NULL
)
{
	if(is.null(plot_index_experiments))plot_index_experiments=1:dim(cnolist$valueCues)[1];
	if(is.null(plot_index_cues))plot_index_cues=1:dim(cnolist$valueCues)[2];
    if(is.null(plot_index_signals))plot_index_signals=1:dim(cnolist$valueSignals[[1]])[2];
	if(is.null(indices))indices=indexFinder(cnolist,model);
	if(is.null(adjMatrix))adjMatrix=incidence2Adjacency(model);
	if(is.null(ode_parameters))ode_parameters=createLBodeContPars(model);

	sim_data=getLBodeDataSim(cnolist,model,
			ode_parameters,indices,time,verbose,transfer_function,
			reltol,atol,maxStepSize,maxNumSteps,maxErrTestsFails);

	times=cnolist$timeSignals;
	sim_data=lapply(sim_data,function(x)x[plot_index_experiments,plot_index_signals]);
	expResults=lapply(cnolist$valueSignals,function(x)x[plot_index_experiments,plot_index_signals]);

	namesSignals=cnolist$namesSignals[plot_index_signals];
	namesCues=c(cnolist$namesStimuli,cnolist$namesInhibitors);

	valueCues=cbind(cnolist$valueStimuli,cnolist$valueInhibitors);
	valueCues=as.matrix(valueCues);
	valueCues[which(valueCues>0)]=1;
	valueCues=valueCues[plot_index_experiments,plot_index_cues];
	names(valueCues)=namesCues[plot_index_cues];

	plotOptimResults(SimResults=sim_data,expResults=expResults,
			times=times,namesCues=namesCues,namesSignals=namesSignals,
			valueCues=valueCues);
			
  return(sim_data);
}
