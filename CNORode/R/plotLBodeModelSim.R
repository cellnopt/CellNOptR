plotLBodeModelSim <-function
(
		cnolist,				    model,					    ode_parameters=NULL,
		indices=NULL,			    adjMatrix=NULL,			  	time=1,
		verbose=0, 				    transfer_function=3,		reltol=1e-4,
		atol=1e-3,				    maxStepSize=Inf,		  	maxNumSteps=100000,
		maxErrTestsFails=50,  		large=FALSE,          		nsplit=4
)
{

	if(is.null(indices))indices=indexFinder(cnolist,model);
	if(is.null(adjMatrix))adjMatrix=incidence2Adjacency(model);
	if(is.null(ode_parameters))ode_parameters=createLBodeContPars(adjMatrix,model$namesSpecies);
  states_index=which(as.logical(getStates(adjMatrix)));

	sim_data=getLbodeModelSim(cnolist,model,
			ode_parameters,indices,time,verbose,transfer_function,
			reltol,atol,maxStepSize,maxNumSteps,maxErrTestsFails);

	sim_data=lapply(sim_data,function(x) x[,states_index]);

	times=cnolist$timeSignals;
	cnolist$valueSignals=sim_data;
	cnolist$namesSignals=model$namesSpecies[states_index];
	cnolist$namesCues=c(cnolist$namesStimuli,cnolist$namesInhibitors);

	cnolist$valueCues=cbind(cnolist$valueStimuli,cnolist$valueInhibitors);
	cnolist$valueCues=as.matrix(cnolist$valueCues);
	cnolist$valueCues[which(cnolist$valueCues>0)]=1;
	names(cnolist$valueCues)=cnolist$namesCues;

	if(large)
	{
		plotCNOlistLarge(cnolist,nsplit);
	}
	else
	{
		plotCNOlist(cnolist);
	}
	

  return(sim_data);
}