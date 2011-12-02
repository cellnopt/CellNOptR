get_logic_based_ode_data_simulation<-function
(
		cnolist,				    model,					ode_parameters=NULL,
		indices=NULL,			    time=1,					verbose=0,
		transfer_function=3,		reltol=1e-4,			atol=1e-3,
		maxStepSize=Inf,		  	maxNumSteps=100000,		maxErrTestsFails=50
)
{
	adjMat=incidence2Adjacency(model);
	if(is.null(indices))indices <- indexFinder(cnolist,model,verbose=FALSE);
	if(is.null(ode_parameters))ode_parameters=makeParameterList(adjMat,model$namesSpecies);
	sim_function=get_simulation_function(cnolist,model,adjMat,
			indices1=indices, odeParameters1=ode_parameters$parValues, time1=time,verbose1=verbose,
			transfer_function1=transfer_function,reltol1=reltol,atol1=atol,maxStepSize1=maxStepSize,
			maxNumSteps1=maxNumSteps,maxErrTestsFails1=maxErrTestsFails);
	sim=sim_function(cnolist,model,ode_parameters$parValues);
	sim=lapply(sim,function(x) x[,indices$signals]);
	return(sim);
}