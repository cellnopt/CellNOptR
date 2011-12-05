getLBodeModelSim<-function
(
		cnolist,				model,					ode_parameters=NULL,
		indices=NULL,			time=1,					verbose=0,
		transfer_function=3,	reltol=1e-4,			atol=1e-3,
		maxStepSize=Inf,		maxNumSteps=100000,		maxErrTestsFails=50
)
{
	adjMat=incidence2Adjacency(model);
	if(is.null(indices))indices <- indexFinder(cnolist,model,verbose=FALSE);
	if(is.null(ode_parameters))ode_parameters=createLBodeContPars(adjMat,model$namesSpecies);
	sim_function=getLBodeSimFunction(cnolist,model,adjMat,
			indices, ode_parameters1$parValues, time,verbose,
			transfer_function,reltol,atol,maxStepSize,
			maxNumSteps,maxErrTestsFails);
	return(sim_function(cnolist,model,ode_parameters$parValues));
}