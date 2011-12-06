parEstimationLBodeSSm <-function
(
		cnolist,				model,					ode_parameters=NULL,
		indices=NULL,			maxeval=Inf,			maxtime=100,			
		ndiverse=NULL,			dim_refset=NULL, 		local_solver=NULL,      
		time=1,					verbose=0, 				transfer_function=3,	
		reltol=1e-4,			atol=1e-3,				maxStepSize=Inf,		
		maxNumSteps=100000,		maxErrTestsFails=50,	nan_fac=1
)
{
	adjMat=incidence2Adjacency(model);
	if(is.null(ode_parameters)){
		ode_parameters=createLBodeContPars(model,random=TRUE);
	}
	if(is.null(indices))indices <- indexFinder(cnolist,model,verbose=FALSE);
	
	problem=list();
	problem$f<-getLBodeContObjFunction(cnolist,	model,ode_parameters,indices,
	time,verbose,transfer_function,reltol,atol,maxStepSize,maxNumSteps,maxErrTestsFails);
	problem$x_L <- ode_parameters$LB[ode_parameters$index_opt_pars];
	problem$x_U <- ode_parameters$UB[ode_parameters$index_opt_pars];
	problem$x_0<- ode_parameters$parValues[ode_parameters$index_opt_pars];
	problem$int_var =0;
	problem$bin_var =0;
	opts=list();
	opts$maxeval=maxeval;
	opts$maxtime=maxtime;
	if(!is.null(local_solver))opts$local_solver=local_solver;
	if(!is.null(ndiverse))opts$ndiverse=ndiverse;      
	if(!is.null(dim_refset))opts$dim_refset=dim_refset;  
	results=essR(problem,opts);
	ode_parameters$parValues[ode_parameters$index_opt_pars]=results$xbest;
	ode_parameters$ssm_results=results;
	return(ode_parameters);	
}

