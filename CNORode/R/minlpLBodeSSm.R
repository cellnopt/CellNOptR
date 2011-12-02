minlpLBodeSSm <-
function
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
	if(is.null(ode_parameters))
	{
		ode_parameters=createLBodeContPars(adjMat,model$namesSpecies,random=TRUE);
	}
	
	n_cont=length(ode_parameters$index_opt_pars);
	n_int=dim(model$interMat)[2];
	problem=list();
	problem$f<-getLBodeMINLPObjFunction(cnolist,model,ode_parameters,indices,time,
			verbose,transfer_function,reltol,atol,maxStepSize,maxNumSteps,maxErrTestsFails);
	problem$x_L <- c(ode_parameters$LB[ode_parameters$index_opt_pars],matrix(0,1,n_int));
	problem$x_U <- c(ode_parameters$UB[ode_parameters$index_opt_pars],matrix(1,1,n_int));
	problem$x_0<- c(ode_parameters$parValues[ode_parameters$index_opt_pars],as.integer(round(runif(12))));
	problem$int_var=0;
	problem$bin_var =n_int;
	opts=list();
	opts$maxeval=maxeval;
	opts$maxtime=maxtime;
	if(!is.null(ndiverse))opts$ndiverse=ndiverse;      
	if(!is.null(dim_refset))opts$dim_refset=dim_refset;
	return(essR(problem,opts));	
}

