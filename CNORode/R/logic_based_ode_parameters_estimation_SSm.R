logic_based_ode_parameters_estimation_SSm <-function
(
		cnolist,			model,			ode_parameters=NULL,
		maxeval=Inf,		maxtime=100,	ndiverse=NULL,
		dim_refset=NULL
)
{
	library(eSSmR)
	library(Rsolnp)
	
	if(is.null(ode_parameters))
	{
		adjMat=incidence2Adjacency(model);
		ode_parameters=makeParameterList(adjMat,model$namesSpecies,random=TRUE);
	}
	
	problem=list();
	problem$f<-get_logic_ode_continuous_objective_function(cnolist,model,ode_parameters,indices);
	problem$x_L <- ode_parameters$LB[ode_parameters$index_opt_pars];
	problem$x_U <- ode_parameters$UB[ode_parameters$index_opt_pars];
	problem$x_0<- ode_parameters$parValues[ode_parameters$index_opt_pars];
	problem$int_var =0;
	problem$bin_var =0;
	opts=list();
	opts$maxeval=maxeval;
	opts$maxtime=maxtime;
	if(!is.null(ndiverse))opts$ndiverse=ndiverse;      
	if(!is.null(dim_refset))opts$dim_refset=dim_refset;  
	opts$local_solver="SOLNP";
	return(essR(problem,opts));	
}

