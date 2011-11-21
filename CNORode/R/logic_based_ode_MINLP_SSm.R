logic_based_ode_MINLP_SSm <-
function(cnolist,model,ode_parameters=NULL,maxeval=Inf,maxtime=100,ndiverse=NULL,dim_refset=NULL,nan_fac=1)
{
	library(eSSmR)
	library(Rsolnp)
	
	if(is.null(ode_parameters))
	{
		adjMat=incidence2Adjacency(model)
		ode_parameters=makeParameterList(adjMat,model$namesSpecies,random=TRUE);
	}
	
	n_cont=length(ode_parameters$index_opt_pars);
	n_int=dim(model$interMat)[2];
	problem=list();
	problem$f<-get_logic_ode_MINLP_objective_function(cnolist,model,ode_parameters,indices,nan_fac);
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

