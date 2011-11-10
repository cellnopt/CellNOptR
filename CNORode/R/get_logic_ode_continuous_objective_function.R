get_logic_ode_continuous_objective_function <-
function
(
		cnolist,				model,					ode_parameters,	
		indices,				time=1,					verbose=0, 
		transfer_function=3,	reltol=1e-4,			atol=1e-3,
		maxStepSize=Inf,		maxNumSteps=100000,		maxErrTestsFails=50
)
{
	adjMatrix=incidence2Adjacency(model);
	sim_function=get_simulation_function(cnolist,model,adjMatrix1=adjMatrix,
			indices1=indices, odeParameters1=ode_parameters$parValues, time1=time,verbose1=verbose,
			transfer_function1=transfer_function,reltol1=reltol,atol1=atol,maxStepSize1=maxStepSize,
			maxNumSteps1=maxNumSteps,maxErrTestsFails1=maxErrTestsFails)

	logic_ode_continuous_objective_function<-function
	(
			x,							cnolist1=cnolist,		model1=model,
			adjMatrix1=adjMatrix,		indices1=indices,		ode_parameters1=ode_parameters,
			sim_function1=sim_function
	)
	{
		ode_parameters1$parValues[ode_parameters1$index_opt_pars]=x;
		sim=sim_function1(cnolist1,model1,ode_parameters1$parValues);
		sim<-unlist(lapply(sim,function(x) x[,indices1$signals]));
		measured_values=unlist(cnolist1$valueSignals);                                                    
		NaNs=which(is.na(sim));
		not_NaNs=which(!is.na(sim));
		error=sum((sim[not_NaNs]-measured_values[not_NaNs])^2);
		if(is.na(error))error=0;
		res=(error+length(NaNs))/length(sim);
		return(res);
	}
	return(logic_ode_continuous_objective_function);
}

