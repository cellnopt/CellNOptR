get_logic_ode_MINLP_objective_function <-
function(cnolist1,model1,ode_parameters1,indices1,
		time1=1,verbose1=0, transfer_function1=3,reltol1=1e-4,atol1=1e-3,maxStepSize1=Inf,
		maxNumSteps1=100000,maxErrTestsFails1=50)
{
	adjMatrix1=incidence2Adjacency(model1);
	n_cont1=length(ode_parameters1$index_opt_pars);
	n_int1=dim(model1$interMat)[2];
	sim_function1<-get_simulation_function(cnolist1,model1,adjMatrix=adjMatrix1,
			indices=indices1,odeParameters=ode_parameters1$parValues, time=time1,verbose=verbose1,
			transfer_function=transfer_function1,reltol=reltol1,atol=atol1,maxStepSize=maxStepSize1,
			maxNumSteps=maxNumSteps1,maxErrTestsFails=maxErrTestsFails1)
	
	logic_ode_MINLP_objective_function<-
			function(x,cnolist=cnolist1,model=model1,adjMatrix=adjMatrix1,n_cont=n_cont1,n_int=n_int1,
					indices=indices1,ode_parameters=ode_parameters1,sim_function=sim_function1)
	{
		temp_model=model;
		x_cont=as.double(x[seq(1,n_cont)]);
		x_int=as.integer(x[seq((n_cont+1),(n_cont+n_int))]);
		if(length(which(as.logical(x_int)))<2)return(1);
		temp_model$interMat=model$interMat[,which(as.logical(x_int))];
		temp_model$notMat=model$notMat[,which(as.logical(x_int))];
		ode_parameters$parValues[ode_parameters$index_opt_pars]=x_cont;
		sim=sim_function(cnolist,temp_model,ode_parameters$parValues);
		sim<-unlist(lapply(sim,function(x) x[,indices$signals]));
		measured_values=unlist(cnolist$valueSignals);                                                    
		NaNs=which(is.na(sim));
		not_NaNs=which(!is.na(sim));
		error=sum((sim[not_NaNs]-measured_values[not_NaNs])^2);
		if(is.na(error))error=0;
		res=(error+length(NaNs))/length(sim);
		return(res);
	}
	return(logic_ode_MINLP_objective_function);
}

