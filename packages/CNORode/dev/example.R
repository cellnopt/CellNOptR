# authors: TC
# This example provide an example and function to call a C function
# from CNORinterface shared library.


# You can test using :
# R --no-restore --no-save  < example.R
#
# or within a R shell by copy/paste the code that follows


# First read the MIDAS and SIF file to get some data

makeParameterList<-
		function
				(
				adjMat,namesSpecies,
				LB_n=1,LB_k=0.4,LB_tau=1e-1,
				UB_n=3,UB_k=0.6,UB_tau=10,
				default_n=3,default_k=0.5,default_tau=1,
				LB_in=c(),UB_in=c(),opt_n=TRUE,
				opt_k=TRUE,opt_tau=TRUE,random=FALSE
)
{
	isState=getStates(adjMat);
	numStates=as.integer(sqrt(length(adjMat)));
	parNames=c();
	parValues=c();
	index_k=c();
	index_n=c();
	index_tau=c();
	LB=c();
	UB=c();
	index_opt_pars=c();
	
	count=0;
	for(j in 1:numStates)
	{ 
		if(isState[j])
		{
			jCol=adjMat[,j];
			inputs=which(as.logical(jCol));
			numInputs=length(inputs);
			for(i in 1:numInputs)
			{
				count=count+1;
				parNames[count]=paste(namesSpecies[inputs[i]],"_n_",namesSpecies[j],sep="");
				parValues[count]=default_n;
				index_k=c(index_k,count);
				LB[count]=LB_n;
				UB[count]=UB_n;
				if(opt_n)index_opt_pars=c(index_opt_pars,count);
				
				count=count+1;
				parNames[count]=paste(namesSpecies[inputs[i]],"_k_",namesSpecies[j],sep="");
				parValues[count]=default_k;
				index_n=c(index_k,count);
				LB[count]=LB_k;
				UB[count]=UB_k;
				if(opt_k)index_opt_pars=c(index_opt_pars,count);
			}
			count=count+1;
			parNames[count]=paste("tau_",namesSpecies[j],sep="");
			parValues[count]=default_tau;
			index_tau=c(index_tau,count);
			LB[count]=LB_tau;
			UB[count]=UB_tau;
			if(opt_tau)index_opt_pars=c(index_opt_pars,count);
		}
	}
	if(length(LB_in)==length(LB))LB=LB_in;
	if(length(UB_in)==length(UB))UB=UB_in;
	if(random)parValues=LB[index_opt_pars]+((UB[index_opt_pars]-LB[index_opt_pars])*runif(length(index_opt_pars)));
	
	parList=list(parNames=parNames,parValues=parValues,
			index_opt_pars=index_opt_pars,index_n=index_n,
			index_k=index_k,index_tau=index_tau,LB=LB,UB=UB);
	
	return(parList);
}

getStates=function(adjacency)
{
	nSpecies=dim(adjacency)[1];
	count=0;
	res=matrix(0,1,nSpecies);
	for(j in 1:nSpecies)
	{
		for(i in 1:nSpecies)
		{
			if(adjacency[i,j])
			{
				res[j]=1;
			}
		}
	}
	return(res);
}

incidence2Adjacency=function(model)
{ 
	incidence=model$interMat;
	nNodes=dim(incidence)[1];
	nEdges=dim(incidence)[2];
	adjacency=matrix(0,nNodes,nNodes);
	
	for(j in 1:nEdges)
	{
		for(i in 1:nNodes)
		{
			if(incidence[i,j]==1)
			{
				node1=i;
				for(k in 1:nNodes)
				{
					if(incidence[k,j]==-1)
					{
						adjacency[k,node1]=1
					}
				}
			}
		}
	}
	return(adjacency);  
}

make_default_ode_pars<-function(interMat)
{
	adjMat=getAdjMat(interMat);
	n_nodes=dim(adjMat)[1];
	ode_parameters=c();
	
	for (j in 1:n_nodes) 
	{
		isState=FALSE;
		for (i in 1:n_nodes) 
		{
			if(adjMat[i,j])
			{
				ode_parameters=c(ode_parameters,3);
				ode_parameters=c(ode_parameters,0.5);
				isState=TRUE;
			}	
		}
		if(isState)ode_parameters=c(ode_parameters,1);
	}
	return(ode_parameters)
}



get_logic_ode_continuous_objective_function<-function(cnolist1,model1,ode_parameters1,indices1,
		time1=1,verbose1=0, transfer_function1=3,reltol1=1e-4,atol1=1e-3,maxStepSize1=Inf,maxNumSteps1=100000,
		maxErrTestsFails1=50)
{
	adjMatrix1=incidence2Adjacency(model1);
	sim_function1=get_simulation_function(cnolist1,model1,adjMatrix=adjMatrix1,
			indices=indices1, odeParameters=ode_parameters1$parValues, time=time1,verbose=verbose1,
			transfer_function=transfer_function1,reltol=reltol1,atol=atol1,maxStepSize=maxStepSize1,
			maxNumSteps=maxNumSteps1,maxErrTestsFails=maxErrTestsFails1)

	logic_ode_continuous_objective_function<-
			function(x,cnolist=cnolist1,model=model1,adjMatrix=adjMatrix1,
					indices=indices1,ode_parameters=ode_parameters1,sim_function=sim_function1)
	{
		ode_parameters$parValues[ode_parameters$index_opt_pars]=x;
		sim=sim_function(cnolist,model,indices,ode_parameters$parValues);
		sim<-unlist(lapply(sim,function(x) x[,indices$signals]));
		measured_values=unlist(cnolist$valueSignals);                                                    
		NaNs=which(is.na(sim));
		not_NaNs=which(!is.na(sim));
		error=sum((sim[not_NaNs]-measured_values[not_NaNs])^2);
		if(is.na(error))error=0;
		res=(error+length(NaNs))/length(sim);
		return(res);
	}
	return(logic_ode_continuous_objective_function);
}


get_simulation_function<-function(cnolist1,model1,adjMatrix1,
		indices1, odeParameters1,time1=1,verbose1=0,transfer_function1=3,
		reltol1=1e-4,atol1=1e-3,maxStepSize1=Inf,maxNumSteps1=100000,maxErrTestsFails1=50)
{
				
	simulate_logic_based_ode_model <- 
	function(cnolist, sif, odeParameters,indices=indices1, adjMatrix=adjMatrix1, time=time1,verbose=verbose1, 
	transfer_function=transfer_function1,reltol=reltol1,atol=atol1,maxStepSize=maxStepSize1,
	maxNumSteps=maxNumSteps1,maxErrTestsFails=maxErrTestsFails1)
	{	 
		interMat <- as.integer(as.vector(t(sif$interMat)))
		notMat <- as.integer(as.vector(t(sif$notMat)))
		adjMatrix <- as.integer(as.vector(t(adjMatrix)))
		nRows <- as.integer(dim(sif$interMat)[1])
		
		nCols <- as.integer(dim(sif$interMat)[2])
		verbose=as.integer(verbose);
		# ode
	
		nPars <- as.integer(length(odeParameters))
	
		# cnolist
		timeSignals <- as.double(cnolist$timeSignals)
		valueInhibitors <- as.double(t(cnolist$valueInhibitors))
		valueSignals <- as.double(t(cnolist$valueSignals[[time]]))
		valueStimuli <- as.double(t(cnolist$valueStimuli))
		nTimes = as.integer(length(cnolist$timeSignals))
		
		# [[1]] allows to access to the first object in the list and retrieve its dimensions
		nExperiments = as.integer(dim(cnolist$valueSignals[[time]]))
		
		#indices
		nSignals <- as.integer(length(indices$signals))
		
		indexSignals <- as.integer(as.vector(indices$signals))
		nStimuli <- as.integer(length(indices$stimulated))
		indexStimuli <- as.integer(as.vector(indices$stimulated))
		nInhibitors <- as.integer(length(indices$inhibited))
		indexInhibitors <- as.integer(as.vector(indices$inhibited))
		transfer_function=as.integer(transfer_function);
		reltol=as.double(reltol);
		atol=as.double(atol);
		maxStepSize=as.double(maxStepSize);
		maxNumSteps=as.integer(maxNumSteps);
		maxErrTestsFails=as.integer(maxErrTestsFails);
		res = .Call("sim_logic_ode",interMat,notMat,adjMatrix,nRows,nCols,nPars,timeSignals,
				valueInhibitors,valueSignals,valueStimuli,nTimes,nExperiments,nSignals,
				indexSignals,nStimuli,indexStimuli,nInhibitors,indexInhibitors,odeParameters,
				verbose,transfer_function,reltol,atol,maxStepSize,maxNumSteps,maxErrTestsFails);
		return(res);
	}
	return(simulate_logic_based_ode_model);
}

get_logic_ode_MINLP_objective_function<-function(cnolist1,model1,ode_parameters1,indices1,
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
		sim=sim_function(cnolist,temp_model,indices,ode_parameters$parValues);
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

logic_based_ode_parameters_estimation_SSm<-
function(cnolist,model,ode_parameters=NULL,maxeval=Inf,maxtime=100,ndiverse=NULL,dim_refset=NULL)
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

logic_based_ode_MINLP_SSm<-
		function(cnolist,model,ode_parameters=NULL,maxeval=Inf,maxtime=100,ndiverse=NULL,dim_refset=NULL)
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
	problem$f<-get_logic_ode_MINLP_objective_function(cnolist,model,ode_parameters,indices);
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



#library(CellNOptR)
##detach("package:testPackage");

#setwd("tests/test4");
#install.packages("CNORode_1.0.zip",repos=NULL);
#setwd("tests/test1");
#library("CNORode")

#s = readSif('model.sif')
#m = readMIDAS('initialData.csv')
#cnolist = makeCNOlist(m, subfield=FALSE)

#indices <- indexFinder(cnolist, s, verbose = TRUE)
#modelNCNOindices <- findNONC(s, indices, verbose = TRUE)
#s <- cutNONC(s, modelNCNOindices);

#results=logic_based_ode_parameters_estimation_SSm(cnolist,s)
#logic_based_ode_MINLP_SSm(cnolist,s,ndiverse=10,dim_refset=6)



