#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv2 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-2.0.html
#
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/softare/cno
#
##############################################################################
# $Id$
minlpLBodeSSm <-
function
(
	cnolist,				model,					ode_parameters=NULL,
	int_x0=NULL,			indices=NULL,			maxeval=Inf,			
	maxtime=100,			ndiverse=NULL,			dim_refset=NULL, 		
	local_solver=NULL,      time=1,					verbose=0, 				
	transfer_function=3,	reltol=1e-4,			atol=1e-3,				
	maxStepSize=Inf,		maxNumSteps=100000,		maxErrTestsFails=50,	
	nan_fac=1
)
{	

    library(essR)
	adjMat=incidence2Adjacency(model);
	if(is.null(ode_parameters))
	{
		ode_parameters=createLBodeContPars(model,random=TRUE);
	}
	
	n_cont=length(ode_parameters$index_opt_pars);
	n_int=dim(model$interMat)[2];
	
	if(is.null(int_x0)){
		int_x0=as.integer(round(runif(n_int)));
	}
	
	problem=list();
	problem$f<-getLBodeMINLPObjFunction(cnolist,model,ode_parameters,indices,time,
			verbose,transfer_function,reltol,atol,maxStepSize,maxNumSteps,maxErrTestsFails);
	problem$x_L <- c(ode_parameters$LB[ode_parameters$index_opt_pars],matrix(0,1,n_int));
	problem$x_U <- c(ode_parameters$UB[ode_parameters$index_opt_pars],matrix(1,1,n_int));
	problem$x_0<- c(ode_parameters$parValues[ode_parameters$index_opt_pars],int_x0);
	problem$int_var=0;

	problem$bin_var=n_int;
	opts=list();
	opts$maxeval=maxeval;
	opts$maxtime=maxtime;
	if(!is.null(ndiverse))opts$ndiverse=ndiverse;      
	if(!is.null(dim_refset))opts$dim_refset=dim_refset;
	optimization_res=essR_optim(problem,opts);
	
	ode_parameters$bitString=optimization_res$xbest[(n_cont+1):(n_cont+n_int)];
	index_reactions=which(as.logical(ode_parameters$bitString));
	
	model$reacID=model$reacID[index_reactions];
	model$interMat=model$interMat[,index_reactions];
	model$notMat=model$notMat[,index_reactions];
	
	ode_parameters$model=model;
	
	ode_parameters$parValues[ode_parameters$index_opt_pars]=optimization_res$xbest[1:n_cont];
	ode_parameters$ssm_results=optimization_res;
	return(ode_parameters);
	
}

