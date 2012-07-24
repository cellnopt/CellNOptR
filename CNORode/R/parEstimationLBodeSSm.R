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
   tryCatch({library(MEIGOR)}, error=function(e){print("MEIGOR (essR) package not found.
	SSm not available. Install the package and load it or try the Genetic Algorithm
	optimiser instead.");return(ode_parameters);});


	adjMat=incidence2Adjacency(model);
	if(is.null(ode_parameters)){
		ode_parameters=createLBodeContPars(model,random=TRUE);
	}
	if(is.null(indices))indices <- indexFinder(cnolist,model,verbose=FALSE);
	
	#Check if essR is installed
	dummy_f<-function(x){
		return(0);
	}
	problem<-list(f=dummy_f,x_L=rep(0),x_U=c(1));
	opts<-list();
	opts$maxeval=0;
	opts$maxtime=0;

    val=essR(problem,opts)

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

