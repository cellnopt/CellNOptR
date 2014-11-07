#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2013 - EBI
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id: minlpLBodeSSm.R 3184 2013-01-21 13:50:31Z cokelaer $
minlpLBodeSSm <-
function
(
    cnolist,                model,                    ode_parameters=NULL,
    int_x0=NULL,            indices=NULL,            maxeval=Inf,
    maxtime=100,            ndiverse=NULL,            dim_refset=NULL,
    local_solver=NULL,      time=1,                    verbose=0,
    transfer_function=3,    reltol=1e-4,            atol=1e-3,
    maxStepSize=Inf,        maxNumSteps=100000,        maxErrTestsFails=50,
    nan_fac=1
)
{

    if (class(cnolist)=="CNOlist"){cnolist = compatCNOlist(cnolist)}
    adjMat=incidence2Adjacency(model);
    if(is.null(ode_parameters))
    {
        ode_parameters=createLBodeContPars(model,random=TRUE);

    }

    ode_parameters$model=model;
    #Check if essR is installed
    dummy_f<-function(x){
        return(0);
    }
    problem<-list(f=dummy_f,x_L=rep(0),x_U=c(1));
    opts<-list();
    opts$maxeval=0;
    opts$maxtime=0;

    val=tryCatch({essR(problem,opts)}, error=function(e){print("essR package not found.
    SSm not available. Install the package and load it or try the Genetic Algorithm
    optimiser instead.");return(ode_parameters);});
    ######################################################################################

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
    optimization_res=essR(problem,opts);

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

