#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2013 - EMBL-EBI
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
# $Id: $




parEstimationLBodeGA<-function (cnolist, model, ode_parameters = NULL, indices = NULL, 
    mutationChance=NA, popSize=200,iters=100,elitism=NA, time = 1,monitor=TRUE,
	verbose = 0, transfer_function = 3, reltol = 1e-04, atol = 0.001,
	maxStepSize = Inf, maxNumSteps = 1e+05, maxErrTestsFails = 50, nan_fac = 1,
    initial_state=0.1) 
{

    if (class(cnolist)=="CNOlist"){cnolist = compatCNOlist(cnolist)}
    adjMat = incidence2Adjacency(model)
    if (is.null(ode_parameters)) {
        ode_parameters = createLBodeContPars(model, random = TRUE)
    }
    if (is.null(indices)) 
        indices <- indexFinder(cnolist, model, verbose = FALSE)
    problem = list()
    f_obj <- getLBodeContObjFunction(cnolist, model, ode_parameters, 
        indices, time, verbose, transfer_function, reltol, atol, 
        maxStepSize, maxNumSteps, maxErrTestsFails, initial_state)
    x_L <- ode_parameters$LB[ode_parameters$index_opt_pars]
    x_U <- ode_parameters$UB[ode_parameters$index_opt_pars]
    x_0 <- ode_parameters$parValues[ode_parameters$index_opt_pars]
    
    if(monitor){
            plot.new();
            monitor_func <- function(obj){
                # plot the population
                xlim = c(0, obj$iters);
                par(mfrow=c(1,1));
                nIter=obj$iters-length(which(is.na(obj$best)));
                popSize=obj$popSize;
                heading = paste("Number of Evaluations=",popSize*nIter);
                plot(obj$best,main=heading,xlab="Number of Iterations", ylab="Objective Fuction",col='red')
                lines(obj$best, xlim=xlim,col="blue");
                print(paste("Iteration:",nIter,"  Best_f:",min(obj$best,na.rm=TRUE),"  N_Evals:",popSize*nIter,sep="")); 
            }
    }
    else{
        monitor_func=NULL;
    }

    res=rbga(x_L, x_U,popSize=popSize, iters=iters,monitorFunc=monitor_func, evalFunc=f_obj,
     showSettings=FALSE, verbose=FALSE,elitism=elitism,mutationChance=mutationChance);

    best_individual_index=which(res$evaluations==min(res$evaluations));
    best_individual=res$population[best_individual_index,];

    ode_parameters$parValues[ode_parameters$index_opt_pars] = best_individual;
    ode_parameters$res=res;
	
	
    return(ode_parameters)
}
