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
# $Id$
getLBodeContObjFunction<-
function
(
        cnolist,                model,                    ode_parameters,
        indices=NULL,           time=1,                   verbose=0,
        transfer_function=3,    reltol=1e-4,            atol=1e-3,
        maxStepSize=Inf,        maxNumSteps=100000,        maxErrTestsFails=50,
        nan_fac=1, useVariances=F, initial_state=0.1
)
{

    if (class(cnolist)=="CNOlist"){cnolist = compatCNOlist(cnolist)}
    adjMatrix=incidence2Adjacency(model);

    if(is.null(indices))indices=indexFinder(cnolist,model,verbose=FALSE);

    sim_function=getLBodeSimFunction(cnolist,model,adjMatrix1=adjMatrix,
            indices1=indices, odeParameters1=ode_parameters$parValues, time1=time,verbose1=verbose,
            transfer_function1=transfer_function,reltol1=reltol,atol1=atol,maxStepSize1=maxStepSize,
            maxNumSteps1=maxNumSteps,maxErrTestsFails1=maxErrTestsFails,
            initial_state1=initial_state)

    logic_ode_continuous_objective_function<-function
    (
            x,                            cnolist1=cnolist,        model1=model,
            adjMatrix1=adjMatrix,        indices1=indices,        ode_parameters1=ode_parameters,
            sim_function1=sim_function,        nan_fac1=nan_fac,
            useVariances1=useVariances, initial_state1=initial_state
    )
    {
        ode_parameters1$parValues[ode_parameters1$index_opt_pars]=x;
        sim=sim_function1(cnolist1,model1,ode_parameters1$parValues);
        sim<-as.vector(unlist(lapply(sim,function(x) x[,indices1$signals])));
        measured_values<-as.vector(unlist(lapply(cnolist1$valueSignals,function(x)x)));

        if (useVariances == T){
            variances<-as.vector(unlist(lapply(cnolist1$valueVariances,function(x)x)));
            variances = sapply(variances, function(x){min(x, 0.5)}) # make sure variance are <=0.5
            weights = 1. - sqrt(variances*2)
        }else{
            weights = 1. - 0. * measured_values
        }

        NaNs_sim = which(is.na(sim));
        not_NaNs_data = which(!is.na(measured_values));
        not_NaNs_sim = which(!is.na(sim));
        not_NaNs = intersect(not_NaNs_sim,not_NaNs_data);
        NaNs = intersect(not_NaNs_data,NaNs_sim);

        error = sum(  weights[not_NaNs] *  (sim[not_NaNs]-measured_values[not_NaNs])^2);

        norm_weights = sum(weights[not_NaNs])/length(not_NaNs)
        if (is.na(norm_weights)==F){
            error = error/norm_weights
        }
        res=(error+length(NaNs)*nan_fac1)/length(not_NaNs_data);
        return(res);
    }
    return(logic_ode_continuous_objective_function);
}

