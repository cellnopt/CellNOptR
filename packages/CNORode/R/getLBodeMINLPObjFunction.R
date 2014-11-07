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
# $Id: getLBodeMINLPObjFunction.R 3184 2013-01-21 13:50:31Z cokelaer $
getLBodeMINLPObjFunction <-function
(
        cnolist,                model,                    ode_parameters,
        indices=NULL,            time=1,                    verbose=0,
        transfer_function=3,    reltol=1e-4,            atol=1e-3,
        maxStepSize=Inf,        maxNumSteps=100000,        maxErrTestsFails=50,
        nan_fac=1
)
{
    if (class(cnolist)=="CNOlist"){cnolist = compatCNOlist(cnolist)}

    adjMatrix=incidence2Adjacency(model);

    if(is.null(indices))indices=indexFinder(cnolist,model,verbose=FALSE);

    n_cont=length(ode_parameters$index_opt_pars);
    n_int=dim(model$interMat)[2];

    sim_function<-getLBodeSimFunction(cnolist,model,adjMatrix1=adjMatrix,
            indices1=indices,odeParameters1=ode_parameters$parValues, time1=time,verbose1=verbose,
            transfer_function1=transfer_function,reltol1=reltol,atol1=atol,maxStepSize1=maxStepSize,
            maxNumSteps1=maxNumSteps,maxErrTestsFails1=maxErrTestsFails)

    logic_ode_MINLP_objective_function<-function
    (
            x,                        cnolist1=cnolist,                    model1=model,
            adjMatrix1=adjMatrix,    n_cont1=n_cont,                        n_int1=n_int,
            indices1=indices,        ode_parameters1=ode_parameters,        sim_function1=sim_function,
            nan_fac1=nan_fac
    )
    {
        temp_model=model1;
        x_cont=as.double(x[seq(1,n_cont1)]);
        x_int=as.integer(x[seq((n_cont1+1),(n_cont1+n_int1))]);
        if(length(which(as.logical(x_int)))<2)return(1);
        temp_model$interMat=model1$interMat[,which(as.logical(x_int))];
        temp_model$notMat=model1$notMat[,which(as.logical(x_int))];
        ode_parameters1$parValues[ode_parameters1$index_opt_pars]=x_cont;
        sim=sim_function1(cnolist1,temp_model,ode_parameters1$parValues);
                sim<-as.vector(unlist(lapply(sim,function(x) x[,indices1$signals])));
        measured_values<-as.vector(unlist(lapply(cnolist1$valueSignals,function(x)x)));
        NaNs_sim=which(is.na(sim));
        not_NaNs_data=which(!is.na(measured_values));
        not_NaNs_sim=which(!is.na(sim));
        not_NaNs=intersect(not_NaNs_sim,not_NaNs_data);
        NaNs=intersect(not_NaNs_data,NaNs_sim);
        error=sum((sim[not_NaNs]-measured_values[not_NaNs])^2);
        res=(error+length(NaNs)*nan_fac1)/length(not_NaNs_data);
        return(res);
    }
    return(logic_ode_MINLP_objective_function);
}

