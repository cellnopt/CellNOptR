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
# $Id: parEstimationLBode.R 3184 2013-01-21 13:50:31Z cokelaer $


parEstimationLBode<-function (cnolist, model, method="ga",
    ode_parameters = NULL, indices = NULL, paramsGA=NULL, paramsSSm=NULL)
{


    if (class(cnolist)=="CNOlist"){cnolist = compatCNOlist(cnolist)}
    if (method == "essm"){
        if (is.null(paramsSSm)){
            paramsSSm = defaultParametersSSm()
        }
        ode_parameters = parEstimationLBodeSSm(cnolist, model,
            ode_parameters=ode_parameters, indices=indices,
            maxeval=paramsSSm$maxeval,
            maxtime=paramsSSm$maxtime,
            ndiverse=paramsSSm$ndiverse,
            dim_refset=paramsSSm$dim_refset,
            local_solver=paramsSSm$local_solver,
            time=paramsSSm$time,
            verbose=paramsSSm$verbose,
            transfer_function=paramsSSm$transfer_function,
            reltol=paramsSSm$reltol,
            atol=paramsSSm$atol,
            maxStepSize=paramsSSm$maxStepSize,
            maxNumSteps=paramsSSm$maxNumSteps,
            maxErrTestsFails=paramsSSm$maxErrTestsFails,
            nan_fac=paramsSSm$nan_fac)


    }
    else if(method=="ga"){
        if (is.null(paramsGA)){
            paramsGA = defaultParametersGA()
        }
        ode_parameters = parEstimationLBodeGA(cnolist, model,
            ode_parameters=ode_parameters,
            indices=indices,
            mutationChance=paramsGA$mutationChance,
            popSize=paramsGA$popSize,
            iters=paramsGA$iters,
            elitism=paramsGA$elitism,
            time=paramsGA$time,
            monitor=paramsGA$monitor,
            verbose=paramsGA$verbose,
            transfer_function=paramsGA$transfer_function,
            reltol=paramsGA$reltol,
            atol=paramsGA$atol,
            maxStepSize=paramsGA$maxStepSize,
            maxNumSteps=paramsGA$maxNumSteps,
            maxErrTestsFails=paramsGA$maxErrTestsFails,
            nan_fac=paramsGA$nan_fac)
    }
    else{
        stop ("method argument must be either 'ga' or 'essm'." )
    }
    return(ode_parameters)
}

