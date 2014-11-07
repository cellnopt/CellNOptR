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
# $Id: getLBodeModelSim.R 3881 2013-08-06 10:52:07Z cokelaer $
getLBodeModelSim<-function
(
        cnolist,                model,                    ode_parameters=NULL,
        indices=NULL,            timeSignals=NULL,        time=1,
        verbose=0,                transfer_function=3,    reltol=1e-4,
        atol=1e-3,                maxStepSize=Inf,        maxNumSteps=100000,
        maxErrTestsFails=50
)
{

    # this is a temporary solution to deal with CNOlist class
    if (class(cnolist)=="CNOlist"){cnolist = compatCNOlist(cnolist)}


    adjMat=incidence2Adjacency(model);
    if(is.null(indices))indices <- indexFinder(cnolist,model,verbose=FALSE);
    if(is.null(ode_parameters))ode_parameters=createLBodeContPars(model);
    if(!is.null(timeSignals))cnolist$timeSignals=timeSignals;

    # TC aug 2013. it looks like cnolist and model parameters are useless here ?? could set to anything ??
    # smae for ode_parameters$parValues
    # need to be fixed at some point

    
    #sim_function=getLBodeSimFunction(cnolist,model,adjMat,
    #        indices, ode_parameters$parValues, time,verbose,
    #        transfer_function,reltol,atol,maxStepSize,
    #        maxNumSteps,maxErrTestsFails);

    sim_function=getLBodeSimFunction(NULL,NULL,adjMat,
            indices, NULL, time,verbose,
            transfer_function,reltol,atol,maxStepSize,
            maxNumSteps,maxErrTestsFails);

    # TC aug 2013 !!! only parValues are used in ode_parameters 
    sim_data=sim_function(cnolist,model,ode_parameters$parValues)

    sim_data=lapply(sim_data, function(x) as.matrix(x))
    return(sim_data);
}
