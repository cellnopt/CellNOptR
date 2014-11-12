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



simulate <-function(cnolist,model,ode_parameters=NULL,
        indices=NULL, adjMatrix=NULL, time=1, verbose=0, transfer_function=3, 
        reltol=1e-4, atol=1e-3, maxStepSize=Inf, maxNumSteps=100000,
        maxErrTestsFails=50)
{
    if (class(cnolist)=="CNOlist"){cnolist = compatCNOlist(cnolist)}
    if(is.null(indices))
        indices=indexFinder(cnolist,model);
    if(is.null(adjMatrix))
        adjMatrix=incidence2Adjacency(model);
    if(is.null(ode_parameters))
        ode_parameters=createLBodeContPars(model);

    timeSignals=NULL;

    sim_data=getLBodeDataSim(cnolist,model,
            ode_parameters,indices,timeSignals,time,verbose,
            transfer_function,reltol,atol,maxStepSize,maxNumSteps,
            maxErrTestsFails);



    return(sim_data)
}
