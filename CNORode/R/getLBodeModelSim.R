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
getLBodeModelSim<-function
(
		cnolist,				model,					ode_parameters=NULL,
		indices=NULL,			timeSignals=NULL,		time=1,					
		verbose=0,				transfer_function=3,	reltol=1e-4,			
		atol=1e-3,				maxStepSize=Inf,		maxNumSteps=100000,		
		maxErrTestsFails=50
)
{
	adjMat=incidence2Adjacency(model);
	if(is.null(indices))indices <- indexFinder(cnolist,model,verbose=FALSE);
	if(is.null(ode_parameters))ode_parameters=createLBodeContPars(model);
	if(!is.null(timeSignals))cnolist$timeSignals=timeSignals;
	sim_function=getLBodeSimFunction(cnolist,model,adjMat,
			indices, ode_parameters$parValues, time,verbose,
			transfer_function,reltol,atol,maxStepSize,
			maxNumSteps,maxErrTestsFails);
	sim_data=sim_function(cnolist,model,ode_parameters$parValues)
	sim_data=lapply(sim_data, function(x) as.matrix(x))
	return(sim_data);
}
