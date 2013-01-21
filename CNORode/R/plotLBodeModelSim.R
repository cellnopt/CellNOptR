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
plotLBodeModelSim <-function
(
		cnolist,				    model,					    ode_parameters=NULL,
		indices=NULL,			    adjMatrix=NULL,			  	timeSignals=NULL,
		time=1,						verbose=0, 				    transfer_function=3,		
		reltol=1e-4,				atol=1e-3,				    maxStepSize=Inf,
		maxNumSteps=100000,			maxErrTestsFails=50,  		large=FALSE,          		
		nsplit=4
)
{

    if (class(cnolist)=="CNOlist"){cnolist = compatCNOlist(cnolist)}
	if(is.null(indices))indices=indexFinder(cnolist,model);
	if(is.null(adjMatrix))adjMatrix=incidence2Adjacency(model);
	if(is.null(ode_parameters))ode_parameters=createLBodeContPars(model);
	if(!is.null(timeSignals))cnolist$timeSignals=timeSignals;

  	states_index=which(as.logical(getStates(adjMatrix)));

	sim_data=getLBodeModelSim(cnolist,model,
			ode_parameters,indices,timeSignals,time,verbose,transfer_function,
			reltol,atol,maxStepSize,maxNumSteps,maxErrTestsFails);

	sim_data=lapply(sim_data,function(x) x[,states_index]);
	times=cnolist$timeSignals;
	cnolist$valueSignals=sim_data;
    cnolist$valueVariances = sim_data # variance should be set to NA but does not matter here 
	cnolist$namesSignals=model$namesSpecies[states_index];
	cnolist$namesCues=c(cnolist$namesStimuli,cnolist$namesInhibitors);

	cnolist$valueCues=cbind(cnolist$valueStimuli,cnolist$valueInhibitors);
	cnolist$valueCues=as.matrix(cnolist$valueCues);
	cnolist$valueCues[which(cnolist$valueCues>0)]=1;
	names(cnolist$valueCues)=cnolist$namesCues;


	if(large)
	{
		plotCNOlistLarge(CNOlist(cnolist),nsplit);
	}
	else
	{
		plotCNOlist(CNOlist(cnolist));
	}


  return(sim_data);
}
