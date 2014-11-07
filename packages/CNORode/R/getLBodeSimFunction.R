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
# $Id: getLBodeSimFunction.R 4446 2014-03-12 15:35:32Z cokelaer $
getLBodeSimFunction <-function
(
		cnolist1,				model1,					adjMatrix1,
		indices1, 				odeParameters1,			time1=1,
		verbose1=0,				transfer_function1=3,	reltol1=1e-4,	
		atol1=1e-3,				maxStepSize1=Inf,		maxNumSteps1=100000,
		maxErrTestsFails1=50, initial_state1=0.1
)
{
				
	simulate_logic_based_ode_model <- function
	(
			cnolist,	 						sif, 									odeParameters,
			indices=indices1,					adjMatrix=adjMatrix1, 					time=time1,
			verbose=verbose1, 					transfer_function=transfer_function1,	reltol=reltol1,
			atol=atol1,							maxStepSize=maxStepSize1,				maxNumSteps=maxNumSteps1,
			maxErrTestsFails=maxErrTestsFails1, initial_state=initial_state1
	)
	{	 
		interMat <- as.integer(as.vector(t(sif$interMat)))
		notMat <- as.integer(as.vector(t(sif$notMat)))
		adjMatrix <- as.integer(as.vector(t(adjMatrix)))
		nRows <- as.integer(dim(sif$interMat)[1])
		
		nCols <- as.integer(dim(sif$interMat)[2])
		verbose=as.integer(verbose);
	
		nPars <- as.integer(length(odeParameters))
	
		# cnolist
		timeSignals <- as.double(cnolist$timeSignals)
		valueInhibitors <- as.double(t(cnolist$valueInhibitors))
		valueSignals <- as.double(t(cnolist$valueSignals[[time]]))
		valueStimuli <- as.double(t(cnolist$valueStimuli))
		nTimes = as.integer(length(cnolist$timeSignals))
		
		# [[1]] allows to access to the first object in the list and retrieve its dimensions
		nExperiments = as.integer(dim(cnolist$valueSignals[[time]]))
		
		#indices
		nSignals <- as.integer(length(indices$signals))
		
		indexSignals <- as.integer(as.vector(indices$signals))
		nStimuli <- as.integer(length(indices$stimulated))
		indexStimuli <- as.integer(as.vector(indices$stimulated))
		nInhibitors <- as.integer(length(indices$inhibited))
		indexInhibitors <- as.integer(as.vector(indices$inhibited))
		transfer_function=as.integer(transfer_function);
		reltol=as.double(reltol);
		atol=as.double(atol);
		maxStepSize=as.double(maxStepSize);
		maxNumSteps=as.integer(maxNumSteps);
		maxErrTestsFails=as.integer(maxErrTestsFails);
		break_at_1st_fail=as.integer(0);
	    initial_state = initial_state

		res = .Call("sim_logic_ode",interMat,notMat,adjMatrix,nRows,nCols,nPars,timeSignals,
				valueInhibitors,valueSignals,valueStimuli,nTimes,nExperiments,nSignals,
				indexSignals,nStimuli,indexStimuli,nInhibitors,indexInhibitors,odeParameters,
				verbose,transfer_function,reltol,atol,maxStepSize,maxNumSteps,maxErrTestsFails,
				break_at_1st_fail, initial_state);
		return(res);
	}
	return(simulate_logic_based_ode_model);
}

