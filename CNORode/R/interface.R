interface <-
function(cnolist, sif, indices, odeParameters, time=1,verbose=0)
{	 
	interMat <- as.integer(as.vector(t(sif$interMat)))
	notMat <- as.integer(as.vector(t(sif$notMat)))
	nRows <- as.integer(dim(sif$interMat)[1])
	nCols <- as.integer(dim(sif$interMat)[2])
	verbose=as.integer(verbose);
	# ode 
	nPars <- as.integer(length(odeParameters))
	
	# cnolist
	timeSignals <- as.double(cnolist$timeSignals)
	valueInhibitors <- as.double(t(cnolist$valueInhibitors))
	valueSignals <- as.double(as.vector(cnolist$valueSignals[[time]]))
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
	
	res = .Call("sim_logic_ode",interMat,notMat,nRows,nCols,nPars,timeSignals,
			valueInhibitors,valueSignals,valueStimuli,nTimes,nExperiments,nSignals,
			indexSignals,nStimuli,indexStimuli,nInhibitors,indexInhibitors,
			odeParameters,verbose);
	
	return(res);
	
}

