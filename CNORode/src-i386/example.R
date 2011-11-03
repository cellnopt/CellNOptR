# authors: TC
# This example provide an example and function to call a C function
# from CNORinterface shared library.


# You can test using :
# R --no-restore --no-save  < example.R
#
# or within a R shell by copy/paste the code that follows


# First read the MIDAS and SIF file to get some data

library(CellNOptR)
s = readSif('../../CNOR/CellNOptR/inst/ToyModel/ToyPKNMMB.sif')
m = readMIDAS('../../CNOR/CellNOptR/inst/ToyModel/ToyDataMMB.csv')
cnolist = makeCNOlist(m, subfield=FALSE)
indices = indexFinder(cnolist, s)



# Second load the CNORinterface to C functions and create an interface to ease the wrapping
dyn.load("CNORinterface.so")

odeParameters <- c(2.203934e+000,	6.262312e-001,	1.359371e+000,	8.331686e-001,	8.413469e-001,	1.088145e+000,	5.515802e-001,	4.646723e-001,	1.968970e+000,	6.239473e-001,	8.709632e-001,	4.088349e+000,	4.171849e-001,	4.874472e+000,	8.243882e-001,	1.261108e-001,	2.448862e+000,	6.793143e-001,	2.056699e+000,	4.625627e-001,	5.942217e-001,	4.564655e+000,	4.420256e-001,	3.202982e+000,	1.512126e-001,	1.608986e-001,	2.505269e+000,	3.879579e-001,	3.688188e-001,	1.786170e+000,	7.341170e-001,	2.064854e-001,	4.098042e+000,	2.482593e-001,	9.666349e-001,	4.161417e+000,	4.353775e-001,	2.235785e-001)



interface <- function(cnolist, sif, indices, odeParameters, time=1){

  # 
  interMat <- as.integer(as.vector(sif$interMat))
  notMat <- as.integer(as.vector(sif$notMat) )
  nRows <- as.integer(dim(sif$interMat)[1])
  nCols <- as.integer(dim(sif$interMat)[2])
  


  # ode 
  nPars <- as.integer(length(odeParameters))

  # cnolist
  timeSignals <- as.double(cnolist$timeSignals)
  valueInhibitors <- as.double(cnolist$valueInhibitors )
  valueSignals <- as.double(as.vector(cnolist$valueSignals[[time]]))
  valueStimuli <- as.double(cnolist$valueStimuli)
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
  
  res = .C("CNORinterface", 
    interMat=interMat,
    notMat=notMat,   nRows=nRows,  nCols=nCols,
    valeuInhibitors=valueInhibitors, indexInhibitors=indexInhibitors, nInhibitors=nInhibitors,
    valueSignals=valueSignals, indexSignals=indexSignals, nSignals=nSignals,
    valueStimuli=valueStimuli, indexstimuli=indexStimuli, nStimuli=nStimuli,
    timeSignals=timeSignals, nTimes=nTimes,
    nExperiments=nExperiments,
    odeParams=as.double(odeParameters), nPars=nPars, output=
    )
  return(res)
}


# Finally, call the function

res = interface(cnolist, s, indices,  odeParameters)


