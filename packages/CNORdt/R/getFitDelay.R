# This file is part of the CNO software
# 
# Copyright (c) 2011-2012 - EBI
# 
# File author(s): CNO developers (cno-dev@ebi.ac.uk)
# 
# Distributed under the GPLv2 License.  See accompanying file LICENSE.txt or copy at
# http://www.gnu.org/licenses/gpl-2.0.html
# 
# CNO website: http://www.ebi.ac.uk/saezrodriguez/software.html
# 
# $Id: getFitDelay.R 4535 2014-04-23 12:07:41Z aidanmac $

getFitDelay <- function(CNOlist, model, simList, indexList, sizeFac=1e-04, NAFac=1, nInTot, boolUpdates) {
  
  if ((class(CNOlist)=="CNOlist") == FALSE) {
    CNOlist = CellNOptR::CNOlist(CNOlist)
  }
  
  # dimensions, time points
  nTimes = length(CNOlist@timepoints)
  sigs = dim(CNOlist@signals[[1]])
  simT0 = simulatorT0(CNOlist, model, simList, indexList)
  
  # interpolate experimental data so it can be compared to boolean simulation  
  splineStore = list()
  splineAdd = 1
  
  for (nExper in 1:dim(CNOlist@signals[[1]])[1]) {
    for (nSig in 1:dim(CNOlist@signals[[1]])[2]) {
      yTest = c()
      for (a in 1:nTimes) {
        yTest = c(yTest, CNOlist@signals[[a]][nExper, nSig])
      }
      if (!is.na(yTest[1])) {
        cS = splinefun(CNOlist@timepoints, yTest)
        splineStore[splineAdd] = list(cS)
      } else {
        cS = splinefun(CNOlist@timepoints, rep(0, nTimes))
        splineStore[splineAdd] = list(cS)
      }
      splineAdd = splineAdd + 1
    }
  }
  
  ySilico = array(dim=c(dim(CNOlist@signals[[1]])[1], length(indexList$signals), boolUpdates))
  # produce x-coordinates by taking the range of experimental time and dividing by boolUpdates
  xCoords = seq(0, by=boolUpdates/CNOlist@timepoints[length(CNOlist@timepoints)], length.out=boolUpdates)
  
  # interpolate data according to boolUpdates
  count = 1
  for(nExper in 1:dim(ySilico)[1]) {
    for(nSig in 1:dim(ySilico)[2]) {			             
      yOut = splineStore[[count]](xCoords)
      ySilico[nExper,nSig,] = yOut;
      count = count+1
    }
  }
  
  # where is the feedback where strong interactions can be enforced?
  strongWeak = rep(0,length(model$reacID))
  negEdges = feedbackWrapper(model) # the index of edges that can be changed
  strongWeak[negEdges] = 1
  
  # optimization
  delayFinder <- function(optimStringF) {
    
    delayThresh = optimStringF[1:length(model$reacID)] # integer section of optimString
    strongWeak[negEdges] = optimStringF[(length(model$reacID)+1):length(optimStringF)] # binary section of optimString
    simResults = simulatorDelay(CNOlist, model, simList, indexList,
                                boolUpdates=boolUpdates, delayThresh=delayThresh,
                                strongWeak=strongWeak, prevSim=simT0)
    simResults = convert2array(simResults, dim(CNOlist@signals[[1]])[1], length(model$namesSpecies), boolUpdates)
    simResults = simResults[,indexList$signals,]
    
    errorVector = as.vector(ySilico) - as.vector(simResults)
    sse = sum(errorVector^2)
    return(sse)			
  }
  
  optimString = c(rep(1,length(model$reacID)), rep(0,length(negEdges))) # initialize parameter vector to optimize
  problem = list(
    f=delayFinder, # define function to minimize
    x_O = optimString, # define initial value of vector
    x_L = c(rep(0,length(model$reacID)), rep(0,length(negEdges))), # lower bounds on parameters
    x_U = c(rep(10,length(model$reacID)), rep(1,length(negEdges))), # upper bounds on parameters
    int_var = length(model$reacID) + length(negEdges)
  )
  
  opts = list(maxtime=15, local_solver=0, dim_refset=6, ndiverse=50) # essR options
  est = essR(problem, opts)	# optimize
  
  bestDelay = est$xbest[1:length(model$reacID)] # extract the delay vector
  bestSW = est$xbest[(length(model$reacID)+1):length(est$xbest)] # extract the strongWeak vector
  strongWeak[negEdges] = bestSW # apply strongWeak best to full vector
  
  simBest = simulatorDelay(CNOlist, model, simList, indexList,
                           boolUpdates=boolUpdates, delayThresh=bestDelay,
                           strongWeak=strongWeak, prevSim=simT0)
  simBest = convert2array(simBest, dim(CNOlist@signals[[1]])[1], length(model$namesSpecies), boolUpdates)
  simBest = simBest[,indexList$signals,]
  plotCNOlist(plotData(CNOlistPB, simBest))
  
  diff <- (simBest - ySilico)
  r <- diff^2
  deviationPen <- sum(r[!is.na(r)])/nTimes
  NAPen <- NAFac * length(which(is.na(simBest)))
  nDataPts <- dim(CNOlist@signals[[1]])[1] * dim(CNOlist@signals[[1]])[2] * nTimes
  nInputs <- length(which(model$interMat == -1))
  
  # nInTot: number of inputs of expanded model nInputs: number of inputs of cut model
  sizePen <- (nDataPts * sizeFac * nInputs)/nInTot
  
  score <- deviationPen + NAPen + sizePen
  return(list(score=score, estimate=est$xbest, xCoords=xCoords, yInter=ySilico, simResults=simBest))
}
