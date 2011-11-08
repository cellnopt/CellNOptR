# authors: TC
# This example provide an example and function to call a C function
# from CNORinterface shared library.


# You can test using :
# R --no-restore --no-save  < example.R
#
# or within a R shell by copy/paste the code that follows


# First read the MIDAS and SIF file to get some data

makeParameterList=function(adjMat,namesSpecies)
{
	isState=getStates(adjMat);
	numStates=as.integer(sqrt(length(adjMat)));
	parNames=c();
	parValues=c();
	count=0;
	for(j in 1:numStates)
	{ 
		if(isState[j])
		{
			jCol=adjMat[,j];
			inputs=which(as.logical(jCol));
			numInputs=length(inputs);
			for(i in 1:numInputs)
			{
				count=count+1;
				parNames[count]=paste(namesSpecies[inputs[i]],"_n_",namesSpecies[j],sep="");
				parValues[count]=3;
				count=count+1;
				parNames[count]=paste(namesSpecies[inputs[i]],"_k_",namesSpecies[j],sep="");
				parValues[count]=0.5;
			}
			count=count+1;
			parNames[count]=paste("tau_",namesSpecies[j],sep="");
			parValues[count]=1;
		}
	}
	parList=list(parNames=parNames,parValues=parValues);
	
	return(parList);
}

getStates=function(adjacency)
{
	nSpecies=dim(adjacency)[1];
	count=0;
	res=matrix(0,1,nSpecies);
	for(j in 1:nSpecies)
	{
		for(i in 1:nSpecies)
		{
			if(adjacency[i,j])
			{
				res[j]=1;
			}
		}
	}
	return(res);
}

incidence2Adjacency=function(model)
{ 
	incidence=model$interMat;
	nNodes=dim(incidence)[1];
	nEdges=dim(incidence)[2];
	adjacency=matrix(0,nNodes,nNodes);
	
	for(j in 1:nEdges)
	{
		for(i in 1:nNodes)
		{
			if(incidence[i,j]==1)
			{
				node1=i;
				for(k in 1:nNodes)
				{
					if(incidence[k,j]==-1)
					{
						adjacency[k,node1]=1
					}
				}
			}
		}
	}
	return(adjacency);  
}

interface <- function(cnolist, sif, indices, odeParameters, time=1,verbose=0)
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

library(CellNOptR)
##detach("package:testPackage");
#install.packages("CNORode_1.0.zip",repos=NULL);
library("CNORode");
setwd("C:/Users/davidh/workspace_java/CNOR_ode/tests/test3")
s = readSif('model.sif')
m = readMIDAS('initialData.csv')
cnolist = makeCNOlist(m, subfield=FALSE)
cnolist$timeSignals=seq(0,10)

indices <- indexFinder(cnolist, s, verbose = TRUE)
modelNCNOindices <- findNONC(s, indices, verbose = TRUE)
s <- cutNONC(s, modelNCNOindices);

adjMat=incidence2Adjacency(s)
odeParameters=makeParameterList(adjMat,s$namesSpecies)

# Finally, call the function
res = interface(cnolist, s, indices,odeParameters$parValue,verbose=TRUE);
out<-lapply(res,function(x) x[,indices$signals])
plotCNOlist(cnolist)
windows()
cnolist$valueSignals=out;
plotCNOlist(cnolist);

