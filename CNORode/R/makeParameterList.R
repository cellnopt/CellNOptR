makeParameterList <-
function(adjMat,namesSpecies)
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
				parValues[count]=1;
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

