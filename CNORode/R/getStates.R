getStates <-
function(adjacency)
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

