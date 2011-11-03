incidence2Adjacency <-
function(model)
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

