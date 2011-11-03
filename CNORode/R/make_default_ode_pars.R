make_default_ode_pars <-
function(interMat)
{
	adjMat=getAdjMat(interMat);
	n_nodes=dim(adjMat)[1];
	ode_parameters=c();
	
	for (j in 1:n_nodes) 
	{
		isState=FALSE;
		for (i in 1:n_nodes) 
		{
			if(adjMat[i,j])
			{
				ode_parameters=c(ode_parameters,3);
				ode_parameters=c(ode_parameters,0.5);
				isState=TRUE;
			}	
		}
		if(isState)ode_parameters=c(ode_parameters,1);
	}
	return(ode_parameters)
}

