#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv2 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-2.0.html
#
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/software.html
#
##############################################################################
# $Id$
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

