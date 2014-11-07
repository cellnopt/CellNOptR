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
# $Id: getStates.R 3157 2013-01-09 16:04:09Z cokelaer $
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

