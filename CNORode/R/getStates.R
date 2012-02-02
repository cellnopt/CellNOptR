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

