/*
 * getAdjacencyMatrix.c
 *
 *  Created on: 7 Oct 2011
 *      Author: davidh
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int*** get_support_truth_tables(int n,int *nInputs)
{
	int i,j,k;
	int*** support_truth_tables=(int***)malloc(n*sizeof(int**));

	for (i = 0; i < n; ++i)
	{
		support_truth_tables[i]=(int**)malloc(pow(2,nInputs[i])*sizeof(int*));
		for (j = 0; j < pow(2,nInputs[i]); ++j)
		{
			support_truth_tables[i][j]= decimal2binary(j,nInputs[i]);
		}
	}
	return(support_truth_tables);
}

