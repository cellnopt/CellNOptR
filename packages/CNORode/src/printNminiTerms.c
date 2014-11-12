/*
 * printAdjMat.c
 *
 *  Created on: 7 Oct 2011
 *      Author: davidh
 */

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
void printNminiterms(int*** miniTerms,int* nInputs, int* nMiniterms, int nRows)
{
	int i,j,k;
	printf("\n");
	for (i = 0; i < nRows; i++)
	{
		printf("Number of miniterms:%d\n",nMiniterms[i]);
		printf("Number of n inputs:%d\n",nInputs[i]);
		printf("Species %d\n",i);
		for (j = 0; j <nInputs[i]; j++)
		{
			for (k = 0; k < nMiniterms[i]; ++k)
			{
				printf("%d\t",miniTerms[i][j][k]);
			}
			printf("\n");
		}

	}
}
