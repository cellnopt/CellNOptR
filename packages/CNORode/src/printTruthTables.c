/*
 * printAdjMat.c
 *
 *  Created on: 7 Oct 2011
 *      Author: davidh
 */

#include <stdio.h>
#include <stdlib.h>
#include <R.h>


void printTruthTables(int** truthTables,int* nBits, int nRows)
{
	int i,j;
	printf("-----------------------------\n");
		for (i = 0; i < nRows; i++)
		{
			for (j = 0; j < nBits[i]; ++j)
			{
				printf("%d \n",truthTables[i][j]);
			}
			printf("------------------------\n");
		}

}
