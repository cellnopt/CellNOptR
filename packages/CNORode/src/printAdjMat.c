/*
 * printAdjMat.c
 *
 *  Created on: 7 Oct 2011
 *      Author: davidh
 */
#include <R.h>
#include <stdio.h>
#include <stdlib.h>

void printAdjMat(int** adjMat,int nNodes)
{
	int i,j;

	for (i = 0; i < nNodes; i++)
	{
		for (j = 0; j < nNodes; j++)
		{
			printf("%d\t",adjMat[i][j]);
		}
		printf(";\n");
	}
}
