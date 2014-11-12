/*
 * printInterMat.c
 *
 *  Created on: 7 Oct 2011
 *      Author: davidh
 */

#include <stdio.h>
#include <stdlib.h>
#include <R.h>

void printInterMat(int** interMat,int nRows,int nCols)
{
	int i,j;
	for (i = 0; i < nRows; ++i)
	{
		for (j = 0; j < nCols; ++j)
		{
			printf("%d\t",interMat[i][j]);
		}
		printf("\n");
	}

}
