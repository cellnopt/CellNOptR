/*
 * findStates.c
 *
 *  Created on: 7 Oct 2011
 *      Author: davidh
 */

#include <stdio.h>
#include <stdlib.h>
#include <R.h>

int *findStates(int **adjMatrix, int n)
{
    int *stateVec=malloc(n*sizeof(int));
    int i,j;
    for (j = 0; j <n; j++)
    {
        stateVec[j]=0;
        for(i=0;i<n;i++)
        {
            if(adjMatrix[i][j])
            {
                stateVec[j]=1;
            }
        }
     }

   return(stateVec);
}
