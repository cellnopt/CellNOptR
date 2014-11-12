#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>

int** get_input_index(int** adjMat,int n,int* numInputs)
{
    int **input_index=(int**)malloc(n*sizeof(int*));
    int i,j,count;
    for (j = 0; j < n; ++j)
    {
    	input_index[j]=(int*)malloc(numInputs[j]*sizeof(int));
    	count=0;
    	for (i = 0; i < n; ++i)	if(adjMat[i][j])input_index[j][count++]=i;
	}
    return(input_index);
}
