#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int** truth_table_index(int n,int** truth_tables, int* numBits)
{
    int **truth_table_index=(int**)malloc(n*sizeof(int*));
    int i,j,count_bits;
    for (i = 0; i < n; ++i)
    {
    	count_bits=0;
    	for (j = 0; j < numBits[i]; ++j)
    	{
    		if(truth_tables[i][j])count_bits++;
    	}
    	truth_table_index[i]=(int*)malloc(count_bits*sizeof(int));
    	count_bits=0;
    	for (j = 0; j < numBits[i]; ++j)
    	{
    		if(truth_tables_index[i][j])truth_table_index[i][count_bits++]=j;
    	}
	}
    return(input_index);
}
