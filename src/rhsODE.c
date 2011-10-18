
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "CVODES/include/cvodes/cvodes.h"           /* prototypes for CVODES fcts. and consts. */
#include "CVODES/include/sundials/sundials_types.h" /* definition of type realtype */
#include "CNOStructure.h"
#include "CVODES/include/nvector/nvector_serial.h"  /* serial N_Vector types, fcts., and macros */

#define Ith(v,i) ( NV_DATA_S(v)[i] )

double normHill(double x,double n,double k);

int* decimal2binary(int decimal_value,int nBits);

int rhsODE(realtype t, N_Vector y, N_Vector ydot, void *data)
{
	CNOStructure* myData=(CNOStructure*) data;

	    int j,i,k;
	    int countPar=0;
	    double tempProd;
		double kHill,nHill;
	    double* hillFuncValues;
	    int countState=0;
	    int inputCount;
	    double test;
	    int* binary_value;

	    //Loop through every column j in the Graph adjacency matrix
	    for (j = 0; j <(*myData).nStates; j++)
	    {
	        if((*myData).isState[j])
	        {
	        	hillFuncValues= (double *)malloc((*myData).numInputs[j]*sizeof(double));
	        	Ith(ydot,countState)=0;

	           inputCount=0;
	           for (i = 0; i < (*myData).nRows; ++i)
	           {
	        	   if((*myData).adjacencyMatrix[i][j])
	        	   {
	        		  kHill=(*myData).odeParameters[countPar++];
	        		  nHill=(*myData).odeParameters[countPar++];

	        		   if((*myData).isState[i])
	        		   {
	        			   hillFuncValues[inputCount++]=
	        					   normHill(Ith(y,inputCount),kHill,nHill);
	        		   }
	        		   else
	        		   {
	        			   hillFuncValues[inputCount++]=
	        					   normHill((*myData).state_array[j],kHill,nHill);
	        		   }
	        	   }
	           }

	           for (i = 0; i < (*myData).numBits[j]; ++i)
	           {
	        	   if ((*myData).truthTables[j][i])
	        	   {
	        		   tempProd=1;
	        		   binary_value = decimal2binary(i-1,(*myData).numInputs[i]);
	        		   for (k = 0; k < (*myData).numInputs[j]; k++)
	        		   {
	        			   if(binary_value[i]==0)
	        			   {
	        				   tempProd*=(1-hillFuncValues[k]);
	        			   }
	        			   else tempProd*=hillFuncValues[k];
	        		   }
					free(binary_value);
					Ith(ydot,countState)+=tempProd;
	        	   }
	           }
	           free(hillFuncValues);
	           printf("%f\n",Ith(ydot,countState));
	           countState++;
	        }

	    }
	    return(0);
	}

