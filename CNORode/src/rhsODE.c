
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <cvodes/cvodes.h>          /* prototypes for CVODES fcts. and consts. */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include "CNOStructure.h"
#include <nvector/nvector_serial.h>/* serial N_Vector types, fcts., and macros */

#define Ith(v,i) ( NV_DATA_S(v)[i] )

int* decimal2binary(int decimal_value,int nBits);

int rhsODE(realtype t, N_Vector y, N_Vector ydot, void *data)
{
		int i,j,k;
		CNOStructure* myData=(CNOStructure*) data;
	    int countPar=0;
	    double tempProd;
		double kHill,nHill;
	    double hillFuncValues[(*myData).maxNumInputs];
	    int countState=0;
	    int inputCount;

	    //Loop through every column j in the Graph adjacency matrix
	    for (j = 0; j <(*myData).nRows; j++)
	    {
	        if((*myData).isState[j])
	        {
	          // hillFuncValues= (double*)malloc((*myData).numInputs[j]*sizeof(double));
	           Ith(ydot,countState)=0;
	           inputCount=0;
	           for (i = 0; i < (*myData).nRows; ++i)
	           {
	        	   if((*myData).adjacencyMatrix[i][j])
	        	   {
	        		  nHill=(*myData).odeParameters[countPar++];
	        		  kHill=(*myData).odeParameters[countPar++];

	        		   if((*myData).isState[i])
	        		   {
	        			   hillFuncValues[inputCount++]=
	        					   (*myData).transfer_function(Ith(y,(*myData).state_index[i]),nHill,kHill);
	        		   }
	        		   else
	        		   {
	        			   hillFuncValues[inputCount++]=
	        					   (*myData).transfer_function((*myData).state_array[i],nHill,kHill);
	        		   }
	        	   }
	           }

	           //For every bit in the truth table
	           for (i = 0; i < (*myData).numBits[j]; ++i)
	           {
	        	   if ((*myData).truthTables[j][i])
	        	   {
	        		   tempProd=1;
	        		   for (k = 0; k < (*myData).numInputs[j]; k++)
	        		   {

	        			   if((*myData).support_truth_tables[j][i][k]==0)
	        			   {
	        				   tempProd*=(1-hillFuncValues[k]);
	        			   }
	        			   else tempProd*=hillFuncValues[k];
	        		   }

	        		   Ith(ydot,countState)+=tempProd;
	        	   }
	           }
	           Ith(ydot,countState)=
	        		   (Ith(ydot,countState)-Ith(y,countState))
	        		   	   *(*myData).odeParameters[countPar++]
	        		   	   	   *(1-(*myData).inhibitor_array[j]);

	           countState++;
	        }
	    }

	    return(0);
	}

