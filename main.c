/* 
 * File:   main.c
 * Author: David
 *
 * Created on 18 de Junho de 2011, 17:41
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double normHill(double x,double n,double k);

typedef struct
{
    int n;
    int** adjMatrix;
    int* isState;
    int* boolTable;
    int* numInputs;
    int* numBits;
    int* inputValues;
    double* stateValues;
    int numStates;
    double* parameters;
} myODEdata;
/*
 * 
 */

int interpretMatrix(void* data)
{
    myODEdata* myData=(myODEdata*) data;
    //int** adjMatrix=data.adjMatrix;
    double res[(*myData).n];
    int j,i,x,y,number,countTableLine,tempCountTableLine,x1,y1,first,countInput;
    int countPar=0;
    double tempProd;
    //Loop through every column j in the matrix
    for (j = 0; j <(*myData).n; j++)
    {
        //if this column is from state and not an input
        res[j]=0;
        if((*myData).isState[j])
        {
           double hillFuncValues[(*myData).numInputs[j]];
           first=1;
           //start the truthTable line with 0;
           countTableLine=-1;
           //start the iterative sum with 0
           //Get the integer number representing a boolean function
           number=(*myData).boolTable[j];
           //Parse it to a truth Table
           x = y = 0;
           for(y =(*myData).numBits[j]-1; y >= 0; y--)
           {
               //Each bit here is a new truthTable Line
                countTableLine++;
                x = number / (1 << y);
                number = number - x * (1 << y);

                /*Example numBits[j]=4 boolTale[j]=3
                 * x1 x2 y
                 * 0  0  0 
                 * 0  1  1
                 * 1  0  1 
                 * 1  0  0
                 * 
                 * x=0 when i=0,x=1 when i=1; x=1 when i=2, x=0 i=3.
                 */
                tempCountTableLine=countTableLine;
                if(x)
                {
                    tempProd=1;
                    y1=(*myData).numInputs[j]-1;
                    //Initiate the product with the value 1
                    countInput=-1;
                    for (i = 0; i <(*myData).n; i++)
                    {
                        if((*myData).adjMatrix[i][j])
                        {
                            countInput++;
                            x1 = tempCountTableLine / (1 << y1);
                            tempCountTableLine = tempCountTableLine - x1 * (1 << y1);
                            y1--;

                            //if it is a state find the value from the state vector
                            if((*myData).isState[i])
                            {
                                if(first)
                                {
                                    hillFuncValues[countInput]=normHill(
                                            (*myData).stateValues[i],
                                            (*myData).parameters[countPar++],
                                            (*myData).parameters[countPar++]);
                                }
                                if(x1)
                                {
                                    tempProd*=hillFuncValues[countInput];
                                }
                                else
                                {
                                    tempProd*=(1-hillFuncValues[countInput]);
                                }
                            }
                                //otherwise get the value from the input input vector
                            else
                            {
                                if(first)
                                {
                                    hillFuncValues[countInput]=normHill(
                                            (*myData).inputValues[i],
                                            (*myData).parameters[countPar++],
                                            (*myData).parameters[countPar++]);
                                }
                                if(x1)
                                {
                                    tempProd*=hillFuncValues[countInput];
                                }
                                else
                                {
                                    tempProd*=(1-hillFuncValues[countInput]);
                                }
                            }
                          }
                        }
                        first=0;
                        res[j]+=tempProd;
                    }
                }
          
           res[j]*=(*myData).parameters[countPar++];
           printf("%f\n",res[j]);
           free(hillFuncValues);
           }
        }
     /*for(i=0;i<(*myData).n;i++){
    if((*myData).isState[i]){printf("%f\n",res[i]);}}
    return(0);*/
}

double normHill(double x,double n,double k)
{
    return(x);
}

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

int main(int argc, char** argv)
{
    int i,j,n;
    n=5;
    int **array=(int **)malloc(n*sizeof(int *));
    for(i = 0; i < n; i++){
        array[i] = malloc(n * sizeof(int));}

    int adjMatrix[5][5]={{0,0,1,1,0},{0,0,1,0,0},{0,0,0,0,1},{0,0,0,0,0},{0,0,0,0,0}};
    for (i = 0; i <n; i++){
        for (j = 0; j <n; j++){
            array[i][j]=adjMatrix[i][j];}}
    int *boolFunc=(int*)malloc(n*sizeof(int));
    boolFunc[0]=-1;boolFunc[1]=-1;boolFunc[2]=7;boolFunc[3]=1;boolFunc[4]=1;
    int *numInputs=(int*)malloc(n*sizeof(int));
    int *numBits=(int*)malloc(n*sizeof(int));
    numInputs[0]=-1;numInputs[1]=-1;numInputs[2]=2;numInputs[3]=1;numInputs[4]=1;
     int *isState=findStates(array,n);
    for (i = 0; i <n; i++){
        if(isState[i])
        {
            numBits[i]=pow(2,numInputs[i]);
        }
        else{numBits[i]=-1;}
    }
    
    int *inputValues=(int*)malloc(n*sizeof(int));
    inputValues[0]=1;inputValues[1]=1;inputValues[2]=-1;inputValues[3]=-1;inputValues[4]=-1;
    double *stateValues=(double*)malloc(n*sizeof(double));
    stateValues[0]=0;stateValues[1]=0;stateValues[2]=0.5;stateValues[3]=0.5;stateValues[4]=0.5;
    int numStates=0;
    for(i=0;i<n;i++){
    if(isState[i]){numStates++;}}
    double *parameters=(double*)malloc(11*sizeof(double));

    parameters[0]=3;parameters[1]=0.5;parameters[2]=3;parameters[3]=0.5;parameters[4]=1;
    parameters[5]=3;parameters[6]=0.5;parameters[7]=1;
    parameters[8]=3;parameters[9]=0.5;parameters[10]=1;
    
    myODEdata data;
    data.adjMatrix =array;
    data.isState=isState;
    data.n = n;
    data.boolTable=boolFunc;
    data.numInputs=numInputs;
    data.numBits=numBits;
    data.inputValues=inputValues;
    data.stateValues=stateValues;
    data.numStates=numStates;
    data.parameters=parameters;
    myODEdata *ptrData=malloc(sizeof(data));
    (*ptrData).n=n;
    (*ptrData).adjMatrix=array;
    (*ptrData).isState=isState;
    (*ptrData).boolTable=boolFunc;
    (*ptrData).numInputs=numInputs;
    (*ptrData).numBits=numBits;
    (*ptrData).inputValues=inputValues;
    (*ptrData).stateValues=stateValues;
    (*ptrData).numStates=numStates;
    (*ptrData).parameters=parameters;

    //void* ptr=&data;
    interpretMatrix(ptrData);
    //interpretMatrix(ptr);
    return (EXIT_SUCCESS);
}
