#include "include/nvector_serial.h"
#include "include/sundials_dense.h"
#include "include/sundials_band.h"
#include <stdio.h>

/*	Macro for accessing vector values. These macros work on a system of 
	indices starting with 1. To access the first entry of the vector x, 
	use the form Ith(x, 1)
*/
#define Ith(v,i) ( NV_DATA_S(v)[i - 1] )

/* IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   DENSE_ELEM macro in dense.h. DENSE_ELEM numbers rows and columns of a
   dense matrix starting from 0. 
*/
#define DenseIJth(A,i,j) DENSE_ELEM(A,i-1,j-1) 

// CVODES Functions

// Sample right hand side function
int rhs(realtype t, N_Vector y, N_Vector ydot, void *data)
{
  realtype y1, y2, y3;
  
  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  
  double L = 49.3;
  double a21 = 0.011; double a12 = 0.012;
  double a31 = 0.0039; double a13 = 0.000035;
  double a01 = 0.021; double a02 = 0.016;
  
  Ith(ydot, 1) = -(a01 + a21 + a31)*y1 + a12*y2 + a13*y3 + L;
  Ith(ydot, 2) = a21*y1 - (a02 + a12)*y2;
  Ith(ydot, 3) = a31*y1 - a13*y3;

  return(0);
}

// Sample dense Jacobian function
int densejac(long int N, DenseMat J, realtype t,
               N_Vector y, N_Vector fy, void *data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y1, y2, y3;
  
  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  
  DenseIJth(J,1,1) = RCONST(-0.04);
  DenseIJth(J,1,2) = RCONST(1.0e4)*y3;
  DenseIJth(J,1,3) = RCONST(1.0e4)*y2;
  DenseIJth(J,2,1) = RCONST(0.04); 
  DenseIJth(J,2,2) = RCONST(-1.0e4)*y3-RCONST(6.0e7)*y2;
  DenseIJth(J,2,3) = RCONST(-1.0e4)*y2;
  DenseIJth(J,3,2) = RCONST(6.0e7)*y2;

  return(0);
}


/*
 * g routine. Compute functions g_i(t,y) for i = 0,1. 
 */

int g(realtype t, N_Vector y, realtype *gout, void *data)
{
  realtype y1, y3;

  y1 = Ith(y,1); y3 = Ith(y,3);
  gout[0] = y1 - RCONST(0.0001);
  gout[1] = y3 - RCONST(0.01);

  return(0);
}
