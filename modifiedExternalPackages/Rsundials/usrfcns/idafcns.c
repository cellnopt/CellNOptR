#include "include/nvector_serial.h"
#include "include/sundials_dense.h"

/* IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   DENSE_ELEM macro in dense.h. DENSE_ELEM numbers rows and columns of a
   dense matrix starting from 0. 
*/
#define DenseIJth(A,i,j) DENSE_ELEM(A,i-1,j-1) 



// Sample residual function
int resrob(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *rdata)
{
  realtype *yval, *ypval, *rval;

  yval = NV_DATA_S(yy); 
  ypval = NV_DATA_S(yp); 
  rval = NV_DATA_S(rr);
  
  rval[0] = -0.04*yval[0] + 1.0e4*yval[1]*yval[2] - ypval[0];
  rval[1] = 0.04*yval[0] - 1.0e4*yval[1]*yval[2] - 3.0e7*yval[1]*yval[1] - ypval[1];
  rval[2] = yval[0] + yval[1] + yval[2] - 1;

  return(0);
}

// Sample Jacobian function
int jacrob(long int Neq, realtype tt, N_Vector yy, N_Vector yp,
           N_Vector resvec, realtype cj, void *jdata, DenseMat JJ,
           N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
  realtype *yval;
  
  yval = NV_DATA_S(yy);

  DenseIJth(JJ,1,1) = RCONST(-0.04) - cj;
  DenseIJth(JJ,2,1) = RCONST(0.04);
  DenseIJth(JJ,3,1) = 1;
  DenseIJth(JJ,1,2) = RCONST(1.0e4)*yval[2];
  DenseIJth(JJ,2,2) = RCONST(-1.0e4)*yval[2] - RCONST(6.0e7)*yval[1] - cj;
  DenseIJth(JJ,3,2) = 1;
  DenseIJth(JJ,1,3) = RCONST(1.0e4)*yval[1];
  DenseIJth(JJ,2,3) = RCONST(-1.0e4)*yval[1];
  DenseIJth(JJ,3,3) = 1;

  return(0);
}

// Sample root finding function
int grob(realtype t, N_Vector yy, N_Vector yp, realtype *gout,
                void *g_data)
{
  realtype *yval, y1, y3;

  yval = NV_DATA_S(yy); 
  y1 = yval[0]; y3 = yval[2];
  gout[0] = y1 - RCONST(0.0001);
  gout[1] = y3 - RCONST(0.01);

  return(0);
}
