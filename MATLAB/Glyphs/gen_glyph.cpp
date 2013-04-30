/*=================================================================
 *
 * GEN_GLYPH.cpp	Sample .MEX file corresponding to gen_glyph.M
 *	        Solves for Reynolds stress tensor glyph given stress tensor
 *
 * The calling syntax is:
 *
 *		[X,Y,Z,C] = gen_glyph(Rij,dim)
 *
 *  You may also want to look at the corresponding M-code, yprime.m.
 *  Rij is the stress tensor to pass in
 *  dim is the dimension of the sphere used to generate the glyph
 * This is a MEX-file for MATLAB.  
 * Copyright 1984-2006 The MathWorks, Inc.
 *
 *=================================================================*/
/* $Revision: 1.10.6.4 $ */
#include <math.h>
#include "matrix.h"
#include "mex.h"

/* Input Arguments */
#define	Rij_IN	prhs[0]
#define	dim_IN	prhs[1]

/* Output Arguments */
#define	X_OUT	plhs[0]
#define	Y_OUT	plhs[1]
#define	Z_OUT	plhs[2]
#define	C_OUT	plhs[3]

#define XO(i,j) xo[i+dim*j]
#define YO(i,j) yo[i+dim*j]
#define ZO(i,j) zo[i+dim*j]
#define CO(i,j) co[i+dim*j]
#define rij(i,j) rij[i+3*j]


#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

static void make_glyph(const double *rij, const int dim,double *xo,double *yo,double *zo,double *co,double *X,double *Y,double *Z) ;

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    
    /* Check for proper number of arguments */
    if (nrhs != 2) { 
	mexErrMsgTxt("Two input arguments required."); 
    } else if (nlhs != 4) {
	mexErrMsgTxt("Four output arguments required."); 
    } 
    
    /* Check the dimensions of Rij, must be 3x3 */
    mwSize m,n; 

    m = mxGetM(Rij_IN); 
    n = mxGetN(Rij_IN);
    if (!mxIsDouble(Rij_IN) || mxIsComplex(Rij_IN) || m != 3 || n != 3) { 
	mexErrMsgTxt("Rij must be a 3x3 matrix."); 
    }

    /* Check the properties of dim */ 
    if (!mxIsDouble(dim_IN) || mxIsComplex(dim_IN) || (int)mxGetScalar(dim_IN) < 1) { 
	mexErrMsgTxt("dim must be an integer greater than 0."); 
    } 
    int dim = (int)mxGetScalar(dim_IN);

    /* Create matrix for the return arguments */ 
    X_OUT = mxCreateDoubleMatrix(dim, dim, mxREAL); 
    Y_OUT = mxCreateDoubleMatrix(dim, dim, mxREAL); 
    Z_OUT = mxCreateDoubleMatrix(dim, dim, mxREAL); 
    C_OUT = mxCreateDoubleMatrix(dim, dim, mxREAL); 

    /* Assign pointers to the various parameters */
    double *xo, *yo, *zo, *co, *rij; 

    xo = mxGetPr(X_OUT);
    yo = mxGetPr(Y_OUT);
    zo = mxGetPr(Z_OUT);
    co = mxGetPr(C_OUT);
    rij = mxGetPr(Rij_IN);

    /* load in Rij matrix */
//     for (int i=0; i<m; i++) {
//       for (int j=0; j<n; j++) {
//         rij[m][n] = mxGetPr(Rij_IN)[i + j*m];
//       }
//     }

    /* generate sphere with dim x dim points */
    mxArray *lhs[3], *x;
    x = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(x) = (double)(dim-1);
    mexCallMATLAB(3, lhs, 1, &x, "sphere");

    double *X,*Y,*Z;
    X = mxGetPr(lhs[0]);
    Y = mxGetPr(lhs[1]);
    Z = mxGetPr(lhs[2]);

    /* Do the actual computations in a subroutine */
    make_glyph(rij,dim,xo,yo,zo,co,X,Y,Z);
    return;
    
}

static void make_glyph(const double *rij, const int dim,double *xo,double *yo,double *zo,double *co,double *X,double *Y,double *Z) {

  double x,y,z,mag;

  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      x = X[i + dim*j];
      y = Y[i + dim*j];
      z = Z[i + dim*j];

      mag = x*x*rij(0,0) + y*y*rij(1,1) + z*z*rij(2,2) + 2*x*y*rij(0,1) + 2*x*z*rij(0,2) + 2*y*z*rij(1,2);
 //     mag = 1.0;
      XO(i,j) = mag * x;
      YO(i,j) = mag * y;
      ZO(i,j) = mag * z;
      CO(i,j) = mag;
    }
  }
}
