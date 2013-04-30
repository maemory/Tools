/*=================================================================
 *
 * GEN_GLYPH_magnitude.cpp	Sample .MEX file corresponding to gen_glyph_magnitude.M
 *	        Solves for Reynolds stress tensor glyph given stress tensor and points at which to evaluate magnitude
 *
 * The calling syntax is:
 *
 *		[X,Y,Z,C] = gen_glyph(Rij,points,dim)
 *
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
#define	Rij_IN	    prhs[0]
#define	Xp_IN       prhs[1]
#define	Yp_IN       prhs[2]
#define	Zp_IN       prhs[3]

/* Output Arguments */
#define	X_OUT	plhs[0]
#define	Y_OUT	plhs[1]
#define	Z_OUT	plhs[2]
#define	C_OUT	plhs[3]

#define XO(i,j) xo[i+dimi*j]
#define YO(i,j) yo[i+dimi*j]
#define ZO(i,j) zo[i+dimi*j]
#define CO(i,j) co[i+dimi*j]
#define rij(i,j) rij[i+3*j]


#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

static void make_glyph(const double *rij,const int dimi,const int dimj,double *xo,double *yo,double *zo,double *co,double *X,double *Y,double *Z) ;

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    
    /* Check for proper number of arguments */
    if (nrhs != 4) { 
	mexErrMsgTxt("Four input arguments required (Rij, X, Y, Z)."); 
    } else if (nlhs != 4) {
	mexErrMsgTxt("Four output arguments required."); 
    } 
    
    /* Check the dimensions of Rij, must be 3x3 */
    mwSize m,n,mx,nx,my,ny,mz,nz; 

    m = mxGetM(Rij_IN); 
    n = mxGetN(Rij_IN);
    if (!mxIsDouble(Rij_IN) || mxIsComplex(Rij_IN) || m != 3 || n != 3) { 
        mexErrMsgTxt("Rij must be a 3x3 matrix."); 
    }

    /* Check the properties of points */
    mx = mxGetM(Xp_IN); 
    nx = mxGetN(Xp_IN);
    if (!mxIsDouble(Xp_IN) || mxIsComplex(Xp_IN) || mx < 1 || nx < 1) { 
        mexErrMsgTxt("X-points must be real numeric values."); 
    }
    
    my = mxGetM(Yp_IN); 
    ny = mxGetN(Yp_IN);
    if (!mxIsDouble(Yp_IN) || mxIsComplex(Yp_IN) || my < 1 || ny < 1) { 
        mexErrMsgTxt("Y-points must be real numeric values."); 
    }

    mz = mxGetM(Zp_IN); 
    nz = mxGetN(Zp_IN);
    if (!mxIsDouble(Zp_IN) || mxIsComplex(Zp_IN) || mz < 1 || nz < 1) { 
        mexErrMsgTxt("Z-points must be real numeric values."); 
    }    
    
    /* Check the dimensions*/ 
    if (mx != my || my != mz || nx != ny || ny != nz) { 
        mexErrMsgTxt("dimensions must be consistent across X,Y,Z."); 
    } 
    int dimi = (int)mx;
    int dimj = (int)nx;
    
    /* Create matrix for the return arguments */ 
    X_OUT = mxCreateDoubleMatrix(dimi, dimj, mxREAL); 
    Y_OUT = mxCreateDoubleMatrix(dimi, dimj, mxREAL); 
    Z_OUT = mxCreateDoubleMatrix(dimi, dimj, mxREAL); 
    C_OUT = mxCreateDoubleMatrix(dimi, dimj, mxREAL); 

    /* Assign pointers to the various parameters */
    double *xo, *yo, *zo, *co, *rij; 

    xo = mxGetPr(X_OUT);
    yo = mxGetPr(Y_OUT);
    zo = mxGetPr(Z_OUT);
    co = mxGetPr(C_OUT);
    rij = mxGetPr(Rij_IN);

    double *X,*Y,*Z;
    X = mxGetPr(Xp_IN);
    Y = mxGetPr(Yp_IN);
    Z = mxGetPr(Zp_IN);

    /* Do the actual computations in a subroutine */
    make_glyph(rij,dimi,dimj,xo,yo,zo,co,X,Y,Z);
    return;
    
}

static void make_glyph(const double *rij, const int dimi,const int dimj, double *xo,double *yo,double *zo,double *co,double *X,double *Y,double *Z) {

  double x,y,z,mag;

  for (int i=0; i<dimi; i++) {
    for (int j=0; j<dimj; j++) {
      x = X[i + dimi*j];
      y = Y[i + dimi*j];
      z = Z[i + dimi*j];

      mag = x*x*rij(0,0) + y*y*rij(1,1) + z*z*rij(2,2) + 2*x*y*rij(0,1) + 2*x*z*rij(0,2) + 2*y*z*rij(1,2);
 //     mag = 1.0;
      XO(i,j) = mag * x;
      YO(i,j) = mag * y;
      ZO(i,j) = mag * z;
      CO(i,j) = mag;
    }
  }
}
