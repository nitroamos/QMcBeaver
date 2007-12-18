/**
   This has been modified to use Array2D
*/

/* 
 * numerics/linalg/svdcmp.cc
 * 
 * Singular value decomposition. 
 * 
 * Copyright (c) 2003--2004 by Wolfgang Wieser (wwieser@gmx.de)
 * 
 * The SVD algorithm is based on a routine by G.E. Forsythe, 
 * M.A. Malcolm, and C.B. Moler which in turn is based on 
 * the original algorithm of G.H. Golub and C. Reinsch. 
 * 
 * This file may be distributed and/or modified under the terms of the
 * GNU General Public License version 2 as published by the Free Software
 * Foundation.
 * 
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 * 
 */ 

#ifndef SVDCMP_H
#define SVDCMP_H

#include "Array2D.h"
#include "math.h"

template<typename T>inline T SignOfNeg(T a,T b)
	{  return(b<0.0 ? fabs(a) : -fabs(a));  }

template<typename T>
T HYPOT(T a, T b){
  double absa = fabs(a);
  double absb = fabs(b);
  if(absa > absb)
    return absa*sqrt(1.0 + absb*absb/absa/absa);
  else
    return (absb == 0.0 ? 0.0 : absb*sqrt(1.0 + absa*absa/absb/absb) );
}

// Singular value decomposition. 
// Decompose the m x n - matrix A into U, V, W so that 
// A = U * W * transpose(V). 
// Thereby, U is an column-orthonormal m x n - matrix and replaces 
//              the input matrix A once the function returns,
//          W is a n x n diagonal matrix stored in the array w and 
//          V is an orthogonal n x n - matrix. 
// The elements W[] are called the "singular values" of A. 
// maxiter: Limit number of iterations. 
// Return value: 
//   0 -> OK
//   1 -> bailed out (no convergence after maxiter-many iterations)
//  -2 -> illegally-sized matrices 
template<typename T>int _SVDecompose(Array2D<T> &A, Array1D<T> &W, Array2D<T> &V,
				     int maxiter);

// Perform SVD (forward/) back substitution. 
// The SCD solution for the (under/overdetermined) equation A * x = b 
// can be obtained by first SVDecomposing A into U,W,V and then 
// using the back substitution to calculate x from b. 
// Return value: 
//   0 -> OK
//  -2 -> illegally-sized matrices 
template<typename T>int _SVDFwBackSubst(const Array2D<T> &U, const Array1D<T> &W,
					const Array2D<T> &V, Array2D<T> &inv);

int SVDecompose(Array2D<double> &a,Array1D<double> &w,Array2D<double> &v,int maxiter);
int SVDecompose(Array2D<float> &a,Array1D<float> &w,Array2D<float> &v,int maxiter);

int SVDFwBackSubst(const Array2D<double> &u, const Array1D<double> &w,
		   const Array2D<double> &v, Array2D<double> &inv);
int SVDFwBackSubst(const Array2D<float> &u, const Array1D<float> &w,
		   const Array2D<float> &v, Array2D<float> &inv);


#endif


