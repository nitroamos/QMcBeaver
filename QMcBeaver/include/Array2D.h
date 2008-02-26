//            QMcBeaver
//
//         Constructed by
//
//     Michael Todd Feldmann
//              and
//   David Randall "Chip" Kent IV
//
// Copyright 2000.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

/**************************************************************************
This SOFTWARE has been authored or contributed to by an employee or 
employees of the University of California, operator of the Los Alamos 
National Laboratory under Contract No. W-7405-ENG-36 with the U.S. 
Department of Energy.  The U.S. Government has rights to use, reproduce, 
and distribute this SOFTWARE.  Neither the Government nor the University 
makes any warranty, express or implied, or assumes any liability or 
responsibility for the use of this SOFTWARE.  If SOFTWARE is modified 
to produce derivative works, such modified SOFTWARE should be clearly 
marked, so as not to confuse it with the version available from LANL.   
 
Additionally, this program is free software; you can distribute it and/or 
modify it under the terms of the GNU General Public License. Accordingly, 
this program is  distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY;  without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A  PARTICULAR PURPOSE.  See the GNU General Public License 
for more details. 
**************************************************************************/

#ifndef Array2D_H
#define Array2D_H

#include <iostream>
#include <iomanip>
#include <typeinfo>
#include "Complex.h"
#include <assert.h>
#include "Array1D.h"
#include "fastfunctions.h"
#include "cppblas.h"

/*
  Kahan's method is very slightly more accurate with
  single precision. It is implemented for curiousity's sake,
  since it's only available when gemm libraries are not linked.
  http://en.wikipedia.org/wiki/Kahan_summation_algorithm
*/
static const bool USE_KAHAN = false;

/*
  These are function prototypes for the LAPACK
  libraries that might be linked in.
  
  You may or may not need to modify the function names
  depending on the local fortran compiler's calling convention.
*/
extern "C"
{
  void FORTRAN_FUNC(dgesv)(int* N, int* NRHS,
	      double *A, int* lda, int* ipiv,
	      double *B, int* ldb,
	      int* info);
  void FORTRAN_FUNC(sgesv)(int* N, int* NRHS,
	      float *A, int* lda, int* ipiv,
	      float *B, int* ldb,
	      int* info);
  void FORTRAN_FUNC(dggev)(const char* jobvl, const char* jobvr, int* n,
	      double *a, int* lda, double *b, int* ldb,
	      double *alphar, double *alphai, double *beta,
	      double *vl, int* ldvl,
	      double *vr, int* ldvr,
	      double *work, int* lwork,
	      int *info);
}

//Array2D is included in all the classes where this change is relevant
#if defined SINGLEPRECISION || defined QMC_GPU
typedef float  qmcfloat;
const static  float REALLYTINY = 1e-35f;
#else
typedef double qmcfloat;
const static double REALLYTINY = 1e-300;
#endif

using namespace std;

/**
   A 2-dimensional template for making arrays.  All of the memory allocation
   and deallocation details are dealt with by the class. Some of the operators
   require the Array2D to contain doubles or floats.
*/

template <class T> class Array2D
{
 private:
  /**
     Number of elements in the array's first dimension.
  */
  
  int n_1;
  
  
  /**
     Number of elements in the array's second dimension.
  */
  
  int n_2;
  
  /**
     Array containing the data.
  */
  T * pArray;
  
  /**
     These three arrays are used in the determinant/inverse
     calculations.
  */
  Array1D<T> diCol;
  Array1D<int> diINDX;
  Array1D<T> diVV;
  
 public:
  /**
     Gets a pointer to an array containing the array elements.  
     The ordering of this array is NOT specified.  
  */
  T* array()
    {
      return pArray;
    }
  
  /**
     Gets the number of elements in the array's first dimension.
     
     @return number of elements in the array's first dimension.
  */
  int dim1() const
    {
      return n_1;
    }
  
  /**
     Gets the number of elements in the array's second dimension.
     
     @return number of elements in the array's second dimension.
  */
  int dim2() const
    {
      return n_2;
    }
  
  /**
     Gets the total number of elements in the array.
     
     @return total number of elements in the array.
  */
  int size() const
    {
      return n_1*n_2;
    }
  
  /**
     Allocates memory for the array.
     
     @param i size of the array's first dimension.
     @param j size of the array's second dimension.
  */
  
  void allocate(int i, int j)
    {
#ifdef QMC_DEBUG
      if( i < 0 || j < 0)
	{
	  cerr << "Error: invalid dimensions in Array2D::allocate\n"
	       << " i = " << i << endl
	       << " j = " << j << endl;
	  assert(0);
	}
#endif
      if( n_1 != i || n_2 != j )
	{
	  deallocate();
	  
	  n_1 = i;
	  n_2 = j;
	  
	  if(n_1 >= 1 && n_2 >= 1)
	    {
	      pArray = new T[ n_1*n_2 ];
	    }
	  else
	    {
	      pArray = 0;
	    }
	}
    }
  
  /**
     Deallocates memory for the array.
  */
  
  void deallocate()
    {
      if(pArray != 0)
	delete [] pArray;
      pArray = 0;
      
      n_1 = 0;
      n_2 = 0;
      
      diCol.deallocate();
      diINDX.deallocate();
      diVV.deallocate();
    }
  
  void setupMatrixMultiply(const Array2D<T> & rhs, Array2D<T> & result,
			   const bool rhsIsTransposed) const
    {
      if(rhsIsTransposed)
	{
#ifdef QMC_DEBUG
	  if(n_2 != rhs.n_2)
	    {
	      cerr << "ERROR: Transposed Matrix multiplication: " << n_1 << "x"
		   << n_2 << " * " << rhs.n_2 << "x" << rhs.n_1 << endl;
	      exit(1);
	    }
	  if(result.dim1() != n_1 || result.dim2() != rhs.n_1)
	    {
	      cerr << "Warning: we had to reallocate result: "
		   << result.dim1() << "x" << result.dim2() << " to "
		   << n_1 << "x" << rhs.n_1 << endl;
	    }
#endif
	  result.allocate(n_1,rhs.n_1);
	}
      else
	{
#ifdef QMC_DEBUG
	  if(n_2 != rhs.n_1)
	    {
	      cerr << "ERROR: Matrix multiplication: " << n_1 << "x"
		   << n_2 << " * " << rhs.n_1 << "x" << rhs.n_2 << endl;
	      exit(1);
	    }
	  if(result.dim1() != n_1 || result.dim2() != rhs.n_2)
	    {
	      cerr << "Warning: we had to reallocate result: "
		   << result.dim1() << "x" << result.dim2() << " to "
		   << n_1 << "x" << rhs.n_2 << endl;
	    }
#endif
	  result.allocate(n_1,rhs.n_2);
	}
    }
  
#ifdef USEBLAS
  
  void gemm(const Array2D<double> & rhs, Array2D<double> & result,
	    const bool rhsIsTransposed) const
    {
      setupMatrixMultiply(rhs,result,rhsIsTransposed);
      
      //MxN = alpha * MxK * KxN + beta * MxN
      double alpha = 1.0;
      double beta  = 0.0;
      
#ifdef USE_CBLAS
      int M        = result.n_1;
      int N        = result.n_2;
      int K        = n_2;
      int lda      = n_2;
      int ldb      = rhs.n_2;
      int ldc      = result.n_2;
      
      CBLAS_TRANSPOSE myTrans = rhsIsTransposed ? CblasTrans : CblasNoTrans;
      cblas_dgemm(CBLAS_ORDER(CblasRowMajor),
		  CBLAS_TRANSPOSE(CblasNoTrans), myTrans,
		  M, N, K,
		  alpha, pArray, lda,
		  rhs.pArray, ldb,
		  beta, result.pArray, ldc);
#else
      /*
	since fortran is column major, we need to mess around a bit
	to trick it into handling our row major data correctly.
	a row major matrix as a transposed column major matrix.
      */
      int M        = result.n_2;
      int N        = result.n_1;
      int K        = n_2;
      int lda      = rhs.n_2;
      int ldb      = n_2;
      int ldc      = result.n_2;
      
      char transa = rhsIsTransposed ? 'T' : 'N';
      char transb = 'N';
      
      FORTRAN_FUNC(dgemm)(&transa, &transb,		    
	     &M, &N, &K,
	     &alpha,
	     rhs.pArray, &lda,
	     pArray, &ldb,
	     &beta,
	     result.pArray, &ldc);
#endif
    }
  
  void gemm(const Array2D<float> & rhs, Array2D<float> & result,
	    const bool rhsIsTransposed) const
    {
      setupMatrixMultiply(rhs,result,rhsIsTransposed);
      
      //MxN = alpha * MxK * KxN + beta * MxN
      float alpha = 1.0;
      float beta  = 0.0;
      
#ifdef USE_CBLAS
      int M        = result.n_1;
      int N        = result.n_2;
      int K        = n_2;
      int lda      = n_2;
      int ldb      = rhs.n_2;
      int ldc      = result.n_2;
      
      CBLAS_TRANSPOSE myTrans = rhsIsTransposed ? CblasTrans : CblasNoTrans;
      cblas_sgemm(CBLAS_ORDER(CblasRowMajor),
		  CBLAS_TRANSPOSE(CblasNoTrans), myTrans,
		  M, N, K,
		  alpha, pArray, lda,
		  rhs.pArray, ldb,
		  beta, result.pArray, ldc);
#else
      /*
	since fortran is column major, we need to mess around a bit
	to trick it into handling our row major data correctly.
	a row major matrix as a transposed column major matrix.
      */
      int M        = result.n_2;
      int N        = result.n_1;
      int K        = n_2;
      int lda      = rhs.n_2;
      int ldb      = n_2;
      int ldc      = result.n_2;
      
      char transa = rhsIsTransposed ? 'T' : 'N';
      char transb = 'N';
      
      FORTRAN_FUNC(sgemm)(&transa, &transb,		    
	     &M, &N, &K,
	     &alpha,
	     rhs.pArray, &lda,
	     pArray, &ldb,
	     &beta,
	     result.pArray, &ldc);
#endif
    }
  
  /**
     This uses ATLAS to calculate a dot product between all the elements in one
     Array2D and all the elements in another Array2D
  */
  double dotAllElectrons(const Array2D<double> & rhs)
    {
#ifdef USE_CBLAS
      return cblas_ddot(n_1*n_2, pArray, 1, rhs.pArray, 1);
#else
      int inc = 1;
      int num = n_1 * n_2;
      return FORTRAN_FUNC(ddot)(&num, pArray, &inc, rhs.pArray, &inc);
#endif
    }
  
  float dotAllElectrons(const Array2D<float> & rhs)
    {
#ifdef USE_CBLAS
      return cblas_sdot(n_1*n_2, pArray, 1, rhs.pArray, 1);
#else
      int inc = 1;
      int num = n_1 * n_2;
      return FORTRAN_FUNC(sdot)(&num, pArray, &inc, rhs.pArray, &inc);
#endif
    }
  
  /**
     This uses ATLAS to calculate a dot product between all the elements in one
     row from one Array2D and one row in another Array2D
  */
  double dotOneElectron(const Array2D<double> & rhs, int whichElectron)
    {
#ifdef USE_CBLAS
      return cblas_ddot(n_1, pArray + whichElectron*n_2, 1, rhs.pArray + whichElectron*n_2, 1);
#else
      int inc = 1;
      return FORTRAN_FUNC(ddot)(&n_1, pArray + whichElectron*n_2, &inc, rhs.pArray + whichElectron*n_2, &inc);
#endif
    }
  
  float dotOneElectron(const Array2D<float> & rhs, int whichElectron)
    {
#ifdef USE_CBLAS
      return cblas_sdot(n_1, pArray + whichElectron*n_2, 1, rhs.pArray + whichElectron*n_2, 1);
#else
      int inc = 1;
      return FORTRAN_FUNC(sdot)(&n_1, pArray + whichElectron*n_2, &inc, rhs.pArray + whichElectron*n_2, &inc);
#endif	
    }

  void operator+=( const Array2D<double> & rhs)
    {
#ifdef USE_CBLAS
      cblas_daxpy(n_1*n_2, 1.0, rhs.pArray,1, pArray, 1);
#else
      int num = n_1*n_2;
      int inc = 1;
      double a = 1.0;
      FORTRAN_FUNC(daxpy)(&num, &a, rhs.pArray, &inc, pArray, &inc);
#endif
    }

#else
  void gemm(const Array2D<T> & rhs, Array2D<T> & result,
	    const bool rhsIsTransposed) const
    {
      setupMatrixMultiply(rhs,result,rhsIsTransposed);
      
      //MxN = MxK * KxN
      int M = n_1;
      int N = rhsIsTransposed ? rhs.n_1 : rhs.n_2;
      int K = n_2;
      T * A = pArray;
      T * B = rhs.pArray;
      T * C = result.pArray;
      
      if(rhsIsTransposed)
	{
	  if(USE_KAHAN)
	    {
	      T Cee, Why, Tee;
	      int i, j, k;
	      
	      for (i = 0; i < M; ++i)
		{
		  const register T *Ai_ = A + i*K;
		  for (j = 0; j < N; ++j)
		    {
		      const register T *B_j = B + j*K;
		      register T cij = Ai_[0] * B_j[0];
		      Cee = 0;
		      for (k = 1; k < K; ++k)
			{
			  /*
			    This is the KSF. The parenthesis
			    are critical to its success.
			  */
			  Why = Ai_[k] * B_j[k] - Cee;
			  Tee = cij + Why;
			  Cee = (Tee - cij) - Why;
			  cij = Tee;
			}
		      C[i*N + j] = cij;
		    }
		}
	    }
	  else//no kahan summation formula
	    {
	      int i, j, k;
	      
	      for (i = 0; i < M; ++i)
		{
		  const register T *Ai_ = A + i*K;
		  for (j = 0; j < N; ++j)
		    {
		      const register T *B_j = B + j*K;
		      register T cij = 0;
		      for (k = 0; k < K; ++k)
			{
			  cij += Ai_[k] * B_j[k];
			}
		      C[i*N + j] = cij;
		    }
		}
	    }
	}
      else// rhs is not transposed, KSF not implemented here
	{
	  for (int i = 0; i < M; ++i)
	    {
	      register T *Ai_ = A + i*K;
	      for (int j = 0; j < N; ++j)
		{
		  register T cij = 0;
		  for (int k = 0; k < K; ++k)
		    {
		      cij += Ai_[k] * B[k*N+j];
		    }
		  C[i*N + j] = cij;
		}
	    }
	}
    }
  
  /**
     Calculates a dot product between all the elements in one
     Array2D and all the elements in another Array2D
  */
  T dotAllElectrons(const Array2D<T> & rhs)
    {
      register T temp = 0;
      for(int i=0; i<n_1*n_2; i++)
        temp += pArray[i]*rhs.pArray[i];
      return temp;
    }
  
  /**
     Calculates a dot product between all the elements in one
     row from one Array2D and one row in another Array2D
  */
  T dotOneElectron(const Array2D<T> & rhs, int whichElectron)
    {
      register T temp = 0;
      T * l = pArray + whichElectron*n_2;
      T * r = rhs.pArray + whichElectron*n_2;
      for(int i=0; i<n_1; i++)
        temp += l[i]*r[i];
      return temp;
    }

  void operator+=( const Array2D<T> & rhs)
    {
      for(int i=0; i<n_1*n_2; i++)
        pArray[i] += rhs.pArray[i];
    }

#endif
  
  Array2D<T> operator*(const Array2D<T> & rhs) const
    {
      Array2D<T> TEMP(n_1,rhs.n_2);
      TEMP = 0;
      gemm(rhs,TEMP,false);
      return TEMP;
    }
  
  void setToIdentity(int i)
    {
      for(int j=0; j<min(n_1,n_2); j++)
        {
          (*this)(i,j) = (T)(0.0);
          (*this)(j,i) = (T)(0.0);
        }
      (*this)(i,i) = (T)(1.0);
    }
  
  void setToIdentity()
    {
      *this = (T)(0.0);
      for(int i=0; i<min(n_1,n_2); i++)
        (*this)(i,i) = (T)(1.0);
    }
  
  bool isIdentity()
    {
      for(int i=0; i<n_1; i++)
        for(int j=0; j<n_2; j++)
          {
            T val       = (*this)(i,j);
            if(i == j) val -= 1.0;
            if(fabs(val) > REALLYTINY) return false;
          }
      return true;
    }
  
  /**
     Transpose the matrix.
  */
  void transpose()
    {
      Array2D<T> old = *this;
      if(n_1 != n_2)
	{
	  n_1 = old.dim2();
	  n_2 = old.dim1();
	}
      for(int i=0; i<n_1; i++)
        for(int j=0; j<n_2; j++)
	  (*this)(i,j) = old(j,i);
      old.deallocate();
    }
  
  /**
     A measurement of how far this matrix is from being symmetric.
     A returned value of 0.0 means it is perfectly symmetric.
  */
  double nonSymmetry()
    {
      if(n_1 != n_2)
        {
          cerr << "Array2D::symmetric just handles square matrices right now\n";
          return -1.0;
        }
      
      if(n_1 < 2)
        return 0.0;
      
      double diff = 0;
      for(int i=0; i<n_1; i++)
        for(int j=0; j<i; j++)
          {
            double term = 0;
            if(fabs((*this)(i,j)) > REALLYTINY)
              term = ((*this)(i,j) - (*this)(j,i))/(*this)(i,j);
            diff += term * term;
          }
      
      //relative error per pair
      diff = 2.0*diff/(n_1 * n_2 - n_1);
      
      //if the matrix is symmetric, it should return zero
      //otherwize, this measurement of symmetry.
      return diff;
    }
  
  /**
     Compare two Array2D objects to see if they're the same.
     @param reallyBad is the cutoff at which we'll print out the two values
     @param print if false, this function will not print anything.
     @return the relative error, averaged over all the comparisons.
  */
  template <class _RHS>
    double compare(const Array2D<_RHS> & rhs, double reallyBad, bool print, bool & same)
    {
      if(n_1 != rhs.dim1() || n_2 != rhs.dim2())
        {
          cerr << "ERROR: dims don't match for comparison: " << n_1 << "x"
	       << n_2 << " * " << rhs.dim2() << "x" << rhs.dim1() << endl;
          return 0.0;
        }
      
      same = true;
      int largestI = 0, largestJ = 0;
      double badCount = 0;
      double currentError=0, largestError = 0;
      double error = 0, error2 = 0, sample_size = 0;
      
      for(int i=0; i<n_1; i++)
        {
          for(int j=0; j<n_2; j++)
            {
              currentError = (double)(get(i,j)-rhs.get(i,j));
              double den = 1.0;
              if(fabs(rhs.get(i,j)) > 1e-200)
                den = rhs.get(i,j);
              currentError = fabs( currentError / den );
	      
              if(currentError > largestError)
                {
                  largestError = currentError;
                  largestI = i; largestJ = j;
                }
	      
              if(currentError > reallyBad)
                {
                  badCount++;
                  same = false;
                  if(print)
                    {
                      cout << scientific;
                      int prec = 8;
                      int width = 15;
                      cout << "Error at (" << i << "," << j << "): this "
			   << setprecision(prec) << setw(width) << get(i,j) << " rhs "
			   << setprecision(prec) << setw(width) << rhs.get(i,j) << " current error " << currentError << endl;
		      
                    }
                }
	      
              if(currentError < reallyBad)
                {
                  error  += currentError;
                  error2 += currentError * currentError;
                  sample_size++;
                }
            }
        }
      
      if(print)
        {
          cout << scientific;
          int prec = 8;
          int width = 15;
          cout << "Largest Error at (" << largestI << "," << largestJ << "): this "
	       << setprecision(prec) << setw(width) << get(largestI,largestJ) << " rhs "
	       << setprecision(prec) << setw(width) << rhs.get(largestI,largestJ) << " current error " << largestError << endl;
        }
      
      error  =  error/sample_size;
      error2 = error2/sample_size;
      double std_dev = sqrt(error2 - error*error);
      
      if(print)
        printf("sample_size %5.2e ave_rel_error %6.3e std_dev %6.3e (%4.2e worse_than %5.3e)\n",
               (double)sample_size, error, std_dev, badCount, reallyBad);
      
      //This cuttoff is very generous for double, but about right for float.
      //You probably shouldn't rely on "same" for anything important.
      if(fabs(error) > 1e-7) same = false;
      
      return error;
    }
  
  /**
     Checks whether two Array2D objects are the same, to
     within a certain tolerance, 1e-8.
  */
  bool operator==(const Array2D & rhs)
    {
      bool same;
      compare(rhs,1e-5,false,same);
      return same;
    }
  
  void operator=(const Array2D & rhs)
    {
      if(n_1 != rhs.dim1() || n_2 != rhs.dim2())
	allocate(rhs.dim1(),rhs.dim2());
      
      memcpy(pArray, rhs.pArray, sizeof(T)*n_1*n_2);
    }
  
  /**
     Sets two arrays equal. memcpy would probably break if T is a class. is this a problem?
     if you want to set to Array2D's containing a class equal, then use the (i,j) operator.
  */
  template <class _RHS>
    void operator=(const Array2D<_RHS> & rhs)
    {
      if(n_1 != rhs.dim1() || n_2 != rhs.dim2())
	allocate(rhs.dim1(),rhs.dim2());
      
      for(int i=0; i<n_1; i++)
	for(int j=0; j<n_2; j++)
	  pArray[map(i,j)] = (T)(rhs.get(i,j));
    }
  
  /**
     Sets all of the elements in an array equal to the same value. There's an obvious shortcut if
     the constant happens to be 0. in retrospect, while memset works well for doubles, floats, etc,
     using this operator would probably crash the program if T was a class.
  */
  void operator=(const T C)
    {
#ifdef QMC_DEBUG
      if(n_1 < 1 || n_2 < 1)
        {
          cerr << "Warning: using Array2D::operator=(const T C)"
	       << " with dims " << n_1 << " and " << n_2 << endl;
        }
#endif
      if(C == 0)
        {
          memset(pArray,0,sizeof(T)*n_1*n_2);
          return;
        }
      
      for(int i=0; i<n_1*n_2; i++)
        pArray[i] = C;
    }
  
  /**
     Returns the product of an array and a scalar.
  */
  Array2D operator*(const T C)
    {
#ifdef QMC_DEBUG
      if(n_1 < 1 || n_2 < 1)
        {
          cerr << "Warning: using Array2D::operator*(const T C)"
	       << " with dims " << n_1 << " and " << n_2 << endl;
        }
#endif
      Array2D<T> TEMP(n_1,n_2);
      for(int i=0; i<n_1*n_2; i++)
        TEMP.pArray[i] = C*pArray[i];
      return TEMP;
    }
  
  /**
     Returns the difference between two Array2Ds
  */
  Array2D<T> operator-(const Array2D<T> & rhs) const
    {
#ifdef QMC_DEBUG
      if(n_1 < 1 || n_2 < 1)
        {
          cerr << "Warning: using Array2D::operator-(const T C)"
	       << " with dims " << n_1 << " and " << n_2 << endl;
        }
#endif
      Array2D<T> TEMP(n_1,n_2);
      for(int i=0; i<n_1*n_2; i++)
	TEMP.pArray[i] = pArray[i] - rhs.pArray[i];
      return TEMP;
    }
  
  /**
     Returns the sum of two Array2Ds
  */
  Array2D<T> operator+(const Array2D<T> & rhs) const
    {
#ifdef QMC_DEBUG
      if(n_1 < 1 || n_2 < 1)
        {
          cerr << "Warning: using Array2D::operator+(const T C)"
	       << " with dims " << n_1 << " and " << n_2 << endl;
        }
#endif
      Array2D<T> TEMP(n_1,n_2);
      for(int i=0; i<n_1*n_2; i++)
	TEMP.pArray[i] = pArray[i] + rhs.pArray[i];
      return TEMP;
    }
  
  /**
     Sets this array equal to itself times a scalar value.
  */
  void operator*=(const T C)
    {
#ifdef QMC_DEBUG
      if(n_1 < 1 || n_2 < 1)
        {
          cerr << "Warning: using Array2D::operator*=(const T C)"
	       << " with dims " << n_1 << " and " << n_2 << endl;
        }
#endif
      for(int i=0; i<n_1*n_2; i++)
        pArray[i] *= C;
    }
  
  /**
     Sets this array equal to itself divided by a scalar value.
  */
  void operator/=(const T C)
    {
#ifdef QMC_DEBUG
      if(n_1 < 1 || n_2 < 1)
        {
          cerr << "Warning: using Array2D::operator/=(const T C)"
	       << " with dims " << n_1 << " and " << n_2 << endl;
        }
      if(C == 0.0)
	{
	  cerr << "Error: divide by zero in Array2D::operator/=" << endl;
	}
#endif
      
      T inv = 1.0/C;
      operator*=(inv);
    }
  
  /**
     Creates an array.
  */
  Array2D()
    {
      pArray = 0; n_1 = 0; n_2 = 0;
    }

  /**
     Creates an array and allocates memory.
     
     @param i size of the array's first dimension.
     @param j size of the array's second dimension.
  */
  Array2D(int i, int j)
    {
      pArray = 0;
      n_1 = 0; n_2 = 0;
      allocate(i,j);
    }
  
  /**
     Create a random matrix.
     It uses the UNIX rand function, which will
     produce a uniform random number, [0,1)
     
     @param shift will be added to the random value
     @param range after adding shift, we'll multiply by range
     @param iseed will be used to initialize srand, if iseed > 0
  */
  Array2D(int i, int j, double shift, double range, long & iseed)
    {
      pArray = 0;
      n_1 = 0; n_2 = 0;
      allocate(i,j);
      
      if(iseed > 0) srand(iseed);
      for(int idx=0; idx<n_1*n_2; idx++)
        {
          // between 0 and 1.0
          double val = (double)(rand())/RAND_MAX;
          val = range * (val + shift);
          pArray[idx] = (T)val;
        }
    }
  
  /**
     Creates an array and sets it equal to another array.
     
     @param rhs array to set this array equal to.
  */
  Array2D( const Array2D<T> & rhs)
    {
      n_1 = 0;
      n_2 = 0;
      pArray = 0;
      allocate(rhs.dim1(),rhs.dim2());
      operator=(rhs);
    }
  
  /**
     This function makes it easy to replace a row with an Array1D
  */
  void setRow(int whichRow, Array1D<T> & rhs)
    {
#ifdef QMC_DEBUG
      if(whichRow >= n_1)
	cerr << "Error: invalid row index in Array2D::setRow\n";
      if(rhs.dim1() > n_2)
	cerr << "Error: rhs length is incorrect in Array2D::setRow\n";
#endif
      memcpy(pArray + n_2*whichRow, rhs.array(), sizeof(T)*n_2);
    }
  
  /**
     This function takes advantage of the row-major format of the matricies
     to copy numRows worth of data from one Array2D to another
  */
  void setRows(int to, int from, int numRows, const Array2D<T> & rhs)
    {
#ifdef QMC_DEBUG
      if(to < 0 || from < 0 || numRows < 0 ||
	 (to+numRows) > n_1 || (from+numRows) > rhs.dim1() ||
	 n_2 != rhs.dim2())
	{
	  cerr << "Error: incompatible dimensions in setRows:\n"
	       << "      to = " << to << endl
	       << "    from = " << from << endl
	       << " numRows = " << numRows << endl
	       << "     n_1 = " << n_1 << endl
	       << "     n_2 = " << n_2 << endl
	       << " rhs.n_1 = " << rhs.dim1() << endl
	       << " rhs.n_2 = " << rhs.dim2() << endl;
	  assert(0);
	}
      assert(pArray);
      assert(rhs.pArray);
      if(numRows == 0)
	cout << "Warning: 0 rows in Array2D::setRows" << endl;
#endif
      memcpy(pArray + n_2*to, rhs.pArray + n_2*from, sizeof(T)*n_2*numRows);
    }

  /*
    Interchange two rows in our data.
  */
  void swapRows(int i, int j)
    {
      Array2D<T> temp = *this;
      setRows(i,j,1,temp);
      setRows(j,i,1,temp);
    }

  /**
     Copies a column from one Array2D to another
  */
  void setColumn(int to, int from, const Array2D<T> & rhs)
    {
      for(int i=0; i<n_1; i++)
	pArray[map(i,to)] = rhs.pArray[rhs.map(i,from)];
    }
  
  /**
     Destroy's the array and cleans up the memory.
  */
  ~Array2D()
    {
      deallocate();
    }
  
  /**
     This particular choice indicates row-major format.
  */
  inline int map(int i, int j) const
    {
#ifdef QMC_DEBUG
      if(i >= n_1 || j >= n_2 ||
	 i < 0 || j < 0)
	{
	  cerr << "Error: bad dimensions in Array2D::map" << endl
	       << "  i = " << i << endl
	       << "  j = " << j << endl
	       << "n_1 = " << n_1 << endl
	       << "n_2 = " << n_2 << endl;
	  assert(0);
	}
      if(pArray == 0)
	{
	  cerr << "Error: pArray == 0 in Array2D::map" << endl;
	  assert(0);
	}
#endif
      return n_2*i + j;
    }
  
  /*
    A pointer to a row.
  */
  T* operator[](int row)
    {
#ifdef QMC_DEBUG
      if(row >= n_1 || row < 0)
        {
          cerr << "Error: bad dimension in Array2D::operator[]" << endl
	       << " row = " << row << endl
	       << " n_1 = " << n_1 << endl;
          assert(0);
        }
#endif
      return & pArray[n_2 * row];
    }
  
  /**
     Accesses element <code>(i,j)</code> of the array.
  */
  T& operator()(int i,int j)
    {
      return pArray[map(i,j)];
    }
  
  /**
     Accesses element <code>(i,j)</code> of the array, except this prevents
     modification of that element.
  */
  T get(int i, int j) const
    {
      return pArray[map(i,j)];
    }
  
  double frobeniusNorm()
    {
      double sum = 0;
      for(int i=0; i<n_1*n_2; i++)
        sum += pArray[i]*pArray[i];
      return sqrt(sum);
    }
  
  double pInfNorm()
    {
      double norm = 0;
      double sum = 0;
      for(int i=0; i<n_1; i++)
        {
          for(int j=0; j<n_2; j++)
            {
              sum += fabs( get(i,j) );
            }
          norm = max(sum,norm);
          sum = 0;
        }
      return norm;
    }
  
  double p1Norm()
    {
      double norm = 0;
      double sum = 0;
      for(int j=0; j<n_2; j++)
        {
          for(int i=0; i<n_1; i++)
            {
              sum += fabs( get(i,j) );
            }
          norm = max(sum,norm);
          sum = 0;
        }
      return norm;
    }
  
  /**
     Try to fit as many columns as possible in a width approximately
     that of a terminal. Also, try to print as many significant figures
     as possible.
     
     These settings are very approximate.
  */
  friend ostream& operator<<(ostream & strm, const Array2D<T> & rhs)
    {
      if(rhs.n_2 < 8)
	rhs.printArray2D(strm, 15, 10, -1, ',', true);
      else if(rhs.n_2 < 11)
	rhs.printArray2D(strm, 8, 10, -1, ',', true);
      else
	rhs.printArray2D(strm, 2, 7, -1, ',', false);
      
      return strm;
    }
  
  /**
     This makes it easier to print out the Array2D
     for input into Mathematica. In Mathematica, you will
     still have to replace all the scientific notational "e" with "*^".
  */
  void printArray2D(ostream & strm, int numSigFig, int columnSpace, int maxJ, char sep, bool sci) const
    {
      strm.precision(4);
      if(maxJ < 0) maxJ = n_2;
      strm << "{";
      for(int i=0; i<n_1;i++)
        {
          if(i==0) strm <<  "{";
          else     strm << " {";
          for(int j=0; j<n_2 && j<maxJ; j++)
            {
              if(sci) strm.setf(ios::scientific);
	      else    strm.unsetf(ios::scientific);

              strm << setprecision(numSigFig)
		   << setw(numSigFig+columnSpace);
	      strm << pArray[map(i,j)];

              if(j<n_2-1) strm << sep;
            }
          if(i<n_1-1)  strm << "},";
          else         strm << "}}";

	  //label the rows...
	  //strm << " " << setw(4) << i;
          strm << endl;
        }
    }
  
  /**
     LU decomposition using the algorithm in numerical recipes for a dense matrix.
     
     
     @param a a \f$N \times N\f$ matrix which is destroyed during the operation.
     The resulting LU decompositon is placed here.
     @param indx a \f$N\f$ dimensional array which records the row permutation 
     from partial pivoting.
     @param d used to give det(a) the correct sign
     @param calcOK returns false if the calculation is singular and true otherwise
  */
  void ludcmp(double *d, bool *calcOK)
    {
      int i,j,k;
      int imax = -1;
      long double big,dum,temp;
      long double Why, Cee, Tee;
      register long double sum;
      long double one = (T)(1.0);
      long double zero = (T)(0.0);
      
      *calcOK = true;
      *d=1.0;
      
      diVV.allocate(n_1);
      
      //this section finds the largest value from each column
      //and puts its inverse in the vv vector.
      for (i=0;i<n_1;i++)
        {
          big=zero;
          for (j=0;j<n_1;j++)
            if ((temp=(T)(fabs((double)get(i,j)))) > big)
              big=temp;
	  
          //checks if the column is full of zeros
          if (big == zero)
            {
              // cerr << "Singular matrix in routine ludcmp*********" << endl;
              *calcOK = false;
              return;
            }
	  
          diVV(i)= (T)(one/big);
        }
      
      //loop over columns
      for (j=0;j<n_1;j++)
        {
	  
          //part 1: i < j
          for (i=0;i<j;i++)
            {
              sum=-1.0*get(i,j);
	      
              if(USE_KAHAN)
                {
                  Cee = 0;
                  for (k=0;k<i;k++)
                    {
                      Why = get(i,k)*get(k,j) - Cee;
                      Tee = sum + Why;
                      Cee = (Tee - sum) - Why;
                      sum = Tee;
                    }
                }
              else
                {
		  
                  for (k=0;k<i;k++)
                    sum += get(i,k)*get(k,j);
		  
                }
	      
              pArray[map(i,j)] = (T)(-1.0*sum);
            }
	  
          //part 2: i >= j
          big=zero;
          for (i=j;i<n_1;i++)
            {
              sum=-1.0*get(i,j);
	      
              if(USE_KAHAN)
                {
                  Cee = 0;
                  for (k=0;k<j;k++)
                    {
                      Why = get(i,k)*get(k,j) - Cee;
                      Tee = sum + Why;
                      Cee = (Tee - sum) - Why;
                      sum = Tee;
                    }
                }
              else
                {
		  
                  for (k=0;k<j;k++)
                    sum += get(i,k)*get(k,j);
		  
                }
	      
              pArray[map(i,j)] = (T)(-1.0*sum);
	      
              //find the best row to pivot
              if ( (dum=diVV(i)*(T)(fabs((double)sum))) >= big)
                {
                  big=dum;
                  imax=i;
                }
            }
	  
          //if our row isn't the best to pivot, we need to
          //change to the best
          if (j != imax)
            {
	      
              //swap rows
              for (k=0;k<n_1;k++)
                {
                  dum=get(imax,k);
                  pArray[map(imax,k)]=get(j,k);
                  pArray[map(j,k)] = (T)(dum);
                }
              //a.swapRows(imax,j);
	      
              //indicate permutation change
              *d = -(*d);
	      
              //update vv
              diVV(imax)=diVV(j);
            }
          diINDX(j)=imax;
	  
          //make sure we don't divide by zero
          if (get(j,j) == zero)
            pArray[map(j,j)]=(T)(REALLYTINY);
	  
          //each element needs to be divided by it's diagonal
          if (j != n_1)
            {
              dum=one/(get(j,j));
              for (i=j+1;i<n_1;i++)
                pArray[map(i,j)] *= (T)(dum);
            }
	  
        }
    }
  
  /**
     LU backsubstitution using the algorithm in numerical
     recipes for a dense matrix. 
     
     @param a the LU decomposition of a matrix produced by ludcmp
     @param indx a \f$N\f$ dimensional array which records the row permutation 
     from partial pivoting generated by ludcmp
     @param b the \f$N\f$ dimensional array right hand side of the system of 
     equations to solve
  */
  void lubksb(Array1D<T> &b)
    {
      int ii=-1, ip;
      long double sum;
      long double Cee, Why, Tee;
      
      for (int i=0;i<n_1;i++)
        {
          ip=diINDX(i);
          sum= -1.0*b(ip);
          b(ip)=b(i);
	  
          if (ii>=0)
            {
              //we do forward-subsitution to find "y"
	      
              if(USE_KAHAN)
                {
                  Cee = 0;
                  for (int j=ii;j<=i-1;j++)
                    {
                      Why = get(i,j)*b(j) - Cee;
                      Tee = sum + Why;
                      Cee = (Tee - sum) - Why;
                      sum = Tee;
                    }
                }
              else
                {
		  
                  for (int j=ii;j<=i-1;j++)
                    sum += get(i,j)*b(j);
                }
	      
            }
          else if (sum)
            {
              //our first non-zero element
              ii=i;
            }
	  
          b(i) = (T)(-1.0*sum);
        }
      
      for (int i=n_1-1;i>=0;i--)
        {
          sum=-1.0*b(i);
	  
          //back-substitution to find "x"
          if(USE_KAHAN)
            {
              Cee = 0;
              for (int j=i+1;j<n_1;j++)
                {
                  Why = get(i,j)*b(j) - Cee;
                  Tee = sum + Why;
                  Cee = (Tee - sum) - Why;
                  sum = Tee;
                }
            }
          else
            {
	      
              for (int j=i+1;j<n_1;j++)
                sum += get(i,j)*b(j);
            }
	  
          b(i) = (T)(-1.0*sum/get(i,i));
        }
    }
  
  /**
     This method is supposed to improve the numerical
     precision of the inverse. Call this only *after*
     you've already decomposed the matrix into L and U.
  */
  void mprove(const Array2D<T> & A, Array2D<T> & inv)
    {
      int i,j;
      Array1D<T> r = Array1D<T>(n_1);
      
      for(int k=0; k<n_1; k++)
        {
	  
          for(int i=0; i<n_2; i++)
            {
              long double sdp = i == k ? -1.0 : 0.0;
              for(int j=0; j<n_2; j++)
                sdp += (long double) A.get(i,j) * (long double) inv(k,j);
              r(i) = sdp;
            }
	  
          lubksb(r);
          for(int i=0; i<n_1; i++) inv(k,i) -= r(i);
        }
    }
  
  /**
     Calculates the inverse and determinant of a matrix using a 
     dense LU solver.  This 
     method scales as \f$O(1 N^3)\f$.
     
     @param a a \f$N \times N\f$ matrix 
     @param inv inverse of a is returned here
     @param det determinant of a is returned here
     @param calcOK returns false if the calculation is singular and true otherwise
  */
#if defined USELAPACK || defined USECLPK
  void determinant_and_inverse(Array2D<double> &inv, double& det, bool *calcOK)
    {
      diINDX.allocate(n_1);
      
      inv.setToIdentity();
      
      int info = 0;
      int N    = n_1;
      int NRHS = n_1;
      int LDA  = n_1;
      int LDB  = n_1;
      
      //Warning: this will destroy pArray
      FORTRAN_FUNC(dgesv)(&N, &NRHS, pArray, &LDA, diINDX.array(), inv.array(), &LDB, &info);
      
      if(info == 0) *calcOK = true;
      else *calcOK = false;
      
      det = 1.0;
      for(int i=0; i<n_1; i++)
	{
	  if(diINDX(i) != i+1)
	    det *= -1;
	  det *= get(i,i);
	}
    }
  
  void determinant_and_inverse(Array2D<float> &inv, double& det, bool *calcOK)
    {
      diINDX.allocate(n_1);
      
      inv.setToIdentity();
      
      int info = 0;
      int N    = n_1;
      int NRHS = n_1;
      int LDA  = n_1;
      int LDB  = n_1;
      
      //Warning: this will destroy pArray
      FORTRAN_FUNC(sgesv)(&N, &NRHS, pArray, &LDA, diINDX.array(), inv.array(), &LDB, &info);
      
      if(info == 0) *calcOK = true;
      else *calcOK = false;

      det = 1.0;
      for(int i=0; i<n_1; i++)
	{
	  if(diINDX(i) != i+1)
	    det *= -1;
	  det *= get(i,i);
	}
    }
#else
  /**
     this is O(1*N^3)
     
     I originally had several of the arrays here as function-static variables, but
     this was causing the program to randomly crash in other parts of the code. I
     don't know why.
  */
  void determinant_and_inverse(Array2D<T> &inv, double& det, bool *calcOK)
    {
      double d;
      T one = (T)1.0;
      T zero = (T)0.0;
      
      diINDX.allocate(n_1);
      diCol.allocate(n_1);
      inv.allocate(n_1,n_1);
      
      ludcmp(&d,calcOK);
      
      if( *calcOK )
        {
          for(int j=0; j<n_1; j++)
            {
              diCol = zero;
              diCol(j) = one;
              lubksb(diCol);
              for(int i=0; i<n_1; i++)
                inv(i,j) = diCol(i);
            }
        }
      
      for(int i=0; i<n_1; i++)
        d *= get(i,i);
      det = d;
    }
#endif
  
  /*
    Call this function on D^(-1) to update it given
    the new row in D.
    
    The return value will be |D new| / |D|.
    
    This is slightly different from the Sherman-Morrison
    formula because here we are replacing a row, as opposed
    to adding values to the row.
    
    @param row which we are updating the inverse with
    @param newD the matrix which has the new row
    @param rowInD which row in newD to use for the update
    @return the ratio of new to old determinants
  */
  template <class _RHS>
    double inverseUpdateOneRow(int row, Array2D<_RHS> & newD, int rowInD)
    {
      int d = n_1;
#ifdef QMC_DEBUG
      if( rowInD >= newD.dim1() || row >= n_1 ||
	  rowInD < 0 || row < 0 ||
	  n_2 != newD.dim2())
        {
          cout << "Error: dimensions don't match in inverse update." << endl;
          cout << "   n_1 = " << n_1 << endl;
          cout << "   n_2 = " << n_2 << endl;
          cout << "   row = " << row << endl;
          cout << " D.n_1 = " << newD.dim1() << endl;
          cout << " D.n_2 = " << newD.dim2() << endl;
          cout << " D row = " << rowInD << endl;
          assert(0);
        }
#endif
      double detRatio = 0;
      for(int r=0; r<d; r++)
        detRatio += (T)(newD.get(rowInD,r)) * pArray[map(r,row)];
      double q = 1.0/detRatio;
      
      for(int c=0; c<d; c++)
        {
          //We need to preserve this column for the moment
          if(c == row)
            continue;
	  
          double val = 0.0;
          for(int l=0; l<d; l++)
            val += (T)(newD.get(rowInD,l)) * pArray[map(l,c)];
          val *= q;
	  
          for(int r=0; r<d; r++)
            pArray[map(r,c)] -= val * pArray[map(r,row)];
        }
      
      //Lastly, we update the column
      for(int r=0; r<d; r++)
        pArray[map(r,row)] *= q;
      
      return detRatio;
    }

  /*
    Call this function on D^(-1) to update it given
    the new column in D.
    
    The return value will be |D new| / |D|.
    
    This is slightly different from the Sherman-Morrison
    formula because here we are replacing a column, as opposed
    to adding values to the column.
    
    @param column which we are updating the inverse with
    @param newD the matrix which has the new column
    @param columnInD which column in newD to use for the update
    @return the ratio of new to old determinants
  */
  template <class _RHS>
    double inverseUpdateOneColumn(int column, Array2D<_RHS> & newD, int columnInD)
    {
      int d = n_2;
#ifdef QMC_DEBUG
      if( columnInD >= newD.dim2() || column >= n_2 ||
	  columnInD < 0 || column < 0 ||
	  n_1 != newD.dim1())
        {
          cout << "Error: dimensions don't match in inverse update." << endl;
          cout << "     n_1 = " << n_1 << endl;
          cout << "     n_2 = " << n_2 << endl;
          cout << "  column = " << column << endl;
          cout << "   D.n_1 = " << newD.dim1() << endl;
          cout << "   D.n_2 = " << newD.dim2() << endl;
          cout << "D column = " << columnInD << endl;
          assert(0);
        }
#endif
      double detRatio = 0;
      for(int c=0; c<d; c++)
        detRatio += (T)(newD.get(c,columnInD)) * pArray[map(column,c)];
      double q = 1.0/detRatio;
      
      for(int r=0; r<d; r++)
        {
          //We need to preserve this row for the moment
          if(r == column)
            continue;
	  
          double val = 0.0;
          for(int l=0; l<d; l++)
            val += (T)(newD.get(l,columnInD)) * pArray[map(r,l)];
          val *= q;
	  
          for(int c=0; c<d; c++)
            pArray[map(r,c)] -= val * pArray[map(column,c)];
        }
      
      //Lastly, we update the row
      for(int c=0; c<d; c++)
        pArray[map(column,c)] *= q;
      
      return detRatio;
    }
  
  /**
     Interpret this Array2D has a system of equations,
     and by using elementary row operations, convert it
     to its row reduced echelon form.

     If rank(n_1) = n_2, then this should produce the identity
     matrix.

     If rank(n_1) > n_2, then you'll get row(s) of zeros since the system
     is overdetermined.
     
     If rank(n_1) < n_2, then you'll get an identity matrix the size of
     rank(n_1), and columns containing data for the remaining unconstrained
     variables.

     Note: this function does not have the same improvements that rref(int) has.
  */
  void rref()
    {
      double max_data, temp;
      int max_index, i, j, k, l;
      int order = min(n_1,n_2);
      Array1D<int> pivots(order);
      double tolerance = 1e-20;
      bool allowPivot = true;
      
      l = 0;
      for (k = 0; k <  order; k++)
	{
	  max_data  = pArray[map(l,k)];
	  max_index = l;
	  
	  if(allowPivot)
	    {
	      //try to figure out if there is a target row
	      //we'd rather operate on
	      for (i = l; i < n_1; i++)
		{
		  if (fabs(pArray[map(i,k)]) > fabs(max_data))
		    {
		      max_data  = pArray[map(i,k)];
		      max_index = i;
		    }
		}
	    }
	  
	  pivots(k) = max_index;
	  if (fabs(max_data) > tolerance)
	    {
	      //whether or not we pivot, we still want to scale
	      //our target row so that the first element is 1.0
	      for (j = k; j < n_2; j++)
		{
		  temp                     = pArray[map(max_index,j)] / max_data;
		  pArray[map(max_index,j)] = pArray[map(l,j)];
		  pArray[map(l,j)]         = temp;
		}
	      
	      //and we will subtract from all lower rows
	      //the target row, multiplied by a scale factor to give
	      //each row a 0.0 for this column.
	      for (i = l + 1; i < n_1; i++)
		{
		  double scale = pArray[map(i,l)];
		  for (j = k; j < n_2; j++)
		    pArray[map(i,j)] -= pArray[map(l,j)] * scale;
		}
	      
	      l++;
	    }
	}
      
      /*
	At this point, we have 1.0 along the diagonal, and 0.0 below
	the diagonal. Now we need to make the other
	half of that triangle to be zero.
      */
      for (k = 0; k < order-1; k++)
	{
	  if (pArray[map(k+1,k+1)] != 0)
	    {
	      for (i = k; i >= 0; i--) 
		{
		  double scale = pArray[map(i,k+1)];		
		  for (j = k+1; j < n_2; j++) 
		    pArray[map(i,j)] -= pArray[map(k+1,j)] * scale;
		}
	    }
	}
    }

  /**
     Interpret this Array2D as a system of equations,
     and solve for one variable (column index). This will
     turn the column into all 0.0 and one 1.0, causing
     this variable to be entirely dependend upon the other
     variables.

     We can eliminate at most as many independent variables as the rank
     of the matrix, rank(n_1) <= n_1. If you solve
     for more, then the code will arbitrarily choose others to "unsolve".

     Note, this will replace column t with 0.0s and one 1.0.
     The position of the 1.0 is arbitrary.

     @param t the variable to make completely dependent
  */
  void rref(int t)
    {
      if(isDependent(t) != -1) 
	{
	  //this variable is already dependent, so there's nothing to do
	  return;
	}

      /*
	Nothing will be damaged by bad calls, but it might be annoying to print
	messages.
      */
      bool printWarnings = true;
      if(t < 0 || t >= n_2)
	{
	  if(printWarnings)
	    {
	      cerr << "Error: invalid variable elimination. ";
	      cerr << "Number variables = " << n_2 << " but ";
	      cerr << "t = " << t << "." << endl;
	    }
	  return;
	}

      /*
	We just have to make sure that given enough calls
	to this function, all rows will be visited.
      */
      static int s = -1;
      s++;
      if(s >= n_1) s=0;

      double tolerance = 1e-10;

      /*
	We're going to divide by this scale, so
	we need to find a non zero value.
      */
      double scale = pArray[map(s,t)];
      int best = s;
      bool willUnsolve = false;
      for(int r=n_1-1; r>=0; r--)
	{
	  scale  = pArray[map(r,t)];
	  
	  if(fabs(scale) < tolerance)
	    continue;
	  else
	    best = r;

	  willUnsolve = constraintUsed(r) == -1 ? false : true;
	  if(constraintUsed(r) == -1)
	    {
	      best = r;
	      break;
	    }
	}

      if(best != s)
	swapRows(best,s);
      scale = pArray[map(s,t)];

      /*
	If we're given a column with all zeros, then there's nothing we
	can do for this variable.
      */
      if(fabs(scale) <= tolerance)
	{
	  if(printWarnings)
	    {
	      cerr << "Warning: none of the constraints use variable t = " << t
		   << " so we can't eliminate it!" << endl;
	    }
	  return;
	}

      /*
	Scale the constraint row so that the variable column has
	a 1.0 on that row.
	
	If our scale is already 1.0, don't waste time.
       */
      if(fabs(scale - 1.0) > tolerance)
	{
	  for(int c = 0; c < n_2; c++)
	    pArray[map(s,c)] = pArray[map(s,c)] / scale;
	}
      /*
	Mathematically, this is already 1.0, so we set it
	to eliminate roundoff error.
      */
      pArray[map(s,t)] = 1.0;

      /*
	Subtract from all other rows the target row,
	multiplied by the scale factor which will leave 0.0
	in this column.
      */
      for(int r = 0; r < n_1; r++) 
	{
	  scale = pArray[map(r,t)];

	  if(r == s || fabs(scale) < tolerance) continue;

	  for(int c = 0; c < n_2; c++) 
	    {
	      /*
		This subtraction is highly vulnerable to
		"catastrophic cancellation".
	      */
	      pArray[map(r,c)] -= pArray[map(s,c)] * scale;


	      //So eliminate values that are near to zero
	      if(fabs(pArray[map(r,c)]) < tolerance)
		pArray[map(r,c)] = 0.0;
	    }	  

	  //Mathematically, this is already 0.0
	  pArray[map(r,t)] = 0.0;
	}
    }

  /*
    This function is useful in conjunction with rref(int). It will
    count the number of dependent variables in the data.

    Of course, you can not have more dependent variables than
    rows.
  */
  int numDependent()
    {
      int numDep = 0;
      for(int j=0; j<dim2(); j++)
	if(isDependent(j) > -1)
	  numDep++;
      return numDep;
    }

  /*
    This function is useful in conjunction with rref(int).
    It will determine whether a particular column represents
    a dependent variable. If it does, then it will return
    the row index containing the 1.0.
  */
  int isDependent(int column)
    {
      double tolerance = 1e-50;
      int hasOne  = -1;
      bool zeros  = true;
      for(int i=0; i<dim1(); i++)
	{
	  if(fabs(pArray[map(i,column)] - 1.0) < tolerance && hasOne < 0)
	    hasOne = i;
	  else if(fabs(pArray[map(i,column)]) > tolerance)
	    zeros = false;
	}
      
      /*
	We know the variable is dependent if there is exactly
	one 1.0 in the column and the rest of the column is zero.
      */
      if(hasOne > -1 && zeros)
	return hasOne;

      return -1;
    }

  /**
     This function is useful in conjunction with rref(int).

     It will determine whether a particular row (interpreted
     as a constraint) has been used to make a variable dependent.
     If it has, then it will return the index of the dependent
     variable that uses this constraint.
  */
  int constraintUsed(int r)
    {
      for(int c=0; c<n_2; c++)
	if(fabs(pArray[map(r,c)] - 1.0) < 1e-10)
	  if(isDependent(c) == r)
	    return c;
      return -1;
    }

  /**
     Solves:
     A x = lambda B x
     
     Where lambda is the eigenvalue that corresponds to
     the (right-side) eigenvector x.
     
     The eigenvectors are stored as rows in this->pArray
     
     Notice: A and B are over written.
  */
  void generalizedEigenvectors(Array2D<double> A,
			       Array2D<double> B,
			       Array1D<Complex> & eigval)
    {
#ifdef USELAPACK
      allocate(A.n_1,A.n_1);
      
      char vl  = 'N';
      char vr  = 'V';
      int n    = n_1;
      int lda  = A.n_1;
      int ldb  = B.n_1;
      int ldvl = 1;
      int ldvr = n_1;
      int info = 0;
      
      Array1D<double> alphar(n);
      Array1D<double> alphai(n);
      Array1D<double> beta(n);
      
      Array1D<double> work(8*n);
      int lwork = work.dim1();
      FORTRAN_FUNC(dggev)(&vl,&vr,&n,
             A.array(), &lda, B.array(), &ldb,
             alphar.array(), alphai.array(), beta.array(),
             (double*)0,&ldvl,
             pArray, &ldvr,
             work.array(), &lwork,
             &info);
      
      if(lwork == -1 && info == 0)
        {
          cout << "Notice: dggev recommends " << work(0) << " memory." << endl;
        }
      else if(info == 0)
        {
          eigval.allocate(n);
          for(int i=0; i<n; i++)
            {
              //if(fabs(beta(i)) < 1e-50)
              //  cout << "Eigenvalue " << setw(3) << i << " has small beta " << setw(15) << beta(i) << endl;
              double real = alphar(i) / beta(i);
              double imag = alphai(i) / beta(i);
	      
              eigval(i) = Complex(real,imag);
            }
        }
      else
        {
          cout << "Error: dggev reported info = " << info << endl;
        }
      
      alphar.deallocate();
      alphai.deallocate();
      beta.deallocate();
      work.deallocate();
      
#else
      cout << "Error: dggev requires you to compile with LAPACK!!" << endl;
#endif
    }
  
  /**
     Find the eigenvalues for a real, symmetric matrix.
     
     I found this routine at:
     http://www3.telus.net/thothworks/EigRSvaloCPP0.html
  */
  Array1D<double> eigenvaluesRS()
    {
      Array1D<double> wr;
      
      if(n_1 != n_2)
        {
          cout << "Error: no eigenvalues for non square matrix!" << endl
	       << " n_1 = " << n_1 << endl
	       << " n_2 = " << n_2 << endl
	       << endl;
          return wr;
        }
      
      double check = nonSymmetry();
      if(check > 1e-5)
        {
          cout << "Error: eigenvaluesRS only works with symmetric matrices" << endl
	       << " nonSymmetry = " << check << endl;
          return wr;
        }
      
      wr.allocate(n_1);
      wr = 0.0;
      
      Array1D<double> fv1;
      fv1.allocate(n_1);
      
      Array2D<double> A = *this;
      
      int i, ii, ierr = -1, j, k, l, l1, l2;
      double c, c2, c3, dl1, el1, f, g, h, p, r, s, s2, scale, tst1, tst2;
      // ======BEGINNING OF TRED1 ===================================
      
      ii = n_1 - 1;
      for (i = 0; i < ii; i++)
        {
          wr[i] = A[ii][i];
          A[ii][i] = A[i][i];
        }//End for i
      wr[ii] = A[ii][ii]; // Take last assignment out of loop
      
      for (i = (n_1 - 1); i >= 0; i--)
        {
	  
          l = i - 1;
          scale = h = 0.0;
	  
          if (l < 0)
            {
              fv1[i] = 0.0;
              continue;
            } // End if (l < 0)
	  
          for (j = 0; j <= l; j++)
            scale += fabs(wr[j]);
	  
          if (scale == 0.0)
            {
              for (j = 0; j <= l; j++)
                {
                  wr[j] = A[l][j];
                  A[l][j] = A[i][j];
                  A[i][j] = 0.0;
                }//End for j
	      
              fv1[i] = 0.0;
              continue;
            } // End if (scale == 0.0)
	  
          for (j = 0; j <= l; j++)
            {
              wr[j] /= scale;
              h += wr[j]*wr[j];
            }//End for j
	  
          f = wr[l];
          g = ((f >= 0) ? -sqrt(h) : sqrt(h));
          fv1[i] = g*scale;
          h -= f*g;
          wr[l] = f - g;
	  
          if (l != 0)
            {
	      
              for (j = 0; j <= l; j++)
                fv1[j] = 0.0;
	      
              for (j = 0; j <= l; j++)
                {
                  f = wr[j];
                  g = fv1[j] + f*A[j][j];
                  for (ii = (j + 1); ii <= l; ii++)
                    {
                      g += wr[ii]*A[ii][j];
                      fv1[ii] += f*A[ii][j];
                    } // End for ii
                  fv1[j] = g;
                }//End for j
	      
              // Form p
	      
              f = 0.0;
              for (j = 0; j <= l; j++)
                {
                  fv1[j] /= h;
                  f += fv1[j]*wr[j];
                }//End for j
	      
              h = f/(h*2);
	      
              // Form q
	      
              for (j = 0; j <= l; j++)
                fv1[j] -= h*wr[j];
	      
              // Form reduced A
	      
              for (j = 0; j <= l; j++)
                {
                  f = wr[j];
                  g = fv1[j];
		  
                  for (ii = j; ii <= l; ii++)
                    A[ii][j] -= f*fv1[ii] + g*wr[ii];
		  
                }//End for j
	      
            } // End if (l != 0)
	  
          for (j = 0; j <= l; j++)
            {
              f = wr[j];
              wr[j] = A[l][j];
              A[l][j] = A[i][j];
              A[i][j] = f*scale;
            }//End for j
	  
        }//End for i
      
      // ======END OF TRED1 =========================================
      
      // ======BEGINNING OF TQL1 ===================================
      
      for (i = 1; i < n_1; i++)
        fv1[i - 1] = fv1[i];
      
      fv1[n_1 - 1] = tst1 = f = 0.0;
      
      for (l = 0; l < n_1; l++)
        {
	  
          j = 0;
          h = fabs(wr[l]) + fabs(fv1[l]);
	  
          tst1 = ((tst1 < h) ? h : tst1);
	  
          // Look for small sub-diagonal element
	  
          for (k = l; k < n_1; k++)
            {
              tst2 = tst1 + fabs(fv1[k]);
              if (tst2 == tst1) break; // fv1[n_1-1] is always 0, so there is no exit through the bottom of the loop
            }//End for k
	  
          if (k != l)
            {
	      
              do
                {
		  
                  if (j == 30)
                    {
                      ierr = l;
                      break;
                    } // End if (j == 30)
		  
                  j++;
		  
                  // Form shift
		  
                  l1 = l + 1;
                  l2 = l1 + 1;
                  g = wr[l];
                  p = (wr[l1] - g)/(2.0*fv1[l]);
                  r = pythag(p, 1.0);
                  scale = ((p >= 0) ? r : -r); // Use scale as a dummy variable
                  scale += p;
                  wr[l] = fv1[l]/scale;
                  dl1 = wr[l1] = fv1[l]*scale;
                  h = g - wr[l];
		  
                  for (i = l2; i < n_1; i++)
                    wr[i] -= h;
		  
                  f += h;
		  
                  // q1 transformation
		  
                  p = wr[k];
                  c2 = c = 1.0;
                  el1 = fv1[l1];
                  s = 0.0;
                  c3 = 99.9;
                  s2 = 99.9;
                  // Look for i = k - 1 until l in steps of -1
		  
                  for (i = (k - 1); i >= l; i--)
                    {
                      c3 = c2;
                      c2 = c;
                      s2 = s;
                      g = c*fv1[i];
                      h = c*p;
                      r = pythag(p, fv1[i]);
                      fv1[i + 1] = s*r;
                      s = fv1[i]/r;
                      c = p/r;
                      p = c*wr[i] - s*g;
                      wr[i + 1] = h + s*(c*g + s*wr[i]);
                    }//End for i
		  
                  p = -s*s2*c3*el1*fv1[l]/dl1;
                  fv1[l] = s*p;
                  wr[l] = c*p;
                  tst2 = tst1 + fabs(fv1[l]);
                }
              while (tst2 > tst1); // End do-while loop
	      
            } // End if (k != l)
	  
          if (ierr >= 0) //This check required to ensure we break out of for loop too, not just do-while loop
            break;
	  
          p = wr[l] + f;
	  
          // Order eigenvalues
	  
          // For i = l to 1, in steps of -1
          for (i = l; i >= 1; i--)
            {
              if (p >= wr[i - 1])
                break;
              wr[i] = wr[i - 1];
            }//End for i
	  
          wr[i] = p;
	  
        }//End for l
      
      // ======END OF TQL1 =========================================
      
      A.deallocate();
      fv1.deallocate();
      return wr;
    }
};

#endif
