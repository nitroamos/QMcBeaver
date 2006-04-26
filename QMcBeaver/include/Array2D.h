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
#include "cppblas.h"
#include <assert.h>
#include "Array1D.h"

//Array2D is included in all the classes where this change is relevant
#if defined SINGLEPRECISION || defined QMC_GPU
typedef float  qmcfloat;
const static  float REALLYTINY = 1e-35f;
#else
typedef double qmcfloat;
const static double REALLYTINY = 1e-250;
#endif

#if defined SINGLEPRECISION || defined QMC_GPU
#else
#endif

static const bool USE_KAHAN = false;

using namespace std;

/**
   A 2-dimensional template for making arrays.  All of the memory allocation
   and deallocation details are dealt with by the class. Some of the operators
   require the Array2D to contain doubles or floats.
   
   NOTICE: All Array2D a operator * Array2D b methods now require b to be transposed.
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
    int dim1()
    {
      return n_1;
    }

    /**
    Gets the number of elements in the array's second dimension.

    @return number of elements in the array's second dimension.
    */
    int dim2()
    {
      return n_2;
    }

    /**
    Gets the total number of elements in the array.

    @return total number of elements in the array.
    */
    int size()
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
      delete [] pArray;
      pArray = 0;

      n_1 = 0;
      n_2 = 0;
    }

    /**
       The flag USEATLAS dictates whether Array2D will use the library to speed it's math or not.
       Actually, if the dimensions of the data are small, using ATLAS may actually be slower.
       FYI: DGEMM stands for Double-precision GEneral Matrix-Matrix multiplication
       DGEMM convention: MxN = MxK * KxN; lda, ldb, and ldc are the n_2 of their respective Array2Ds
    */

    void setupMatrixMultiply(const Array2D<T> & rhs, Array2D<T> & result,
                             bool rhsIsTransposed)
    {
      if(rhsIsTransposed)
        {
          if(n_2 != rhs.n_2)
            {
              cerr << "ERROR: Transposed Matrix multiplication: " << n_1 << "x"
              << n_2 << " * " << rhs.n_2 << "x" << rhs.n_1 << endl;
              exit(1);
            }
          result.allocate(n_1,rhs.n_1);
        }
      else
        {
          if(n_2 != rhs.n_1)
            {
              cerr << "ERROR: Matrix multiplication: " << n_1 << "x"
              << n_2 << " * " << rhs.n_1 << "x" << rhs.n_2 << endl;
              exit(1);
            }
          result.allocate(n_1,rhs.n_2);
        }

      result = (T)(0.0);
    }

#ifdef USEATLAS

    void gemm(const Array2D<double> & rhs, Array2D<double> & result,
              bool rhsIsTransposed)
    {
      setupMatrixMultiply(rhs,result,rhsIsTransposed);
      CBLAS_TRANSPOSE myTrans = rhsIsTransposed ? CblasTrans : CblasNoTrans;
      cblas_dgemm(CBLAS_ORDER(CblasRowMajor),
                  CBLAS_TRANSPOSE(CblasNoTrans),myTrans,
                  n_1, rhs.n_1, n_2,
                  1.0, pArray, n_2,
                  rhs.pArray, rhs.n_2,
                  0.0, result.pArray, result.n_2);
    }

    void gemm(const Array2D<float> & rhs, Array2D<float> & result,
              bool rhsIsTransposed)
    {
      setupMatrixMultiply(rhs,result,rhsIsTransposed);
      CBLAS_TRANSPOSE myTrans = rhsIsTransposed ? CblasTrans : CblasNoTrans;
      cblas_sgemm(CBLAS_ORDER(CblasRowMajor),
                  CBLAS_TRANSPOSE(CblasNoTrans),myTrans,
                  n_1, rhs.n_1, n_2,
                  1.0, pArray, n_2,
                  rhs.pArray, rhs.n_2,
                  0.0, result.pArray, result.n_2);
    }

    /**
    This uses ATLAS to calculate a dot product between all the elements in one
    Array2D and all the elements in another Array2D
    */
    double dotAllElectrons(const Array2D<double> & rhs)
    {
      return cblas_ddot(n_1*n_2, pArray, 1, rhs.pArray, 1);
    }

    /**
    This uses ATLAS to calculate a dot product between all the elements in one
    row from one Array2D and one row in another Array2D
    */
    double dotOneElectron(const Array2D<double> & rhs, int whichElectron)
    {
      return cblas_ddot(n_1, pArray + whichElectron*n_2, 1, rhs.pArray + whichElectron*n_2, 1);
    }

    float dotAllElectrons(const Array2D<float> & rhs)
    {
      return cblas_sdot(n_1*n_2, pArray, 1, rhs.pArray, 1);
    }

    float dotOneElectron(const Array2D<float> & rhs, int whichElectron)
    {
      return cblas_sdot(n_1, pArray + whichElectron*n_2, 1, rhs.pArray + whichElectron*n_2, 1);
    }

#else

    void gemm(const Array2D<T> & rhs, Array2D<T> & result,
              bool rhsIsTransposed)
    {
      setupMatrixMultiply(rhs,result,rhsIsTransposed);

      int M = n_1;
      int N = rhs.n_1;
      int K = n_2;
      T * A = pArray;
      T * B = rhs.pArray;
      T * C = result.pArray;
      int i, j, k;

      if(rhsIsTransposed)
        {
          if(USE_KAHAN)
            {
              T Cee;
              T Why;
              T Tee;

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
#endif

    Array2D<T> operator*(const Array2D<T> & rhs)
    {
      Array2D<T> TEMP;
      TEMP = 0;
      gemm(rhs,TEMP,false);
      return TEMP;
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
    void ludcmp(int *indx, double *d, bool *calcOK)
    {
      int i,j,k;
      int imax = -1;
      T big,dum,temp;
      register T sum;
      T *vv = new T[n_1];
      T one = (T)(1.0);
      T zero = (T)(0.0);

      *calcOK = true;
      *d=1.0;

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

          vv[i]=one/big;
        }

      //loop over columns
      for (j=0;j<n_1;j++)
        {

          //part 1: i < j
          for (i=0;i<j;i++)
            {
              sum=get(i,j);
              for (k=0;k<i;k++)
                sum -= get(i,k)*get(k,j);
              pArray[map(i,j)]=sum;
            }

          //part 2: i >= j
          big=zero;
          for (i=j;i<n_1;i++)
            {
              sum=get(i,j);
              for (k=0;k<j;k++)
                sum -= get(i,k)*get(k,j);
              pArray[map(i,j)]=sum;

              //find the best row to pivot
              if ( (dum=vv[i]*(T)(fabs((double)sum))) >= big)
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
                  pArray[map(j,k)]=dum;
                }
              //a.swapRows(imax,j);

              //indicate permutation change
              *d = -(*d);

              //update vv
              vv[imax]=vv[j];
            }
          indx[j]=imax;

          //make sure we don't divide by zero
          if (get(j,j) == zero)
            pArray[map(j,j)]=(T)(REALLYTINY);

          //each element needs to be divided by it's diagonal
          if (j != n_1)
            {
              dum=one/(get(j,j));
              for (i=j+1;i<n_1;i++)
                pArray[map(i,j)] *= dum;
            }

        }
      delete [] vv;
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

    void lubksb(int *indx, Array1D<T> &b)
    {
      int ii=-1, ip;
      T sum;

      for (int i=0;i<n_1;i++)
        {
          ip=indx[i];
          sum=b(ip);
          b(ip)=b(i);

          if (ii>=0)   // in previous version had if(ii<0) Checked with NR. ck
            for (int j=ii;j<=i-1;j++)
              sum -= get(i,j)*b(j);//we do forward-subsitution to find "y"
          else if (sum)//our first non-zero element
            ii=i;

          b(i)=sum;
        }

      for (int i=n_1-1;i>=0;i--)
        {
          sum=b(i);
          for (int j=i+1;j<n_1;j++)
            sum -= get(i,j)*b(j);//back-substitution to find "x"
          b(i)=sum/get(i,i);
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
#if defined USELAPACK
    void determinant_and_inverse(Array2D<double> &inv, double& det, bool *calcOK)
    {
      Array1D<int> pivots = Array1D<int>(n_1);
      inv = 0;
      for(int i=0; i<n_1; i++)
        inv(i,i) = 1.0;

      /* int clapack_dgesv(const enum CBLAS_ORDER Order, const int N, const int NRHS,
                        double *A, const int lda, int *ipiv,
                        double *B, const int ldb);*/
      int info = clapack_dgesv(CBLAS_ORDER(CblasRowMajor),n_1,n_1,
                               pArray,n_1,pivots.array(),inv.array(),n_1);

      if(info == 0) *calcOK = true;
      else *calcOK = false;

      det = 1.0;
      for(int i=0; i<n_1; i++)
        det *= get(i,i);

      pivots.deallocate();
    }

    void determinant_and_inverse(Array2D<float> &inv, double& det, bool *calcOK)
    {
      Array1D<int> pivots = Array1D<int>(n_1);
      inv = 0;
      for(int i=0; i<n_1; i++)
        inv(i,i) = 1.0;

      /* int clapack_dgesv(const enum CBLAS_ORDER Order, const int N, const int NRHS,
                        double *A, const int lda, int *ipiv,
                        double *B, const int ldb);*/
      int info = clapack_sgesv(CBLAS_ORDER(CblasRowMajor),n_1,n_1,
                               pArray,n_1,pivots.array(),inv.array(),n_1);

      if(info == 0) *calcOK = true;
      else *calcOK = false;

      det = 1.0;
      for(int i=0; i<n_1; i++)
        det *= get(i,i);

      pivots.deallocate();
    }
#else
    /**
    this is O(1*N^3)
    WARNING!!!
    this has been modified and returns the transpose of the inverse, not the inverse
    this makes calculate_DerivativeRatios in QMCSlater.cpp look nicer

    I originally had several of the arrays here as function-static variables, but
    this was causing the program to randomly crash in other parts of the code. I
    don't know why.
    */
    void determinant_and_inverse(Array2D<T> &inv, double& det, bool *calcOK)
    {
      int n=dim1();
      double d;
      int *INDX = new int[n];
      Array1D<T> col(n);
      T one = (T)1.0;
      T zero = (T)0.0;

      if(n > col.dim1()) col.allocate(n);

      ludcmp(INDX,&d,calcOK);

      inv.allocate(n,n);

      if( *calcOK )
        {
          for(int j=0;j<n;j++)
            {
              for(int i=0;i<n;i++) col(i) = zero;
              col(j) = one;
              lubksb(INDX,col);
              //lubksbForInverse(a,INDX,col,j);
              //for(int i=0;i<n;i++) inv(i,j) = col(i);
              inv.setRow(j,col);
              //for(int i=0;i<n;i++) inv(j,i) = col(i);
            }
        }

      for(int i=0; i<n; i++) d*= get(i,i);
      det = d;

      delete [] INDX;
      col.deallocate();
    }
#endif

    /**
    Sets two arrays equal. memcpy would probably break if T is a class. is this a problem?
    if you want to set to Array2D's containing a class equal, then use the (i,j) operator.
    */
    void operator=(const Array2D & rhs)
    {
      if(n_1 != rhs.n_1 || n_2 != rhs.n_2) allocate(rhs.n_1,rhs.n_2);
      memcpy(pArray, rhs.pArray, sizeof(T)*n_1*n_2);
    }

    /**
    Sets all of the elements in an array equal to the same value. There's an obvious shortcut if
    the constant happens to be 0. in retrospect, while memset works well for doubles, floats, etc,
    using this operator would probably crash the program if T was a class.
    */
    void operator=(const T C)
    {
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
      Array2D<T> TEMP(n_1,n_2);
      for(int i=0; i<n_1*n_2; i++)
        TEMP.pArray[i] = C*pArray[i];
      return TEMP;
    }


    /**
    Sets this array equal to itself times a scalar value.
    */
    void operator*=(const T C)
    {
      for(int i=0; i<n_1*n_2; i++)
        pArray[i] *= C;
    }


    /**
    Sets this array equal to itself divided by a scalar value.
    */

    void operator/=(const T C)
    {
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
    Creates an array and sets it equal to another array.

    @param rhs array to set this array equal to.
    */
    Array2D( const Array2D<T> & rhs)
    {
      n_1 = 0;
      n_2 = 0;
      pArray = 0;
      allocate(rhs.n_1,rhs.n_2);
      operator=(rhs);
    }

    /**
    This function makes it easy to replace a row with an Array1D
    */
    void setRow(int whichRow, Array1D<T> & rhs)
    {
      if(whichRow >= n_1)  cerr << "Error: invalid row index\n";
      if(rhs.dim1() > n_2) cerr << "Error: rhs length is incorrect\n";
      memcpy(pArray + n_2*whichRow, rhs.array(), sizeof(T)*n_2);
    }

    /**
    This function takes advantage of the row-major format of the matricies
    to copy numRows worth of data from one Array2D to another
    */
    void setRows(int to, int from, int numRows, const Array2D<T> & rhs)
    {
      memcpy(pArray + n_2*to, rhs.pArray + n_2*from, sizeof(T)*n_2*numRows);
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
        return n_2*i + j;
      }

    /**
    Accesses element <code>(i,j)</code> of the array.
    */
    T& operator()(int i,int j)
    {
      assert(pArray != 0);
      return pArray[map(i,j)];
    }

    /**
    Accesses element <code>(i,j)</code> of the array, except this prevents
    modification of that element.
    */
    T get(int i, int j) const
        {
          assert(pArray != 0);
          return pArray[map(i,j)];
        }

    /**
    Prints the array to a stream.
    printf("%11.3e",pix[index +1]);
    printf("   ");
    */
    friend ostream& operator<<(ostream & strm, const Array2D<T> & rhs)
    {
      strm.precision(4);
      strm.setf(ios_base::scientific);
      int maxJ = 45;
      for(int i=0; i<rhs.n_1;i++)
        {
          for(int j=0; j<rhs.n_2 && j<maxJ; j++)
            {
              strm.width(20);
              strm.precision(12);
              strm << rhs.pArray[rhs.map(i,j)];
            }
          strm << endl;
        }
      return strm;
    }
  };

#endif
