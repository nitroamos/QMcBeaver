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

#include "Array2D.h"

static const bool USE_KAHAN = false;

template <class T>
T* Array2D<T>::array()
{
  return pArray;
}

#ifdef USEATLAS
/**This matrix multiplication requires rhs to be transposed.*/
//template <class T>
Array2D<double> Array2D<double>::operator*(const Array2D<double> & rhs)
{
  if(n_2 != rhs.n_2)
    {
      cerr << "ERROR: Matrix multiplication: " << n_1 << "x"
      << n_2 << " * " << rhs.n_2 << "x" << rhs.n_1 << endl;
      exit(1);
    }
  Array2D<double> TEMP(n_1,rhs.n_1);
  TEMP = 0;
  cblas_dgemm(CBLAS_ORDER(CblasRowMajor),
              CBLAS_TRANSPOSE(CblasNoTrans),CBLAS_TRANSPOSE(CblasTrans),
              n_1, rhs.n_1, n_2,
              1.0, pArray, n_2,
              rhs.pArray, rhs.n_2,
              0.0, TEMP.pArray, TEMP.n_2);
  return TEMP;
}

/**This matrix multiplication requires rhs to be transposed.*/
Array2D<float> Array2D<float>::operator*(const Array2D<float> & rhs)
{
  if(n_2 != rhs.n_2)
    {
      cerr << "ERROR: Matrix multiplication: " << n_1 << "x"
      << n_2 << " * " << rhs.n_2 << "x" << rhs.n_1 << endl;
      exit(1);
    }
  Array2D<float> TEMP(n_1,rhs.n_1);
  TEMP = 0;
  cblas_sgemm(CBLAS_ORDER(CblasRowMajor),
              CBLAS_TRANSPOSE(CblasNoTrans),CBLAS_TRANSPOSE(CblasTrans),
              n_1, rhs.n_1, n_2,
              1.0, pArray, n_2,
              rhs.pArray, rhs.n_2,
              0.0, TEMP.pArray, TEMP.n_2);
  return TEMP;
}

/**rhs is NOT transposed*/
Array2D<double> Array2D<double>::multiply(const Array2D<double> & rhs)
{
  if(n_2 != rhs.n_1)
    {
      cerr << "ERROR: Matrix multiplication: " << n_1 << "x"
      << n_2 << " * " << rhs.n_1 << "x" << rhs.n_2 << endl;
      exit(1);
    }
  Array2D<double> TEMP(n_1,rhs.n_2);
  TEMP = 0;
  cblas_dgemm(CBLAS_ORDER(CblasRowMajor),
              CBLAS_TRANSPOSE(CblasNoTrans),CBLAS_TRANSPOSE(CblasNoTrans),
              n_1, rhs.n_2, n_2,
              1.0, pArray, n_2,
              rhs.pArray, rhs.n_2,
              0.0, TEMP.pArray, TEMP.n_2);
  return TEMP;
}

/**rhs is NOT transposed*/
Array2D<float> Array2D<float>::multiply(const Array2D<float> & rhs)
{
  if(n_2 != rhs.n_1)
    {
      cerr << "ERROR: Matrix multiplication: " << n_1 << "x"
      << n_2 << " * " << rhs.n_1 << "x" << rhs.n_2 << endl;
      exit(1);
    }
  Array2D<float> TEMP(n_1,rhs.n_2);
  TEMP = 0;
  cblas_sgemm(CBLAS_ORDER(CblasRowMajor),
              CBLAS_TRANSPOSE(CblasNoTrans),CBLAS_TRANSPOSE(CblasNoTrans),
              n_1, rhs.n_2, n_2,
              1.0, pArray, n_2,
              rhs.pArray, rhs.n_2,
              0.0, TEMP.pArray, TEMP.n_2);
  return TEMP;
}

double Array2D<double>::dotAllElectrons(const Array2D<double> & rhs)
{
  return cblas_ddot(n_1*n_2, pArray, 1, rhs.pArray, 1);
}

double Array2D<double>::dotOneElectron(const Array2D<double> & rhs, int whichElectron)
{
  return cblas_ddot(n_1, pArray + whichElectron*n_2, 1, rhs.pArray + whichElectron*n_2, 1);
}

float Array2D<float>::dotAllElectrons(const Array2D<float> & rhs)
{
  return cblas_sdot(n_1*n_2, pArray, 1, rhs.pArray, 1);
}

float Array2D<float>::dotOneElectron(const Array2D<float> & rhs, int whichElectron)
{
  return cblas_sdot(n_1, pArray + whichElectron*n_2, 1, rhs.pArray + whichElectron*n_2, 1);
}

#else

/**This requires rhs to be transposed*/
template <class T>
Array2D<T> Array2D<T>::operator*(const Array2D<T> & rhs)
{
  if(n_2 != rhs.n_2)
    {
      cerr << "ERROR: Matrix multiplication: " << n_1 << "x"
      << n_2 << " * " << rhs.n_2 << "x" << rhs.n_1 << endl;
      exit(1);
    }
  int M = n_1;
  int N = rhs.n_1;
  int K = n_2;
  Array2D<T> TEMP(M,N);
  T * A = pArray;
  T * B = rhs.pArray;
  T * C = TEMP.pArray;
  int i, j, k;
  if(USE_KAHAN)
    {
      T Cee;
      T Why;
      T Tee;
      T temp;
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
                  //this prevents the compiler from optimizing this line away
                  if(Why    !=100) Tee = cij + Why;
                  if(B_j[k] != 99) temp = Tee - cij;
                  if(Ai_[k] != 98) Cee = temp - Why;
                  cij = Tee;
                }
              C[i*N + j] = cij;
            }
        }
    }
  else
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
  return TEMP;
}

/**rhs is NOT transposed*/
template <class T>
Array2D<T> Array2D<T>::multiply(const Array2D<T> & rhs)
{
  if(n_2 != rhs.n_1)
    {
      cerr << "ERROR: Matrix multiplication: " << n_1 << "x"
      << n_2 << " * " << rhs.n_1 << "x" << rhs.n_2 << endl;
      exit(1);
    }
  int M = n_1;
  int N = rhs.n_2;
  int K = n_2;
  Array2D<T> TEMP(M,N);
  T * A = pArray;
  T * B = rhs.pArray;
  T * C = TEMP.pArray;
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
  return TEMP;
}

template <class T>
T Array2D<T>::dotAllElectrons(const Array2D<T> & rhs)
{
  register T temp = 0;
  for(int i=0; i<n_1*n_2; i++)
    temp += pArray[i]*rhs.pArray[i];
  return temp;
}

template <class T>
T Array2D<T>::dotOneElectron(const Array2D<T> & rhs, int whichElectron)
{
  register T temp = 0;
  T * l = pArray + whichElectron*n_2;
  T * r = rhs.pArray + whichElectron*n_2;
  for(int i=0; i<n_1; i++)
    temp += l[i]*r[i];
  return temp;
}
#endif

template <class T>
Array2D<T> Array2D<T>::operator*(const T C)
{
  Array2D<T> TEMP(n_1,n_2);
  for(int i=0; i<n_1*n_2; i++)
    TEMP.pArray[i] = C*pArray[i];
  return TEMP;
}

template <class T>
void Array2D<T>::operator*=(const T C)
{
  for(int i=0; i<n_1*n_2; i++)
    pArray[i] *= C;
}

template <class T>
void Array2D<T>::operator/=(const T C)
{
  T inv = (T)(1.0)/C;
  operator*=(inv);
}

template <class T>
void Array2D<T>::setRow(int whichRow, Array1D<T> & rhs)
{
  if(whichRow >= n_1)
    {
      cerr << "Error: invalid row index\n";
    }

  if(rhs.dim1() > n_2)
    {
      //this is not necessarily a problem since for the copy, we only need n_2
      //elements...
      //if numAlpha != numBeta, then LU will produce an error in this function
      //cerr << "Error: rhs length is too large in " __FUNCTION__ << "\n";
    }

  if(rhs.dim1() < n_2)
    {
      cerr << "Error: rhs length is too small in " << __FUNCTION__ << "\n";
    }

  memcpy(pArray + n_2*whichRow, rhs.array(), sizeof(T)*n_2);
}

template <class T>
void Array2D<T>::setRows(int to, int from, int numRows, const Array2D<T> & rhs)
{
  memcpy(pArray + n_2*to, rhs.pArray + n_2*from, sizeof(T)*n_2*numRows);
}

template <class T>
void Array2D<T>::setColumn(int to, int from, const Array2D<T> & rhs)
{
  for(int i=0; i<n_1; i++)
    pArray[map(i,to)] = rhs.pArray[rhs.map(i,from)];
}

template <class T>
ostream& operator<<(std::ostream & strm, const Array2D<T> & rhs)
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

template <class T>
void Array2D<T>::ludcmp(int *indx, double *d, bool *calcOK)
{
  int i,j,k;
  int imax = -1;
  T big,dum,temp;
  register T sum;
  static T *vv = new T[n_1];
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

}
#undef REALLYTINY

template <class T>
void Array2D<T>::lubksb(int *indx, Array1D<T> &b)
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

//LAPACK doesn't seem to be much help here...
#if defined USELAPACK
void Array2D<double>::determinant_and_inverse(Array2D<double> &inv, double& det, bool *calcOK)
{
  static Array1D<int> pivots = Array1D<int>(n_1);
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
}

void Array2D<float>::determinant_and_inverse(Array2D<float> &inv, double& det, bool *calcOK)
{
  static Array1D<int> pivots = Array1D<int>(n_1);
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
}
#else
/**
this is O(1*N^3)
WARNING!!!
this has been modified and returns the transpose of the inverse, not the inverse
this makes calculate_DerivativeRatios in QMCSlater.cpp look nicer
*/
template <class T>
void Array2D<T>::determinant_and_inverse(Array2D<T> &inv, double& det, bool *calcOK)
{
  int n=dim1();
  double d;
  static int *INDX = new int[n];
  static Array1D<T> col(n);
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

  //delete [] INDX;
}
#endif

template class Array2D<double>;
template class Array2D<float>;
template class Array2D<int>;
