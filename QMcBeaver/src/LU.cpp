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

#include "LU.h"

// This should be a small number.  NR uses 1e-20 for a float which varies
// between 1.2e-38 and 3.4e38.
// A double varies between 2.2e-308 and 1.8e308 so a TINY value of 1e-250
// should easily be acceptable

#if defined SINGLEPRECISION || defined QMC_GPU
#define TINY 1e-35f
#else
#define TINY 1e-250
#endif

void ludcmp(Array2D <qmcfloat> &a, int *indx, double *d, bool *calcOK)
{
  int n=a.dim1();
  int i,j,k;
  int imax = -1;
  qmcfloat big,dum,temp;
  register qmcfloat sum;
  static qmcfloat *vv = new qmcfloat[n];
  qmcfloat one = 1.0;

  *calcOK = true;
  *d=1.0;

  //this section finds the largest value from each column
  //and puts its inverse in the vv vector.
  for (i=0;i<n;i++)
    {
      big=0.0;
      for (j=0;j<n;j++)
        if ((temp=fabs(a(i,j))) > big)
          big=temp;

      //checks if the column is full of zeros
      if (big == 0.0)
        {
          // cerr << "Singular matrix in routine ludcmp*********" << endl;
          *calcOK = false;
          return;
        }

      vv[i]=one/big;
    }

  //loop over columns
  for (j=0;j<n;j++)
    {

      //part 1: i < j
      for (i=0;i<j;i++)
        {
          sum=a(i,j);
          for (k=0;k<i;k++)
            sum -= a(i,k)*a(k,j);
          a(i,j)=sum;
        }

      //part 2: i >= j
      big=0.0;
      for (i=j;i<n;i++)
        {
          sum=a(i,j);
          for (k=0;k<j;k++)
            sum -= a(i,k)*a(k,j);
          a(i,j)=sum;

          //find the best row to pivot
          if ( (dum=vv[i]*fabs(sum) ) >= big)
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
          for (k=0;k<n;k++)
            {
              dum=a(imax,k);
              a(imax,k)=a(j,k);
              a(j,k)=dum;
            }
          //a.swapRows(imax,j);

          //indicate permutation change
          *d = -(*d);

          //update vv
          vv[imax]=vv[j];
        }
      indx[j]=imax;

      //make sure we don't divide by zero
      if (a(j,j) == 0.0)
        a(j,j)=TINY;

      //each element needs to be divided by it's diagonal
      if (j != n)
        {
          dum=one/(a(j,j));
          for (i=j+1;i<n;i++)
            a(i,j) *= dum;
        }

    }

}
#undef TINY

void lubksb(Array2D <qmcfloat>& a, int *indx, Array1D<qmcfloat> &b)
{
  int n=a.dim1();
  int ii=-1, ip;
  qmcfloat sum;

  for (int i=0;i<n;i++)
    {
      ip=indx[i];
      sum=b(ip);
      b(ip)=b(i);

      if (ii>=0)   // in previous version had if(ii<0) Checked with NR. ck
        for (int j=ii;j<=i-1;j++)
          sum -= a(i,j)*b(j);//we do forward-subsitution to find "y"
      else if (sum)//our first non-zero element
        ii=i;

      b(i)=sum;
    }

  for (int i=n-1;i>=0;i--)
    {
      sum=b(i);
      for (int j=i+1;j<n;j++)
        sum -= a(i,j)*b(j);//back-substitution to find "x"
      b(i)=sum/a(i,i);
    }
}

//I was using this to see if there was an in-code shortcut for calculating inverses. ~Amos
//I'll probably remove it at some point, cause I couldn't figure anything out.
void lubksbForInverse(Array2D <qmcfloat>& a, int *indx, Array1D<qmcfloat> &b, int which)
{
  //i can't find any ways to improve this if i'm specifically looking for an inverse.
  int n=a.dim1();
  //ii is meant to hold the index of the first non-zero element
  int ii=-1;

  int ip;
  qmcfloat sum;
  //when solving for the inverse using A.Ai = 1 where 1 is the identity, Ai is inverse
  //A = L.U
  //L.U.Ai = 1
  //L.Y = 1 so Y = Li, Li(i,j) =
  //U.Ai = Li, Ai = Ui.Li

  //we do forward-subsitution to find "y"
  for (int i=0;i<n;i++)
    {
      //first, unscramble the permutations
      ip=indx[i];
      sum=b(ip);
      b(ip)=b(i);

      if (ii>=0)   // in previous version had if(ii<0) Checked with NR. ck
        for (int j=ii;j<=i-1;j++)
          {
            sum -= a(i,j)*b(j);
          }
      else if (sum)//our first non-zero element
        ii=i;

      b(i)=sum;//diagonal is 1
    }

  //back-substitution to find "x"
  //once we've arrived at this stage, there is nothing
  //special about the fact that we're finding the inverse
  //is there a way to break the recursive nature of this?
  //i think the processor is having a hard time predicting branches
  for (int i=n-1;i>=0;i--)
    {
      sum=b(i);
      for (int j=i+1;j<n;j++)
        {
          sum -= a(i,j)*b(j);
        }
      b(i)=sum/a(i,i);
    }
}

double determinant(Array2D <qmcfloat> A, bool *calcOK)
{
  int n = A.dim1();
  if(n == 0) return 1;
  if(n == 1) return A(0,0);
  int *INDX = new int[n];
  double d;

  ludcmp(A,INDX,&d,calcOK);

  delete [] INDX;
  if(d == 0) return d;

  for(int i=0; i<n; i++)
    d *= A(i,i);

  return d;
}


Array2D <qmcfloat> inverse(Array2D <qmcfloat> a, bool *calcOK)
{
  int n=a.dim1();
  Array2D <qmcfloat> Y;
  Y.allocate(n,n);
  double d;
  int *INDX = new int[n];
  Array1D<qmcfloat> col(n);

  ludcmp(a,INDX,&d,calcOK);

  if( *calcOK )
    {
      for(int j=0;j<n;j++)
        {
          for(int i=0;i<n;i++) col(i) = 0.0;
          col(j) = 1.0;
          lubksb(a,INDX,col);
          for(int i=0;i<n;i++) Y(i,j) = col(i);
        }
    }

  delete [] INDX;
  return Y;
}


//LAPACK doesn't seem to be much help here...
#if defined USELAPACK
void determinant_and_inverse(Array2D<qmcfloat> &a, Array2D<qmcfloat> &inv,
                             double& det, bool *calcOK)
{
  static int dim = a.dim1();
  static Array1D<int> pivots = Array1D<int>(dim);
  inv = 0;
  for(int i=0; i<dim; i++)
    inv(i,i) = 1.0;

  /* int clapack_dgesv(const enum CBLAS_ORDER Order, const int N, const int NRHS,
                    double *A, const int lda, int *ipiv,
                    double *B, const int ldb);*/
#if defined SINGLEPRECISION || defined QMC_GPU
  int info = clapack_sgesv(CBLAS_ORDER(CblasRowMajor),dim,dim,
                           a.array(),dim,pivots.array(),inv.array(),dim);
#else
  int info = clapack_dgesv(CBLAS_ORDER(CblasRowMajor),dim,dim,
                           a.array(),dim,pivots.array(),inv.array(),dim);
#endif

  if(info == 0) *calcOK = true;
  else *calcOK = false;

  det = 1.0;
  for(int i=0; i<dim; i++)
    det *= a(i,i);
}
#else
/**
this is O(1*N^3)
WARNING!!!
this has been modified and returns the transpose of the inverse, not the inverse
this makes calculate_DerivativeRatios in QMCSlater.cpp look nicer
*/
void determinant_and_inverse(Array2D<qmcfloat> &a, Array2D<qmcfloat> &inv,
                             double& det, bool *calcOK)
{
  int n=a.dim1();
  double d;
  static int *INDX = new int[n];
  static Array1D<qmcfloat> col(n);

  ludcmp(a,INDX,&d,calcOK);

  inv.allocate(n,n);

  if( *calcOK )
    {
      for(int j=0;j<n;j++)
        {
          for(int i=0;i<n;i++) col(i) = 0.0;
          col(j) = 1.0;
          lubksb(a,INDX,col);
          //lubksbForInverse(a,INDX,col,j);
          //for(int i=0;i<n;i++) inv(i,j) = col(i);
          inv.setRow(j,col);
          //for(int i=0;i<n;i++) inv(j,i) = col(i);
        }
    }

  for(int i=0; i<n; i++) d*= a(i,i);
  det = d;

  //delete [] INDX;
}
#endif

void linearsolver(Array2D<qmcfloat> &a, Array1D<qmcfloat> &b, bool *calcOK)
{
  int n=a.dim1();
  double d;
  int *INDX = new int[n];

  ludcmp(a,INDX,&d,calcOK);
  lubksb(a,INDX,b);

  delete [] INDX;
}

