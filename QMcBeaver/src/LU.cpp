
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

#define TINY 1e-250

void ludcmp(Array2D <double> &a, int *indx, double *d, bool *calcOK)
{
	int n=a.dim1();
	int i,j,k;
	int imax = -1;
	double big,dum,sum,temp;
	double *vv = new double[n];

	*calcOK = true;

	*d=1.0;
	for (i=0;i<n;i++) 
	{
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a(i,j))) > big) big=temp;
		if (big == 0.0) 
		{
			// cerr << "Singular matrix in routine ludcmp*********" << endl;
			*calcOK = false;
			return;
		}
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) 
	{
		for (i=0;i<j;i++) 
		{
			sum=a(i,j);
			for (k=0;k<i;k++) sum -= a(i,k)*a(k,j);
			a(i,j)=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) 
		{
			sum=a(i,j);
			for (k=0;k<j;k++)
				sum -= a(i,k)*a(k,j);
			a(i,j)=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) 
			{
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) 
			{
				dum=a(imax,k);
				a(imax,k)=a(j,k);
				a(j,k)=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a(j,j) == 0.0) a(j,j)=TINY;
		if (j != n) 
		{
			dum=1.0/(a(j,j));
			for (i=j+1;i<n;i++) a(i,j) *= dum;
		}
	}
	delete [] vv;
}
#undef TINY

void lubksb(Array2D <double>& a, int *indx, Array1D<double> &b)
{
	int n=a.dim1();
	int ii=-1,ip;
	double sum;

	for (int i=0;i<n;i++) 
	{
		ip=indx[i];
		sum=b(ip);
		b(ip)=b(i);
		if (ii>=0)   // in previous version had if(ii<0) Checked with NR. ck
			for (int j=ii;j<=i-1;j++) sum -= a(i,j)*b(j);
		else if (sum) ii=i;
		b(i)=sum;
	}
	for (int i=n-1;i>=0;i--) 
	{
		sum=b(i);
		for (int j=i+1;j<n;j++) sum -= a(i,j)*b(j);
		b(i)=sum/a(i,i);
	}
}


double determinant(Array2D <double> A, bool *calcOK)
{
	int n = A.dim1();
	if(n == 0) return 1;
	if(n == 1) return A(0,0);
	int *INDX = new int[n];
	double d;

	ludcmp(A,INDX,&d,calcOK);

	delete [] INDX;
	if(d == 0) return d;

	for(int i=0; i<n; i++) d *= A(i,i);

	return d;
}


Array2D <double> inverse(Array2D <double> a, bool *calcOK)
{
	int n=a.dim1();
	Array2D <double> Y;
	Y.allocate(n,n);
	double d;
	int *INDX = new int[n];
	Array1D<double> col(n);

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

/**
this is O(1*N^3)
WARNING!!!
this has been modified and returns the transpose of the inverse, not the inverse
this makes calculate_DerivativeRatios in QMCSlater.cpp look nicer
*/
void determinant_and_inverse(Array2D<double> a, Array2D<double> &inv, 
							 double& det, bool *calcOK)
{
	int n=a.dim1();
	double d;
	int *INDX = new int[n];
	Array1D<double> col(n);

	ludcmp(a,INDX,&d,calcOK);

	inv.allocate(n,n);

	if( *calcOK )
	{
		for(int j=0;j<n;j++)
		{
			for(int i=0;i<n;i++) col(i) = 0.0;
			col(j) = 1.0;
			lubksb(a,INDX,col);
			//for(int i=0;i<n;i++) inv(i,j) = col(i);
			for(int i=0;i<n;i++) inv(j,i) = col(i);
		}
	}

	for(int i=0; i<n; i++) d*= a(i,i);
	det = d;

	delete [] INDX;
}


void linearsolver(Array2D<double> &a, Array1D<double> &b, bool *calcOK)
{
	int n=a.dim1();
	double d;
	int *INDX = new int[n];

	ludcmp(a,INDX,&d,calcOK);
	lubksb(a,INDX,b);

	delete [] INDX;
}

