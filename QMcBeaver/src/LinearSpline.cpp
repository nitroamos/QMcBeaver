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

#include "LinearSpline.h"

void LinearSpline::operator=( const LinearSpline & rhs )
{
  n=rhs.n;
  f=rhs.f;

  x_list=rhs.x_list;
  y_list=rhs.y_list;
  a0_list=rhs.a0_list;
  a1_list=rhs.a1_list;
}


void LinearSpline::initializeWithFunctionValues(Array1D<double> &x_input, 
                Array1D<double> &y_input)
{
  //initialize the spline tables of values

  n=x_input.dim1();

  x_list.allocate(n);
  y_list.allocate(n);
  a0_list.allocate(n-1);
  a1_list.allocate(n-1);

  for(int i=0;i<n;i++)
    {
      x_list(i)=x_input(i);
      y_list(i)=y_input(i);
    }

  for(int i=0;i<n-1;i++)
    {
      a0_list(i)=0.0;
      a1_list(i)=0.0;
    }

  //if we look at the intervals xj-1,xj,xj+1
  //xj-1 <---h---> xj <---d---> xj+1 are the spacings

  //The spline function, s(x) has this form in each interval:
  //sj(x)=a0+a1(x-xj)
  //with a different set of a0 and a1 coeff for each interval

  //To get the coeff lists values we need to evaluate the following
  //a0=fj
  //a1=f'j

  double d,fj,fj1;

  //coeff(i) is good for the interval i to i+1
  for(int i=0;i<n-1;i++)
    {
      d    = x_list(i+1)-x_list(i);
      fj   = y_list(i);
      fj1  = y_list(i+1);
      a0_list(i)=fj;
      a1_list(i)=(fj1-fj)/d;
    }
}


void LinearSpline::evaluate(double x)
{
  //find which index in the table by means of bisection
  int klo,khi,k;

  //push any value outside the spline range into the last value of the spline

  if(x>=x_list(n-1))
    {
      x=x_list(n-1);
    }
  else if(x<=x_list(0))
    {
      x=x_list(0);
    }

  klo=0;
  khi=n-1;

  while (khi-klo > 1) 
    {
      k=(khi+klo) >> 1;
      if (x_list(k) > x) khi=k;
      else klo=k;
    }

  if(khi>n || khi<1 || klo>(n-1) || klo<0)
    {
      cerr<<"bad values for klo and khi"<<endl;
      cerr<<"klo= "<<klo<<"\tkhi= "<<khi<<endl;
      cerr<<"x= "<<x<<endl;
      exit(1);
    }

  double r,a0,a1;

  r=x-x_list(klo);
  a0=a0_list(klo);
  a1=a1_list(klo);

  f = a0+a1*r;
  df = a1;
}

double LinearSpline::getFunctionValue()
{
  return f;
}

double LinearSpline::getFirstDerivativeValue()
{
  return df;
}

double LinearSpline::getSecondDerivativeValue()
{
  return 0.0;
}


