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

#include "DistributionInverter.h"

void DistributionInverter::operator=( const DistributionInverter & rhs )
{
  F_inverse=rhs.F_inverse;
}


void DistributionInverter::initialize(Array1D<double> &x_input,
				      Array1D<double> &y_input)
{
  //check to make sure no y_input values are negative
  for(int i=0;i<x_input.size();i++)
    {
      if(y_input(i)<0.0)
	{
	  cout << "We can not use DistributionInverter class for a distribution"
	       << endl;
	  cout << "which is negative.  This makes no sense so we will quit."
	       << endl;
	  cout << "exit(1) in DistributionInverter::initializeWithFunctionValues()"
	       << endl;
	  cout << "i:\t" << i << "\tx_input():\t" << x_input(i) 
	       << "\ty_input():\t" << y_input(i) << endl;
	  exit(1);
	}
    }

  make_F_and_F_inverse(x_input,y_input);
}

void DistributionInverter::make_F_and_F_inverse(Array1D<double> &x_input, 
						Array1D<double> &y_input)
{
  int n=x_input.dim1();
  Array1D<double> X_INP;
  Array1D<double> Y_INP;
  X_INP.allocate(n);
  Y_INP.allocate(n);

  double total_integral=0.0;
  double L,R,A,B,int_AB;

  X_INP(0)=x_input(0);
  Y_INP(0)=total_integral;

  for(int i=0;i<(n-1);i++)
    {
      L=x_input(i);
      R=x_input(i+1);

      A = y_input(i);

      B = y_input(i+1);

      //trap rule for integration
      int_AB = (B+A)/2.0*(R-L);

      total_integral=total_integral+int_AB;
      
      X_INP(i+1)=R;
      Y_INP(i+1)=total_integral;
    }
      
  //normalize the global integral to 1.0
  for(int i=0;i<n;i++)
    {
      Y_INP(i)=Y_INP(i)/total_integral;
    }

  //make sure the integral is monotonically increasing
  for(int i=0;i<(n-1);i++)
    {
      if(TINY>(Y_INP(i+1)-Y_INP(i)))
	{
	  Y_INP(i+1)=Y_INP(i)+TINY;
	}
    }
        
  //these background "TINY" increases may have changed the total integral
  //normalize the global integral to 1.0 one more time
  for(int i=0;i<n;i++)
    {
      Y_INP(i)=Y_INP(i)/Y_INP(n-1);
    }

  F_inverse.initializeWithFunctionValues(Y_INP,X_INP);
}

//given a uniform random number this will return the original distribution
double DistributionInverter::random(long &iseed)
{
  F_inverse.evaluate(ran1(&iseed));
  return F_inverse.getFunctionValue();
}








