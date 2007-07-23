//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2000-2.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

#include "QMCWolfeStepLengthSelector.h"

double QMCWolfeStepLengthSelector::stepLength(
    QMCObjectiveFunction *function, Array1D<double> & position,
    Array1D<double> & searchDirection, Array1D<double> & gradient, 
    Array2D<double> & unused,
    double functionValue)

{
  OF = function;

  double alpha_guess = 2.0;

  return wolfeStepLength(alpha_guess,position,searchDirection,gradient,
			 functionValue);
}

double QMCWolfeStepLengthSelector::calculateLineSearchObjectiveFunction(
                                   Array1D<double> & params,
                                   Array1D<double> & searchDirection,
                                   double stepLength)
{
  Array1D<double> p;

  p = params + (searchDirection*stepLength);

  return OF->evaluate(p).getScore();
} 

double QMCWolfeStepLengthSelector::calculateLineSearchObjectiveFunctionDerivative(
                                  Array1D<double> & params,
                                  Array1D<double> & searchDirection,
                                  double stepLength)
{
  const double epsilon = 1.0e-4;

  double phi = calculateLineSearchObjectiveFunction(params,searchDirection,
						    stepLength);
  double phi_epsilon = calculateLineSearchObjectiveFunction(params,
							   searchDirection,
							   stepLength+epsilon);

  return (phi_epsilon-phi)/epsilon;
}
 
double QMCWolfeStepLengthSelector::cubicInterpolateStep(double a_lo, 
							double a_hi, 
							double phi_0, 
							double phi_a_lo, 
							double phi_a_hi, 
							double phi_prime_0)
{
  // Here is a bunch of nonsense to find the cubic interpolated a_new
  // See the book for details
  
  // make sure the problem isn't singular!
  if( fabs(a_lo) < 1e-8 ) a_lo = 1e-8;
  if( fabs(a_hi) < 1e-8 ) a_hi = 1e-8;
  if( fabs(a_lo - a_hi) < 1e-8) return (a_lo+a_hi)/2.0; 
  
  Array2D<double> M(2,2);
  
  M(0,0) = a_lo * a_lo;
  M(0,1) = -a_hi * a_hi;
  M(1,0) = -a_lo * a_lo * a_lo;
  M(1,1) = a_hi * a_hi * a_hi;
  
  Array1D<double> v(2);
    
  v(0) = phi_a_hi - phi_0 - a_hi * phi_prime_0;
  v(1) = phi_a_lo - phi_0 - a_lo * phi_prime_0;
  
  double factor = a_lo * a_lo * a_hi * a_hi * ( a_hi - a_lo );
  
  Array1D<double> ab(2);

  for(int i=0; i<2; i++)
    {
      ab(i) = 0.0;

      for(int j=0; j<2; j++)
        {
          ab(i) = M(i,j) * v(j);
        }
    }

  ab /= factor;
  
  double desc  = ab(1)*ab(1)-3.0*ab(0)*phi_prime_0;
  
  if( desc < 0 ) return (a_hi+a_lo)/2.0;

  
  double a_new = (-ab(1)+sqrt(desc))/3.0/ab(0);
  
  // if the interpolated point isn't in the range, return the midpoint
  if( a_new < a_lo ) a_new = (a_lo + a_hi)/2.0;
  if( a_new > a_hi ) a_new = (a_lo + a_hi)/2.0;
  
  // if the new point is too close to either end scoot it over
  double interptol = 1e-6;
  
  if( fabs(a_new-a_lo) < interptol * fabs(a_hi-a_lo) )
    {
      a_new = a_lo + interptol;
    }
  
  if( fabs(a_new-a_hi) < interptol * fabs(a_hi-a_lo) )
    {
      a_new = a_hi - interptol;
    }
  
  return a_new;
}


double QMCWolfeStepLengthSelector::zoom(double a_lo, double a_hi, 
					double phi_0, 
					double phi_a_lo, double phi_a_hi, 
					double phi_prime_0,
					Array1D<double> & params, 
					Array1D<double> & searchDirection)
{
  const double zoomTol = 1e-7;   // numerical stability tolerance
  const int maxZoomSteps = 1000; // max steps in converging zoom
  const double c1 = 1.0e-4;      // sufficient decrease Wolfe parameter
  const double c2 = 0.1;         // curvature Wolfe parameter

  double a;

  for(int i=0; i<maxZoomSteps; i++)
    {
      // Check for neumerically unstable problem
      if( fabs(a_lo-a_hi) < zoomTol ) return (a_lo+a_hi)/2.0;
      
      // find a point in [a_lo,a_hi]
      a = cubicInterpolateStep( a_lo, a_hi, phi_0, phi_a_lo, 
                                phi_a_hi, phi_prime_0);
      
      // evaluate phi(a)
      double phi_a = 
        calculateLineSearchObjectiveFunction(params,searchDirection,a);  
      
      if( phi_a > phi_0 + c1 * a * phi_prime_0 || phi_a >= phi_a_lo )
        {
          a_hi     = a;
          phi_a_hi = phi_a;
        }
      else
        {
          double phi_prime_a = 
            calculateLineSearchObjectiveFunctionDerivative(params, 
							   searchDirection,a);
          
          if( fabs(phi_prime_a) <= -c2 * phi_prime_0 )
            {
              return a;
            }
          else if( phi_prime_a * (a_hi - a_lo) >= 0 )
            {
              a_hi     = a_lo;
              phi_a_hi = phi_a_lo;
            }

          a_lo     = a;
          phi_a_lo = phi_a;
        }
    }

  cerr << "WARNING: zoom(...) in QMCWolfeStepLengthSelector did not converge!" << endl;

  return a;  
}


double QMCWolfeStepLengthSelector::wolfeStepLength(double alpha_guess, 
					 Array1D<double> & params,
					 Array1D<double> & searchDirection,
					 Array1D<double> & gradient,
					 double functionValue)
{
  const double alpha_max = 2.0;    // maximum allowed step length
  const int MaxWolfeSteps = 1000;  // max steps to converge step length
  const double c1 = 1.0e-4;        // sufficient decrease Wolfe parameter
  const double c2 = 0.1;           // curvature Wolfe parameter



  // This is all right out of the book so look it up to make sense of it.

  double a_max = alpha_max;
  double phi_max = 
    calculateLineSearchObjectiveFunction(params,searchDirection,alpha_max);  

  double a_0 = 0.0;
  double a_1 = alpha_guess;
  
  double phi_0       = functionValue;
  double phi_prime_0 = gradient*searchDirection;

  double phi_a_0 = phi_0;
  double phi_a_1;
  double phi_prime_a_1;
  
  for( int i=0; i<MaxWolfeSteps; i++ )
    {
      phi_a_1 =         
        calculateLineSearchObjectiveFunction(params,searchDirection,a_1);  
      
      // Check for sufficient descent
      if( phi_a_1 > phi_0 + c1 * a_1 * phi_prime_0 ||
          ( phi_a_1 >= phi_a_0 && i>0 ) )
        {
          return zoom(a_0,a_1,phi_0,phi_a_0,phi_a_1,phi_prime_0,params,
                      searchDirection);
        }
      
      // Calculate the directional derivative at the trial point
      phi_prime_a_1 = 
            calculateLineSearchObjectiveFunctionDerivative(params, 
                               searchDirection,a_1);
      
      // Check the curvature condition
      if( fabs(phi_prime_a_1) <= -c2 * phi_prime_0 )
        {
          return a_1;
        }

      if( phi_prime_a_1 >= 0 )
        {
          return zoom(a_1,a_0,phi_0,phi_a_1,phi_a_0,phi_prime_0,params,
                      searchDirection);
        }
      
      // move a(i) to a(i-1)
      a_0     = a_1;
      phi_a_0 = phi_a_1;
      
      // choose next a in (a_i, a_max)
      a_1 = cubicInterpolateStep(a_0,a_max,phi_0,phi_a_0,phi_max,phi_prime_0);
    }

  cerr << "WARNING: wolfeStepLength(...) in QMCWolfeStepLengthSelector did "
       << "not converge!" << endl;

  return a_1;  
}










