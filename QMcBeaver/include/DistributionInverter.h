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

#ifndef DISTRIBUTION_INVERTER_H
#define DISTRIBUTION_INVERTER_H

#define TINY (1e-12)

#include <iostream>
#include <string>
#include <math.h>

#include "Array1D.h"
#include "Random.h"
#include "LinearSpline.h"

using namespace std;

/**
   Very crude 1-d distribution inverter.
*/

class DistributionInverter
{
public:  
  /**
     Creates a new uninitialized instance of this class.
  */
  DistributionInverter(){;} 

  /**
     Copy the rhs oject into this current object.

     @param rhs object to set this equal to.
  */
  void operator=( const DistributionInverter & rhs );

  /**
     Initialize the object with function values for the distribution 
     random numbers are to be generated with respect to.
  */

  void initialize(Array1D<double> &x_input, 
		  Array1D<double> &y_input);

  /**
     Generate a random number which is distributed with respect to the
     distribution this object was initialized with.

     @param iseed seed for generating the random number.
  */
  double random();
  
private:
  //F_inverse=inverse(F); give me a y and I'll give you an x
  LinearSpline F_inverse; 

  /**
     f is incoming distribution.
     F=integral(f); 
     F_inverse=inverse(F); 
     
     @param x_input is the incoming domain values
     @param y_input is the incoming function(x_input) values
  */
  void make_F_and_F_inverse(Array1D<double> &x_input, 
			    Array1D<double> &y_input);
};


#endif




