#ifndef RadialDistributions_H
#define RadialDistributions_H

#include "Array1D.h"

using namespace std;

/**
  This class contains the radial distribution functions for all atoms up to 
  Argon.  These are the r values that go with an x array that goes from 0 to 1
  in steps of .05.
*/

class RadialDistributions
{
 public:

  static Array1D<double> getRadialArray(int Z, int n);

};

#endif
