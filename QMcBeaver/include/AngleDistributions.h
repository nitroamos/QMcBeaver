#ifndef AngleDistributions_H
#define AngleDistributions_H

#include "Array1D.h"

/** 
  This class includes the distributions in phi and theta for 
  QMCDansWalkerInitialization.  The arrays are for the distributions s, spz1, 
  spz2, spy1, spy2, sp2zx1, sp2zx2, sp2zx3, sp2yx1, sp2yx2, sp2yx3, sp3a1, 
  sp3a2, sp3a3, sp3a4, sp3b1, sp3b2, sp3b3, sp3b4.
  These are the angle values that go with an x_array that goes from 0 to 1 by
  steps of .05.
*/

class AngleDistributions
{
 public:

  static Array1D<double> getPhiArray(int index);

  static Array1D<double> getThetaArray(int index);
  
};


#endif
