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

#ifndef QMCINITIALIZEWALKER_H
#define QMCINITIALIZEWALKER_H

#include "Array2D.h"


/** 
  Interface to algorithms which generate new walkers for a QMC calculation.
  A good algorithm will generate walkers which require little time for the
  Metropolis algorithm to be equilibrated.
  */

class QMCInitializeWalker
{
public:
  /**
    Virtual destructor.
    */
  virtual ~QMCInitializeWalker(){};

  /**
    Generates a new walker.

    @return new walker configuration represented by a \f$N \times 3\f$ matrix
    */
  
  virtual Array2D<double> initializeWalkerPosition() = 0;
};

#endif
