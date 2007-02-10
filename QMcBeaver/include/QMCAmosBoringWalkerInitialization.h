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

#ifndef QMCAmosBoringInitialization_H
#define QMCAmosBoringInitialization_H

#include "QMCInitializeWalker.h"
#include "QMCInput.h"
#include "Random.h"

using namespace std;

class QMCAmosBoringWalkerInitialization : public QMCInitializeWalker
{
public:
  /**
    Create an instance of the clas and initializes it.

    @param input input data for the calculation
   */
  QMCAmosBoringWalkerInitialization( QMCInput *input )
    {
      Input = input;
    }


  Array2D<double> initializeWalkerPosition();

private:
  QMCInput *Input;
};

#endif
