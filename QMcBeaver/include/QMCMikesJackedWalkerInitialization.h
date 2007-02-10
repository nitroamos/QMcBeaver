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

#ifndef QMCMikesJackedWalkerInitialization_H
#define QMCMikesJackedWalkerInitialization_H

#include "QMCInput.h"
#include "QMCInitializeWalker.h"
#include "Random.h"


/**
  This is the algorithm made to initialize walkers.  It 
  is based on figuring out how many electrons should be on 
  each atom followed by putting them in a gaussian around 
  the atom.  This is by far a method which needs a serious
  overhaul.  This was a quick fix to initializing the 
  walkers and the ideas are borrowed from CASINO.  This 
  method of initializing is probably very inefficient.  
  This goes without mentioning how ugly the code is.  This
  is a great place for further future work.  A huge dent
  will likely be made on the "Initialization Catastrophe"
  problem here.
 */

class QMCMikesJackedWalkerInitialization : public QMCInitializeWalker
{
public:
  /**
    Create an instance of the clas and initializes it.

    @param input input data for the calculation
   */
  QMCMikesJackedWalkerInitialization( QMCInput *input );


  Array2D<double> initializeWalkerPosition();

private:
  QMCInput *Input;

  Array2D <double> electrons_and_radii();
  double covalent_radi(int ZZ);
};

#endif
