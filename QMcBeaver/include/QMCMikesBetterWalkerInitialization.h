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

#ifndef QMCMikesBetterWalkerInitialization_H
#define QMCMikesBetterWalkerInitialization_H

#include "QMCInput.h"
#include "QMCInitializeWalker.h"
#include "Random.h"
#include "AtomicOrbitalInverter.h"

using namespace std;

/**
  This is the algorithm made to initialize walkers.  It 
  is based on figuring out how many electrons should be on 
  each atom followed by putting them in atomic orbitals on
  each atom.  This puts the right number of electrons in 
  each type of atomic orbital on each atom.
 */

class QMCMikesBetterWalkerInitialization : public QMCInitializeWalker
{
public:
  /**
    Create an instance of the clas and initializes it.

    @param input input data for the calculation
   */
  QMCMikesBetterWalkerInitialization( QMCInput *input );


  Array2D<double> initializeWalkerPosition();

private:
  QMCInput *Input;
  AtomicOrbitalInverter AOI;

  Array3D<double> initializeBunchOfWalkersPosition();
  Array2D<double> FindBestWalker(Array3D<double> &bunchR);
  void FixConstraints(Array2D<double> &Occupations);
  double ObjectiveFunctionForWalkers(Array3D<double> &bunchR,Array2D<double> &Occupations);
  void GradObjectiveFunctionForWalkers(Array3D<double> &bunchR,
				  Array2D<double> &Occupations,
				  Array2D<double> &GradOccupations);
  void BoundGradOccupations(Array2D<double> &GradOccupations,double bound);
  void MoveOccupations(Array2D<double> &Occupations,
		       Array2D<double> &GradOccupations,
		       double dr);
  double Energy_parallel(double r);
  double Energy_opposite(double r);
  double Energy_el_nuclr(double r, int charge);
};

#endif
