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

#ifndef QMCWALKERDATA_H
#define QMCWALKERDATA_H

#include <iostream>
#include <sstream>

#include "Array1D.h"
#include "Array2D.h"
#include "QMCDouble.h"
#include "QMCInput.h"

using namespace std;

/**
   The QMCWalkerData data type is meant to hold all the information
   calculated by QMCFunction that is useful to QMCWalker. This should
   in effect decouple QMCFunction and QMCWalker from each other enabling
   QMCFunction to be treated a little bit differently without significant
   modifications to QMCWalker.
*/
class QMCWalkerData {
 public:
  QMCWalkerData();
  ~QMCWalkerData();

  void initialize(int numDimensions, int numNucForceDim1, int numNucForceDim2);

  /**
     When we're proposing a trial move, we need to save some data
     from the original position. This function copies only the data
     we know we need to keep around in order to calculate the
     acceptance probability.

     However, if the move is not accepted, then we'll need all of this data
     back again if we're doing one electron updates.

     The question is whether it's faster to save all the data in this class,
     or to recalculate the data only if we need it. If we have a lot of
     electrons, then this represents a lot of data. Furthermore, the acceptance
     probability is probably high.
     
     I think it's faster to "update" the trial data back to the
     original data if the move is not accepted.
  */
  void partialCopy( const QMCWalkerData & rhs );

  /**
     When we're updating all electrons, whichE will be set to -1.
     Otherwise, it will be set to the index of the electron currently
     being updated.
  */
  int whichE;

  /**
     Several parts of the code use the interparticle
     distances and unit vectors, so we only want to calculate them once.

     rij is for electron-electron distances
     riA is for electron-nucleus distances
  */
  Array2D<double> rij;
  Array3D<double> rij_uvec;
  Array2D<double> riA;
  Array3D<double> riA_uvec;

  /**
     These are raw data for the calculation of the determinant
     part of the wavefunction. We need to save all of this if we're
     updating one electron at a time.
  */
  Array1D< Array2D<double> >   Dc_invA;
  Array1D< Array2D<double> >   Dc_invB;
  Array2D<qmcfloat>            D_xxA;
  Array2D<qmcfloat>            D_xxB;
  Array1D< Array2D<qmcfloat> > D_xA;
  Array1D< Array2D<qmcfloat> > D_xB;

  Array1D<double> DcA, DcB;
  Array1D<double> rDc_xxA, rDc_xxB;
  Array3D<double> rDc_xA,  rDc_xB;

  Array1D<double> Chi_DensityA, Chi_DensityB;
  
  /**
     These are the raw data specifically from the Jastrow
     part of the wavefunction.

     We want to remember all the calculated Jastrow data. This
     data is needed for an efficient implementation of one electron
     updates, since we need to remember the Jastrow evaluations for
     all the other electrons.

     Uij is for the electron-electron jastrows
     UiA is for the electron-nucleus jastrows
     UijA is for the 3 particle jastrows
  */

  /**
     \f$\ln(J)=\sum{u_{i,j}(r_{i,j})}\f$
  */
  Array2D<double> Uij;

  /**
     \f$\nabla\ln(J)=\nabla\sum{u_{i,j}(r_{i,j})}\f$
  */
  Array3D<double> Uij_x;

  /**
     \f$\nabla^2\ln(J)=\nabla^2\sum{u_{i,j}(r_{i,j})}\f$
  */
  Array2D<double> Uij_xx;

  Array2D<double> UiA;
  Array3D<double> UiA_x;
  Array2D<double> UiA_xx;

  Array2D<double> UijA;
  Array1D< Array2D<double> > UijA_x1;
  Array1D< Array2D<double> > UijA_x2;
  Array2D<double> UijA_xx;

  /**
     These are intermediate data specifically from
     the determinant part of the wavefunction.
  */
  QMCDouble       D;
  Array2D<double> D_x;
  double          D_xx;
  
  /**
     These are the intermediate data specifically from
     the Jastrow part of the wavefunction.
  */
  double          U;
  Array2D<double> U_x;
  double          U_xx;
  
  /**
     These are the final values for the wavefunction.
  */
  double localEnergy, kineticEnergy, potentialEnergy;
  double neEnergy, eeEnergy;
  double x2, y2, z2;

  QMCDouble psi;
  
  /**
     Gets the ratio of the wavefunction gradient to the wavefunction value at
     the last evaluated electronic configuration.  This is also known as the
     quantum force.
  */
  Array2D<double> gradPsiRatio;
  
  /**
     Gets a modified version of the ratio of the wavefunction gradient to the 
     wavefunction value at the last evaluated electronic configuration.  
     The modifications typically help deal with singularities near nodes,
     and the particular type of modification can be selected.  
     This is also known as the modified quantum force.
  */
  Array2D<double> modifiedGradPsiRatio;
  double modificationRatio;

  /**
     This is set to true if there is any known numerical
     problem with the wavefunction. If true, then it shouldn't
     be added to any averages.
  */
  bool singular;
  
  /**
     This holds derivatives of the Hamiltonian (Hellman-Feynmann theorem)
     evaluated with respect to each nucleus.
  */  
  Array2D<double> nuclearDerivatives;
  
  Array1D<double> chiDensity;

  //correlated sampling energies
  Array1D<double> cs_LocalEnergy;
  //correlated sampling weights
  Array1D<double> cs_Weights;

  /*
    We never need dPsi_dai directly, we only collect
    (1/Psi) * (dPsi_dai). Actually, it's easier
    to calculate the ratio anyway.

    So the data in this array is the derivative of Psi
    with respect to parameter ai, divided by Psi.
  */
  Array1D<double> rp_a;

  /*
    This will be the derivative of the kinetic
    energy with respect to parameter ai.
  */
  Array1D<double> p3_xxa;

  bool isSingular();

  double getModifiedLocalEnergy();
  void zero();
  void writeConfigs(Array2D<double> & Positions, double weight);
  friend ostream& operator<<(ostream & strm, const QMCWalkerData & rhs);
  
  /**
     Once we've moved an electron, we want to update the arrays
     in QMCWalkerData which contain the interparticle distances,
     as well as electron-electron unit vectors. The purpose is so that
     we only have to calculate this stuff once, and then everybody
     else just gets the saved data.
  */
  void updateDistances(Array2D<double> & R);

};

#endif





