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
#include "QMCGreensRatioComponent.h"
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

  void initialize(QMCInput * input, int numDimensions,
		  int numNucForceDim1, int numNucForceDim2);

  QMCInput * Input;

  int whichE;
  Array1D< Array2D<double> >             D_invA;
  Array1D< Array2D<double> >             D_invB;
  Array1D< Array2D<qmcfloat> >           Laplacian_DA;
  Array1D< Array2D<qmcfloat> >           Laplacian_DB;
  Array1D< Array1D< Array2D<qmcfloat> > > Grad_DA;
  Array1D< Array1D< Array2D<qmcfloat> > > Grad_DB;

  Array1D<double> PsiA, PsiB;
  Array1D<double> Laplacian_PsiRatioA, Laplacian_PsiRatioB;
  Array3D<double> Grad_PsiRatioA, Grad_PsiRatioB;

  Array1D<double> Chi_DensityA, Chi_DensityB;


  Array2D<double> rij;
  Array3D<double> rij_uvec;
  Array2D<double> riI;
  Array3D<double> riI_uvec;
  
  Array2D<double> Uij;
  Array2D<double> Uij_x;
  Array2D<double> Uij_xx;

  double localEnergy, kineticEnergy, potentialEnergy;
  double neEnergy, eeEnergy;

  double SCF_Laplacian_PsiRatio;
  double lnJ;
  Array2D<double> SCF_Grad_PsiRatio;

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
  
  QMCGreensRatioComponent psi;
  bool isSingular;
  
  double modificationRatio;
  
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
  
  /**
     This holds derivatives of the Hamiltonian (Hellman-Feynmann theorem)
     evaluated with respect to each nucleus.
  */  
  Array2D<double> nuclearDerivatives;
  
  Array1D<double> chiDensity;

  double getModifiedLocalEnergy();
  void zero();
  void writeConfigs(Array2D<double> & Positions, double weight);
  friend ostream& operator<<(ostream & strm, const QMCWalkerData & rhs);
};

#endif





