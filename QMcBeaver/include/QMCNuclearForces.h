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

#ifndef QMCNuclearForces_H
#define QMCNuclearForces_H

#include <iostream>

#include "Array1D.h"
#include "Array2D.h"
#include "LU.h"
#include "QMCInput.h"
#include "random.h"
#include "CubicSpline.h"
#include "QMCWalkerData.h"
#include "fastfunctions.h"

using namespace std;

/**
  The potential energy of the system.
  */

class QMCNuclearForces
{
public:
  /**
    Creates an instance of the class.
    */
  QMCNuclearForces();
  ~QMCNuclearForces(); 
  /**
    Initialize the object.

    @param input data input to control the calculation
    */
  void initialize(QMCInput *input);

  /**


    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix
    */
  void evaluate(Array1D< QMCWalkerData *> &walkerData, 
		Array1D<Array2D<double> * > &xData, int num);

  /**


    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix
    */
  void evaluate(QMCWalkerData &walkerData, Array2D<double> &xData);
        
  void calcCoefficients(int whichNucleus);

  void waveMemorization(int whichNucleus);

  void calculateNuclearContributions();
  
  void getDensities(Array2D<double> & X, Array1D<double> & densities);

  void printCubeLengths(Array2D<double> & points);

  void printPoints(Array2D<double> & points);

  void generateCube(Array2D<double> & cube, double length);

  void randomlyRotate(Array2D<double> & points, double scale);
  
  static int getNumBins()
  {
    return 2;
  }
  
  /**
    Sets two QMCPotential_Energy objects equal.

    @param rhs object to set this object equal to
  */
  void operator=( const QMCNuclearForces & rhs );

 private:
  QMCInput *Input;
  QMCBasisFunction *BF;
  QMCWavefunction  *WF;
  
  double fittingWeightM;
  int numKnots;
  double numSamplesPerArea;
  Array1D<int> numPolyBasisF;//nuc
  Array2D< Array1D<double> > basisCoeffs;//(nuc,q), index
  Array1D<double> rMaxNuc;//nuc
  Array1D<double> radialPoints;//radial index
  Array2D< Array1D<double> > waveValuesHF;// (nuc,kind), radial index
  Array2D< CubicSpline > waveValuesHFSpline;// (nuc,kind), radial index
  Array2D<double> nucleusContributions;
  Array2D<double> alphaOrbitals;
  Array2D<double> betaOrbitals;
  CubicSpline spliner;
  
  /* These are for collecting force density information*/
  double maxBin;
};

#endif
