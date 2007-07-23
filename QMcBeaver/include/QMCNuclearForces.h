/*
  Copyright (c) Amos G. Anderson 2006
  Distributed under GNU general public license (GPL)
  No guarantee or warantee regarding usability or stability is expressed or implied.
  nitroamos@gmail.com
*/

#ifndef QMCNuclearForces_H
#define QMCNuclearForces_H

#include <iostream>

#include "Array1D.h"
#include "Array2D.h"
#include "QMCInput.h"
#include "Random.h"
#include "CubicSpline.h"
#include "QMCWalkerData.h"
#include "fastfunctions.h"

using namespace std;

/**
   This class will measure the nuclear forces for the system. This is essentially
   a port from the fortran code graciously provided by Simone Chiesa as used in
   S. Chiesa, D.M. Ceperley, and S. Zhang, Phys. Rev. Lett. 94, 036404 (2005)
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
     This is the function that will be evaluated for each sampling. It
     will calculate the force for several walkers at the same time.
     
     @param X \f$3N\f$ dimensional configuration of electrons represented by 
     a \f$N \times 3\f$ matrix
  */
  void evaluate(Array1D< QMCWalkerData *> &walkerData, 
		Array1D<Array2D<double> * > &xData, int num);
  
  /**
     This is the function that is probably only going to be called
     during initialization.
     
     @param X \f$3N\f$ dimensional configuration of electrons represented by 
     a \f$N \times 3\f$ matrix
  */
  void evaluate(QMCWalkerData &walkerData, Array2D<double> &xData);

  double getTemperTerm(int nuc, int q, double r);

  /**
     As described in CCZ05, we will be using a polynomial as a part of the tempering
     term used to calculate the force in Eq. 6. This function is called as a part
     of the initialization and calculates the coefficients of the polynomial.

     It will also call waveMemorization if we will be using Eq. 9.
  */
  void calcCoefficients(int whichNucleus);
  
  /**
     This function is only called when we're using the Slater modification
     described in Eq. 9.

     We want to decompose the wavefunction into s, px, py, and pz contributions.
     First, we choose a cutoff radius beyond which we are not interested. Then
     inside this sphere, we record what the contributions are for a series of
     radial "knot" points distributed uniformly inside the cutoff.

     We then create a spline so that we can evaluate these contributions cheaply.

     numKnots is the number of radial points upon which
     the wavefunction will be memorized.
  */
  void waveMemorization(int whichNucleus, int numKnots, double radialCutoff);
  
  /**
     Since the force terms that depend upon 2 nuclei don't change in a calculation,
     we only want to calculate them once.
  */
  void calculateNuclearContributions();
  
  /**
     For a series of electronic coordinates, this function will
     calculate the density.

     @param X \f$3N\f$ dimensional configuration of electrons represented by 
     a \f$N \times 3\f$ matrix
     @param densities the array of densities calculated for each point
  */
  void getDensities(Array2D<double> & X, Array1D<qmcfloat> & densities);
  
  /**
     A debugging function...
  */
  void printCubeLengths(Array2D<double> & points);

  /**
     A debugging function...
  */  
  void printPoints(Array2D<double> & points);
  
  /**
     This will create a cube centered at the origin. The distance
     between any two adjacent points will be length.
  */
  void generateCube(Array2D<double> & cube, double length);

  /**
     This will take the cube centered at the origin 
     and give it a random rotation preserving the distance
     between the vertices as well as the center.

     The parameter scale will then grow or shrink the cube.
  */
  void randomlyRotate(Array2D<double> & points, double scale);
  
  /**
     This is an arbitrary choice, so it can be static. It
     controls the resolution of our binning.
  */
  static int getNumBins()
    {
      return 100;
    }
  
  /**
     Sets two QMCNuclearForces objects equal.
     
     @param rhs object to set this object equal to
  */
  void operator=( const QMCNuclearForces & rhs );
  
 private:
  QMCInput *Input;
  QMCBasisFunction *BF;
  QMCWavefunction  *WF;
  
  /**
     This is the m parameter 'm' used in CCZ05
  */
  double fittingWeightM;

  /**
     We want to know how many sample points per unit area
     on the surface of the shell to use for each knot point.
  */
  double numSamplesPerArea;

  /**
     This is the basis set size 'M' from CCZ05. We want to
     permit it to be different for each atom.
  */
  Array1D<int> numPolyBasisF;//nuc

  /**
     These are the coefficients c_k or a_k for the basis set.
  */
  Array2D< Array1D<double> > basisCoeffs;//(nuc,q), index

  /**
     These are the knot's radial values for each nucleus. The last value
     in this array (independent of method) will be the radial cutoff 'R'.
  */
  Array1D< Array1D<double> > radialPoints;

  /**
     The memorized wavefunction as decomposed into
     s, px, py, and pz contributions.
  */
  Array2D< Array1D<double> > waveValuesHF;// (nuc,kind), radial index

  /**
     This is the interpolating waveValuesHF
  */
  Array2D< CubicSpline > waveValuesHFSpline;// (nuc,kind), radial index

  /**
     The precalculated nucleus contributions to the force.
  */
  Array2D<double> nucleusContributions;

  /**
     Temp variable used by getDensities
  */
  Array2D<qmcfloat> orbitals;

  /**
     Temp variable used by getDensities
  */
  Array2D<qmcfloat> Chi;

  /**
     Eventually, we probably want to remove this guy
  */
  CubicSpline spliner;
};

#endif
