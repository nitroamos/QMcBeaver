#ifndef QMCDansWalkerInitialization_H
#define QMCDansWalkerInitialization_H

#include "QMCInput.h"
#include "QMCInitializeWalker.h"
#include "CubicSpline.h"
#include "AngleDistributions.h"
#include "RadialDistributions.h"

#include "Array1D.h"
#include "Array2D.h"

#include "random.h"

#include <fstream>
#include <string>

class QMCDansWalkerInitialization : public QMCInitializeWalker
{
 public:
  
  /** 
    Creates an instance of the class and initializes it.
    @param input input data for the calculation.
  */

  QMCDansWalkerInitialization(QMCInput * input);

  /**
    Generates a random initial configuration of the electrons.
    @return array of coordinates for the electrons.
  */

  Array2D<double> initializeWalkerPosition();

 private:

  QMCInput * Input;

  /**
    The array of cubic splines for the phi distributions made from the data in 
    the AngleDistributions class.
  */
  
  Array1D<CubicSpline> phiSplines;

  /** 
    This array indicates if a phi spline has been made for a distribution.  The
    element equals 0 if no spline has been made, and 1 if a spline has been
    made.
  */

  Array1D<int> phiSplinesMade;
 
  /**
    The array of cubic splines for the theta distributions made from the data 
    in the AngleDistributions class.
  */

  Array1D<CubicSpline> thetaSplines;

  /**
    This array indicates if a theta spline has been made for a distribution.
    The element equals 0 if no spline has been made, and 1 if a spline has been
    made.
  */

  Array1D<int> thetaSplinesMade;

  /**
    The array of cubic splines for the radial distributions made from the data 
    in the RadialDistributions class.  
    radialSplines[Atomic number, energy level]
  */

  Array2D<CubicSpline> radialSplines;

  /** 
    This array indicates if a radial spline has been made for a distribution.
    The element equals 0 if no spline has been made, and 1 if a spline has been
    made.
  */

  Array2D<int> radialSplinesMade;

  /** 
    The x_array goes from 0 to 1 in steps of .05 and is the same for all 
    distributions.
  */

  Array1D<double> x_array;

  void initializeArrays();

  /** 
    Distributes electrons in an energy level and gives them a random rotation.
    @param Z atomic charge of the nucleus.
    @param n energy level.
    @param nalpha number of alpha electrons.
    @param nbeta number of beta electrons.
    @return cartesian coordinates of the electrons in the energy level.
  */

  Array2D<double> dist_energy_level(int Z, int n, int nalpha, int nbeta);

  /**
    Assigns alpha and beta electrons to each nucleus based on the wavefunction.
    @return array with alpha and beta occupations for each nucleus.
  */

  Array2D<int> assign_electrons_to_nuclei();

  /**
    Distributes electrons around a nucleus.
    @param atomic_charge atomic charge of the nucleus.
    @param n_e total number of electrons.
    @param n_a number of alpha electrons.
    @param n_b number of beta electrons.
    @return cartesian coordinates of the electrons around this nucleus.
  */  

  Array2D<double> dist_center(int atomic_charge, int n_e, int n_a, int n_b);

  /** 
    These functions generate coordinates with respect to the distribution 
    indicated by the integer argument.  If no spline has been made yet for that
    distribution, it is made.
  */

  /** 
    Generates a phi coordinate with respect to the distribution indicated by 
    the index.  If no spline has been made yet for the distribution, it is 
    made.
    @param index of the distribution.
    @return phi coordinate.
  */

  double generatePhiCoordinate(int);
  
  /**
    Generates a theta coordinate with respect to the distribution indicated by
    the index.  If no spline has been made yet for the distribution, it is 
    made.
    @param index of the distribution.
    @return theta coordinate.
  */

  double generateThetaCoordinate(int);
  
  /**
    Generates radial coordinates for each electron with respect to the 
    distribution for the indicated atomic charge and energy level.
    @param Z atomic charge.
    @param n energy level.
    @param nelecs number of electrons.
    @return array of radial coordinates.
  */

  Array1D<double> generateRadialDistances(int,int,int);
};

#endif
