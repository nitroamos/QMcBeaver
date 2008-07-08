#ifndef QMCDansWalkerInitialization_H
#define QMCDansWalkerInitialization_H

#include "QMCInput.h"
#include "QMCInitializeWalker.h"
#include "CubicSpline.h"
#include "AngleDistributions.h"
#include "RadialDistributions.h"

#include "Array1D.h"
#include "Array2D.h"

#include "Random.h"

#include <fstream>
#include <string>

class QMCDansWalkerInitialization : public QMCInitializeWalker
{
 public:
  
  ~QMCDansWalkerInitialization();

  /** 
    Creates an instance of the class and initializes it.
    @param input input data for the calculation.
  */
  QMCDansWalkerInitialization(QMCInput * input);

  /**
    Generates an initial configuration of the electrons.
    @return 2D array with initial electron positions.
  */
  Array2D<double> initializeWalkerPosition();

 private:

  QMCInput * Input;

  /**
    The array of cubic splines made for the phi distributions made 
    from the data in the AngleDistributions class.
  */
  static Array1D<CubicSpline> phiSplines;

  /**
    The phiSplinesMade array indicates if a spline has been made for a 
    distribution.  The element equals 0 if no spline has been made, and 1
    if a spline has been made.
  */
  static Array1D<int> phiSplinesMade;

  /**
    The array of cubic splines made for the theta distributions made 
    from the data in the AngleDistributions class.
  */    
  static Array1D<CubicSpline> thetaSplines;

  /**
    The thetaSplinesMade array indicates if a spline has been made for a 
    distribution.  The element equals 0 if no spline has been made, and 1
    if a spline has been made.
  */
  static Array1D<int> thetaSplinesMade;

  /**
    The array of cubic splines made for the radial distributions made 
    from the data in the RadialDistributions class.
  */ 
  static Array2D<CubicSpline> radialSplines;

  /**
    The radialSplinesMade array indicates if a spline has been made for a 
    distribution.  The element equals 0 if no spline has been made, and 1
    if a spline has been made.
  */
  static Array2D<int> radialSplinesMade;

  /** 
    The x_array goes from 0 to 1 in steps of .05 and is the same for all 
    distributions.
  */
  static Array1D<double> x_array;

  static void initializeArrays();

  static bool arraysInitialized;

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
    From the wavefunction, decides how many alpha and beta electrons are on 
    each nucleus.
    @return 2D array with alpha and beta occupation of each nucleus.
  */
  Array2D<int> assign_electrons_to_nuclei();

  /**
    Distributes electrons around a nucleus.
    @param atomic_charge atomic charge of the nucleus
    @param eff_charge effective charge of the nucleus
    @n_e total number of electrons
    @n_a number of alpha electrons
    @n_b number of beta electrons

    @return cartesian coordinates of the electrons around this nucleus.
  */  
  Array2D<double> dist_center(int atomic_charge, int eff_charge, int n_e, int n_a, int n_b);

  /** 
    Generates a phi coordinate with respect to the distribution indicated by 
    the index.  If no spline has been made yet for that distribution, it is 
    made.
    @param index of the distribution.
    @return phi coordinate.
  */
  double generatePhiCoordinate(int index);

  /** 
    Generates a theta coordinate with respect to the distribution indicated by
    the index.  If no spline has been made for that distribution, it is made. 
    @param index of the distribution.
    @return theta coordinate.
  */
  double generateThetaCoordinate(int index);

  /** 
    Generates radial coordinates for each electron with respect to the 
    distribution for the indicated atomic charge and energy level.
    @param Z atomic charge.
    @param n energy level.
    @param nelecs number of electrons.
    @return array of radial coordinates.
  */
  Array1D<double> generateRadialDistances(int Z,int n,int nelecs);

};

#endif
