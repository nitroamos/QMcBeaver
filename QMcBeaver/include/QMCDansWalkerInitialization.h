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
  */

  QMCDansWalkerInitialization(QMCInput * input);

  /**
    Generates an initial configuration of the electrons.
  */

  Array2D<double> initializeWalkerPosition();

 private:

  QMCInput * Input;

  /**
    These are the arrays of cubic splines made from the distribution data in 
    the AngleDistributions and RadialDistributions classes.
    The SplinesMade arrays indicate if a spline has been made for a 
    distribution yet.  The element equals 0 if no spline has been made, and 1
    if a spline has been made.
  */
  
  Array1D<CubicSpline> phiSplines;
  Array1D<int> phiSplinesMade;

  Array1D<CubicSpline> thetaSplines;
  Array1D<int> thetaSplinesMade;

  Array2D<CubicSpline> radialSplines;
  Array2D<int> radialSplinesMade;

  /** 
    The x_array goes from 0 to 1 in steps of .05 and is the same for all 
    distributions.
  */

  Array1D<double> x_array;

  void initializeArrays();

  /** 
    Distributes electrons in an energy level and gives them a random rotation.
  */

  Array2D<double> dist_energy_level(int, int, int, int);

  /**
    From the wavefunction, decided how many alpha and beta electrons are on 
    each nucleus.
  */

  Array2D<int> assign_electrons_to_nuclei();

  /**
    Distributes electrons around a nucleus.
  */  

  Array2D<double> dist_center(int, int, int, int);

  /** 
    These functions generate coordinates with respect to the distribution 
    indicated by the integer argument.  If no spline has been made yet for that
    distribution, it is made.
  */

  double generatePhiCoordinate(int);
  double generateThetaCoordinate(int);
  Array1D<double> generateRadialDistances(int,int,int);

};

#endif
