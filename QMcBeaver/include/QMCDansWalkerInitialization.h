#ifndef QMCDansWalkerInitialization_H
#define QMCDansWalkerInitialization_H

#include "QMCInput.h"
#include "QMCInitializeWalker.h"

#include "Array1D.h"
#include "Array2D.h"

#include "random.h"

#include <fstream>
#include <string>

// Note: Using dans_walker_initialization requires that radial_dist_arrays and
// angle_dist_arrays be copied from the /QMcBeaver/src/ directory to the 
// working directory.  Radial arrays have been computed for atoms of atomic 
// number up to 18.

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
    An overloaded function to distribute one or more points with respect to a 
    distribution expressed as a 2D array or an element of a 2D array.
  */
  
  Array1D<double> dist_wrt_array(int, Array2D<double>);
  double dist_wrt_array(Array2D<double>, int);

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

};

#endif
