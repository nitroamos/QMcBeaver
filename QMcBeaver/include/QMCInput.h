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

#ifndef QMCInput_H
#define QMCInput_H

#include <string>
#include <fstream>

#include "QMCFlags.h"
#include "QMCMolecule.h"
#include "QMCBasisFunction.h"
#include "QMCWavefunction.h"
#include "QMCConfigIO.h"
#include "QMCJastrowParameters.h"

using namespace std;

/**
  All input to a QMC calculation.  Input includes molecular geometry,
  trial wavefunction, job control information, and other relevant 
  information.
*/

class QMCInput
{
public:
  /**
     Creates an instance of the class.
  */
  QMCInput();

  /**
     Information on how to perform the QMC calculation.
  */
  QMCFlags flags;

  /**
     Geometry and charges of the system being studied.
  */
  QMCMolecule Molecule;

  /**
     Basis set used for the trial wavefunction.
  */
  QMCBasisFunction BF;

  /**
     Trial wavefunction for the calculation.
  */
  QMCWavefunction WF;

  /**
     Jastrow function used in the QMC wavefunction.
  */
  QMCJastrowParameters JP;
  
  /**
     This object is used to output the configuation file.
  */
  QMCConfigIO outputer;

  /**
     The complete set of correlated sampling parameters.
     The set at 0 are for the guiding wavefunction.
  */
  Array1D< Array1D<double> > cs_Parameters;

  /**
     Save parameters related to a parallel MPI calculation 
     in this class so that they are available throughout the
     calculation.

     @param my_rank MPI rank of this processor
     @param nprocs  number of processors used in this calculation
  */
 void setMPIParameters(int my_rank, int nprocs);

 /**
    Load this object's state from a QMC input file and initialize the 
    object.

    @param inputfile QMC input file to load this object's state from.
 */
 void read(string inputfile);

 /**
    This function will make sure that the config is open
    for writing.
 */
 void openConfigFile();

 /**
    The total number of parameters that we'll be optimizing
    across all parts of the wavefunction.
 */
 int getNumberAIParameters();

  /**
    Gets the parameters describing the particle-particle interactions,
    and adds to them the CI coefficients.

    @return parameters for optimizing
  */
  Array1D<double> getAIParameters();

  /**
    Sets the parameters describing the particle-particle interactions,
    as well as the CI coefficients.

    @param params new set of parameters
  */
  void setAIParameters(Array1D<double> & params);

  /**
     Print an array, designating the order that the ai
     parameters are organized by.
  */
  void printAIParameters(ostream & strm,
			 string name,
			 int margin,
			 Array1D<double> & array,
			 bool forcePrintOrbitals);

  /**
     This function makes it easy to print segments of an array,
     where the segments correspond to different kinds of parameters.

     This function is used to print out one segment, with several
     customization options.
  */
  void printArray(ostream & strm, string name, int num, 
		  Array1D<double> & array, int & start, 
		  int margin, int width, int prec, int numPerRow);

  /**
     Prints out the number of parameters we're optimizing in the
     current step, and shows a break down of the numbers.
  */
  void printAISummary();

  /**
     Write this object's state out to a stream. The same format is 
     used as in the QMC input file.
  */
  friend ostream& operator<<(ostream & strm, QMCInput & Input);
};

extern QMCInput globalInput;

#endif
 
