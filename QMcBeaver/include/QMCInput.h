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
    Write this object's state out to a stream. The same format is 
    used as in the QMC input file.
 */
 friend ostream& operator<<(ostream & strm, QMCInput & Input);
};

#endif
 
