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

#ifndef QMCReadAndEvaluateConfigs_H
#define QMCReadAndEvaluateConfigs_H

// Define the energy and log weight values to use when there is an attempt
// to evaluate a singular Jastrow
#define MAXIMUM_ENERGY_VALUE  1.0e20
#define MAXIMUM_LOG_WEIGHT_VALUE  50

#include <iostream>
#include <math.h>

#include "Array1D.h"
#include "QMCInput.h"
#include "QMCProperties.h"
#include "QMCJastrow.h"

using namespace std;

/**
  Calculates properties (QMCProperties) from walkers and related data saved
  to a file during a QMC calculation.
*/

class QMCReadAndEvaluateConfigs
{
public:
  /**
     Creates an instance of the class.
  */
  QMCReadAndEvaluateConfigs();


  /**
     Creates an instance of the class and initializes it.
     @param input data input to control the calculation.
  */
  QMCReadAndEvaluateConfigs(QMCInput *input, int cfgsToSkip);

  /**
     Initializes the object.
     @param input data input to control the calculation.
  */
  void initialize(QMCInput *input, int cfgsToSkip);

  /**
     Calculates properties (QMCProperties) for different parameter sets 
     from walkers and related data saved to a file during a QMC calculation.
     This function is called only by the root node.  The non-root nodes should
     call workerCalculateProperties().
     @param params array of parameters which parameterize the wavefunction.
     @param properties properties calculated from params and the 
     saved configurations.
  */
  void rootCalculateProperties(Array1D < Array1D<double> > &params, 
			       Array1D<QMCProperties> & properties);

  /**
     Calculates properties (QMCProperties) for different parameter sets 
     from walkers and related data saved to a file during a QMC calculation.
     This function is called only by the non-root nodes.  The root node should
     call rootCalculateProperties(params, properties).
  */
  void workerCalculateProperties();

private:
  // The input information
  QMCInput *Input;

  // The Jastrow to be calculated
  QMCJastrow Jastrow;

  // The number of configurations to be skipped.
  int configsToSkip;

  int Nelectrons;
  int Natoms;

  // Values read in from the configuration file
  Array2D<double> R;
  double D1;
  Array2D<double> D2;
  double lnJ;
  double PE;

  // given a set of parameters perform the necessary calcualtions and
  // add the results to the properties.
  void AddNewConfigToProperites(Array1D<double> &Params,
				QMCProperties &Properties);

  // Calculate the local energy of the current configuration with the currently
  // calculated jastrow
  double calc_E_Local_current();

  // Calculate the weight of the current configuration with the currently
  // calculated jastrow
  double calc_log_weight_current();

  // Calculate the properites from the configs for all the parameters in params
  // on the current node
  void locally_CalculateProperties(Array1D < Array1D<double> > &Params, 
	                                  Array1D<QMCProperties> & Properties);

  // perform an mpi reduce operation on the properties collected on each
  // processor
  void MPI_reduce( Array1D <QMCProperties> &local_Properties, 
		   Array1D < QMCProperties> &global_Properties);
};

#endif
