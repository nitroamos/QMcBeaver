//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2000-2.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

#include "QMCCorrelatedSamplingVMCOptimization.h"

void QMCCorrelatedSamplingVMCOptimization::
optimize(QMCInput * input, int configsToSkip)
{
  //put initial Jastrow parameters in as the guess
  Array1D<double> Guess_Jastrow_parameters = input->JP.getParameters();

  if( input->flags.my_rank == 0 )
    {
      // initialize the objective function
      QMCObjectiveFunction ObjFunk;
      ObjFunk.initialize(input, configsToSkip);
      
      // get which optimization algorithm to use
      QMCOptimizationAlgorithm * optAlg = 
	QMCOptimizationFactory::optimizationAlgorithmFactory(ObjFunk, input);
      
      Guess_Jastrow_parameters = 
	optAlg->optimize(Guess_Jastrow_parameters);
      
      delete optAlg;
      optAlg = 0;
      
#ifdef PARALLEL  
      // Send a signal to the worker nodes to quit working on the objective
      // function
      // 1 signals the workers to execute RAEC.workerCalculateProperties
      // 0 signals the termination of the while loop
      
      int WorkSignal = 0;
      MPI_Bcast(&WorkSignal,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
    }
  else
    {
#ifdef PARALLEL
      // If you are a worker, wait for the root to give you work
      QMCReadAndEvaluateConfigs RAEC(input, configsToSkip);
      
      while(true)
	{
	  // Receive a signal from the root node
	  // 1 signals the workers to execute 
	  //      RAEC.workerCalculateProperties
	  // 0 signals the termination of the while loop
	  
	  int WorkSignal;
	  MPI_Bcast(&WorkSignal,1,MPI_INT,0,MPI_COMM_WORLD);
	  
	  if( WorkSignal == 1 ) RAEC.workerCalculateProperties();
	  else break;
	}
#endif
    }
  
#ifdef PARALLEL

  // Send the new Guess Jastrow Parameters to everyone
  MPI_Bcast(Guess_Jastrow_parameters.array(),Guess_Jastrow_parameters.dim1(),
	    MPI_DOUBLE,0,MPI_COMM_WORLD);

#endif

  input->JP.setParameterVector(Guess_Jastrow_parameters);
}




