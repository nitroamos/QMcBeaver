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

void QMCCorrelatedSamplingVMCOptimization::optimize(QMCInput * input,
						    QMCProperties & lastRun,
						    int configsToSkip)
{
  //put initial Jastrow parameters in as the guess
  Array1D<double> Guess_parameters = globalInput.getParameters();

  double value;

  QMCDerivativeProperties dp(&lastRun,0,0);

  Array1D<double> gradient;
  Array2D<double> hessian;

  value = lastRun.energy.getSeriallyCorrelatedVariance();

  if( globalInput.flags.optimize_Psi_criteria == "analytical_energy_variance" )
    {
      gradient = dp.getParameterGradient();
    }

  if( globalInput.flags.optimize_Psi_method == "analytical_energy_variance" ||
      globalInput.flags.optimize_Psi_method == "automatic" )
    {
      hessian  = dp.getParameterHessian();
    }

  if( globalInput.flags.my_rank == 0 )
    {
      // initialize the objective function
      QMCObjectiveFunction ObjFunk;
      ObjFunk.initialize(input, configsToSkip);
      
      // get which optimization algorithm to use
      QMCOptimizationAlgorithm * optAlg = 
	QMCOptimizationFactory::optimizationAlgorithmFactory(ObjFunk, input);
      
      Guess_parameters = 
	optAlg->optimize(Guess_parameters,
			 value,
			 gradient,
			 hessian);
      
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
  MPI_Bcast(Guess_parameters.array(),Guess_parameters.dim1(),
	    MPI_DOUBLE,0,MPI_COMM_WORLD);

#endif

  globalInput.setParameterVector(Guess_parameters);

  if(globalInput.flags.my_rank == 0)
    {
      double penalty = globalInput.JP.calculate_penalty_function();
      if(penalty >= 1e10)
	{
	  clog << endl << endl << endl;
	  clog << "Error: the Jastow's new guess parameters have bad poles (penalty = " << penalty << ")!" << endl;
	  clog << "  Parameters are: " << globalInput.JP.getParameters();
	  clog << "   its poles are: " << globalInput.JP.getPoles();
	  clog << "  Either your guess for initial Jastrow parameters is bad or you need to change your "
	       << "optimization choices." << endl;
	  exit(0);
	} else {
	  clog << "Notice: the new parameters have an acceptable singularity penalty = " << penalty << endl;
	}
    }
}




