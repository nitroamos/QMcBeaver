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
						    QMCFutureWalkingProperties & fwLastRun,
						    int configsToSkip)
{
  //put initial Jastrow parameters in as the guess
  Array1D<double> orig_parameters = globalInput.getAIParameters();
  Array1D<double> Guess_parameters;

  double value;

  QMCDerivativeProperties dp(&lastRun,&fwLastRun,0);

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

      bool acceptable = false;
      int iter = 0;      
      double a_diag_factor = 1.0;
      do {
	Guess_parameters = 
	  optAlg->optimize(orig_parameters,
			   value,
			   gradient,
			   hessian,
			   a_diag_factor);

	globalInput.setAIParameters(Guess_parameters);	
	double penalty = globalInput.JP.calculate_penalty_function();
	if(penalty >= 1e10)
	  {
	    clog << endl << endl << endl;
	    clog << "Error: the Jastow's new guess parameters have bad poles (penalty = ";
	    clog.width(20);
	    clog.precision(10);
	    clog << penalty << ")!" << endl;
	    clog << "  Parameters are: " << globalInput.JP.getJWParameters();
	    clog << "   its poles are: " << globalInput.JP.getPoles();

	    /*
	      This will help us take more delicate steps, although
	      the optimziation will be slower
	    */
	    a_diag_factor *= 4.0;
	    acceptable = false;
	  } else {
	    clog << "Notice: the new parameters have an acceptable singularity penalty = " << penalty << endl;
	    acceptable = true;
	  }
	iter++;
      } while( !acceptable && iter < 10);
      
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

  globalInput.setAIParameters(Guess_parameters);

  double penalty = globalInput.JP.calculate_penalty_function();
  if(penalty >= 1e10)
    {
      clog << "Error: unable to produce acceptable parameters, quitting..." << endl;
      clog << "  Either your guess for initial Jastrow parameters is bad or you need to change your "
	   << "optimization choices." << endl;
      exit(0);
    }
}




