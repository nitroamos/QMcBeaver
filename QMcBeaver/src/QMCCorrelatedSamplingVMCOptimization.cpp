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

int QMCCorrelatedSamplingVMCOptimization::optStep = 1;

void QMCCorrelatedSamplingVMCOptimization::optimize(QMCInput * input,
						    QMCProperties & lastRun,
						    QMCFutureWalkingProperties & fwLastRun,
						    int configsToSkip)
{
  double rel = 0;
  double abs = 0;

  //put initial Jastrow parameters in as the guess
  Array1D<double> orig_parameters = globalInput.getAIParameters();
  Array1D<double> Guess_parameters(orig_parameters.dim1());

  QMCDerivativeProperties dp(&lastRun,&fwLastRun,0);
  
  if(optStep == 1)
    {
      if(fabs(globalInput.WF.getCINorm() - 1.0) > 1.0e-7)
	clog << "Notice: CI norm = " << globalInput.WF.getCINorm() << endl;
      clog.precision(12);
      clog.width(20);
      clog.unsetf(ios::fixed);
      clog.setf(ios::scientific);
      globalInput.printAIParameters(clog,"Best Step params:",20,orig_parameters,false);
      clog << endl << endl;
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
      double penalty, norm;
      do {
	Guess_parameters = 
	  optAlg->optimize(orig_parameters,dp,
			   a_diag_factor,
			   optStep);

	rel = 0;
	abs = 0;
	double num = 0;
	for(int i=0; i<Guess_parameters.dim1(); i++)
	  {
	    num += 1.0;
	    double g = Guess_parameters(i);
	    double o = orig_parameters(i);
	    double v = g - o;
	    abs += v*v;
	    if( fabs(o) > 1.0e-100 )
	      {
		v /= o;
		rel += v*v;	  
	      }
	  }
	abs = sqrt(abs)/num;
	rel = sqrt(rel)/num;

	clog.unsetf(ios::scientific);
	clog << "(Step = " << setw(3) << optStep << "):"
	     << " abs delta = " << setw(20) << abs
	     << " rel delta = " << setw(20) << rel << endl;
	
	globalInput.setAIParameters(Guess_parameters);	
	penalty = globalInput.JP.calculate_penalty_function();
	norm = globalInput.WF.getCINorm();
	acceptable = true;

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
	  }

	iter++;
      } while( !acceptable && iter < globalInput.flags.optimization_max_iterations);

      if(acceptable && fabs(globalInput.WF.getCINorm() - 1.0) > 1.0e-7)
	clog << "(Step = " << setw(3) << optStep << "): CI norm = " << globalInput.WF.getCINorm() << endl;

      if(acceptable && fabs(penalty) > 1.0e-10)
	clog << "Notice: the new parameters have an acceptable singularity penalty = " << penalty << endl;
      
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

  /*
    Some things might have been changed, so all the processors need to be
    apprised of the new data.
  */
  MPI_Bcast(&globalInput.flags.calculate_Derivatives,1,MPI_INT,0,MPI_COMM_WORLD);

  int numCS = globalInput.cs_Parameters.dim1();
  MPI_Bcast(&numCS,1,MPI_INT,0,MPI_COMM_WORLD);

  if(numCS > 0)
    {
      globalInput.cs_Parameters.allocate(numCS);  
      for(int cs=0; cs<numCS; cs++)
	{
	  globalInput.cs_Parameters(cs).allocate(Guess_parameters.dim1());
	  MPI_Bcast(globalInput.cs_Parameters(cs).array(), globalInput.cs_Parameters(cs).dim1(),
		    MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
    }
  else
    {
      globalInput.cs_Parameters.deallocate();
    }

#endif

  globalInput.setAIParameters(Guess_parameters);
  //globalInput.WF.normalizeCI();

  if(globalInput.flags.a_diag > 0 ||
     optStep%2 == 1)
    {
      clog << "(Step = " << setw(3) << optStep << "):    Best Objective Value ";
      clog.precision(12);
      clog.width(20);
      clog.unsetf(ios::scientific);
      clog.unsetf(ios::fixed);
      clog << left << lastRun.energy.getAverage();

      clog << " s^2 = ";
      clog.precision(12);
      clog.width(20);
      clog << left << lastRun.energy.getSeriallyCorrelatedVariance();
      clog << right;

      clog << " ns = ";
      clog.width(10);
      clog << lastRun.energy.getNumberSamples();
      clog << endl << endl;

      /*
      clog.precision(12);
      clog.width(20);
      clog.setf(ios::scientific);
      clog << left << rel << endl << right;
      */
    }

  clog.setf(ios::scientific);
  if(globalInput.flags.a_diag > 0 ||
     optStep%2 == 0)
    {  
      globalInput.printAIParameters(clog,"Best Step params:",20,Guess_parameters,false);
      clog << endl << endl;

      globalInput.JP.print(clog);
    }

  for(int i=0; i<Guess_parameters.dim1(); i++)
    if(IeeeMath::isNaN(Guess_parameters[i]))
      {
	clog << "Error: wavefunction has nans.\n";
	exit(0);
      }
  
  double penalty = globalInput.JP.calculate_penalty_function();
  if(penalty >= 1e10)
    {
      clog << "Error: unable to produce acceptable parameters, quitting..." << endl;
      clog << "  Either your guess for initial Jastrow parameters is bad or you need to change your "
	   << "optimization choices." << endl;
      exit(0);
    }

  optStep++;
}




