//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2002.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

#include "QMCOptimizationFactory.h"

QMCOptimizationAlgorithm * QMCOptimizationFactory::
    optimizationAlgorithmFactory(QMCObjectiveFunction &objFunc, 
				 QMCInput * input)
{
  QMCOptimizationAlgorithm * optAlg = 0;

  if( input->flags.optimize_Psi_method == "Steepest_Descent" )
    {
      /*
	This method will use the identity matrix for its hessian. It's fairly
	reliable, given good derivatives. However, it will probably take you more
	steps to converge because the best parameters are usually not in the Steepest_Descent
	direction.
      */
      QMCLineSearchStepLengthSelectionAlgorithm *stepAlg =
	QMCLineSearchStepLengthSelectionFactory::factory(
				       input->flags.line_search_step_length);

      optAlg = new QMCSteepestDescent(&objFunc, stepAlg, 
			 input->flags.optimization_max_iterations, 
			 input->flags.optimization_error_tolerance);

    }  
  else if( input->flags.optimize_Psi_method == "analytical_energy_variance" ||
	   input->flags.optimize_Psi_method == "automatic")
    {
      /*
	The analytical_energy_variance method will evaluate a hessian. See the
	the PRL 94, 150201 (2005) paper for the formulas.

	This hessian will provide very bad guesses if the initial parameters
	are far from the minimum well, since the Jastrow parameters are nonlinear.
	However, once we're well inside the well, the method should have quadratic
	convergence.

	To compensate, I've added an "automatic" method, which will use
	Steepest_Descent for a couple of steps, and then automatically switch
	to the hessian. I might have the automatic method automatically increase
	the VMC iterations as the optimization progresses.

	I expect the "automatic" method to work the best.
      */
      QMCLineSearchStepLengthSelectionAlgorithm *stepAlg =
	QMCLineSearchStepLengthSelectionFactory::factory(
				       input->flags.line_search_step_length);

      optAlg = new QMCSteepestDescent(&objFunc, stepAlg, 
			 input->flags.optimization_max_iterations, 
			 input->flags.optimization_error_tolerance);

    }  
  else if( input->flags.optimize_Psi_method == "BFGSQuasiNewton" )
    {
      /*
	This method approximates a hessian to use with approximated derivatives.
      */
      QMCLineSearchStepLengthSelectionAlgorithm *stepAlg =
	QMCLineSearchStepLengthSelectionFactory::factory(
				       input->flags.line_search_step_length);

      optAlg = new QMCBFGSQuasiNewtonLineSearch(&objFunc, stepAlg, 
			 input->flags.optimization_max_iterations, 
			 input->flags.optimization_error_tolerance);

    }
  else if( input->flags.optimize_Psi_method == "CKGeneticAlgorithm1" )
    {      
      /*
	This method does not use any hessians or derivatives. It will
	randomly pick parameters given a standard deviation, and evaluate
	a score for that parameter set. The parameters with better scores
	will then be propagated into successive generations using a genetic
	algorithm.

	The method should be very robust, even if it is slow.
      */
      optAlg = new CKGeneticAlgorithm1(&objFunc, 
         input->flags.ck_genetic_algorithm_1_population_size, 
	 input->flags.ck_genetic_algorithm_1_mutation_rate,
         input->flags.ck_genetic_algorithm_1_initial_distribution_deviation);
    }
  else
    {
      cerr << "ERROR: Unknown optimization algorithm ("
	   << input->flags.optimize_Psi_method
	   << ") in QMCOptimizationFactory!" << endl;
      exit(1);
    }

  return optAlg;
}
