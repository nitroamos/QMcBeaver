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
      QMCLineSearchStepLengthSelectionAlgorithm *stepAlg =
	QMCLineSearchStepLengthSelectionFactory::factory(
				       input->flags.line_search_step_length);

      optAlg = new QMCSteepestDescent(&objFunc, stepAlg, 
			 input->flags.optimization_max_iterations, 
			 input->flags.optimization_error_tolerance);

    }  
  else if( input->flags.optimize_Psi_method == "BFGSQuasiNewton" )
    {
      QMCLineSearchStepLengthSelectionAlgorithm *stepAlg =
	QMCLineSearchStepLengthSelectionFactory::factory(
				       input->flags.line_search_step_length);

      optAlg = new QMCBFGSQuasiNewtonLineSearch(&objFunc, stepAlg, 
			 input->flags.optimization_max_iterations, 
			 input->flags.optimization_error_tolerance);

    }
  else if( input->flags.optimize_Psi_method == "CKGeneticAlgorithm1" )
    {      
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
