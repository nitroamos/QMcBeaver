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

#include "CKGeneticAlgorithm1.h"

CKGeneticAlgorithm1::CKGeneticAlgorithm1(QMCObjectiveFunction * function, 
					 int populationsize,
					 double mutationrate,
					 double distributionwidth)
{
  OF                = function;
  PopulationSize    = populationsize;
  MutationRate      = mutationrate;
  InitialDistributionWidth = distributionwidth;

  iseed = 102304;
}

void CKGeneticAlgorithm1::initializePopulation(Array1D<double> &initialguess)
{
  Array1D < Array1D<double> > ParamArray(PopulationSize);

  // Generate the new random population
  for(int i=0; i<PopulationSize; i++)
    {
      Array1D<double> Parameters( initialguess.dim1() );

      // Add a random number to the initial guess
      for(int j=0; j<initialguess.dim1(); j++)
	{
	  Parameters(j) = initialguess(j) + 
	    InitialDistributionWidth*gasdev( &iseed );
	}
      
      ParamArray(i) = Parameters;
    }

  // Evaluate the scores for the population
  Array1D <QMCObjectiveFunctionResult> OFresults = OF->evaluate( ParamArray );

  // Add the parameters and scores to the population
  Population.clear();
  for(int i=0; i<PopulationSize; i++)
    {
      ParameterScorePair PSP(OFresults(i).getScore(),
			     ParamArray(i) );

      Population.add( PSP );
    }
}

void CKGeneticAlgorithm1::mutate(Array1D<double>& parameters)
{
  for(int i=0; i<parameters.dim1(); i++)
    {
      parameters(i) += MutationRate*gasdev( &iseed );
    }
}

Array1D<double> CKGeneticAlgorithm1::crossover(Array1D<double> &Parent1, 
					       Array1D<double> &Parent2)
{
  Array1D<double> result(Parent1.dim1());

  for(int i=0; i<result.dim1(); i++)
    {
      double alpha = ran1( &iseed );
      result(i)    = alpha * Parent1(i) + (1.0-alpha) * Parent2(i);
    }

  return result;
}

int CKGeneticAlgorithm1::selectParent()
{
  // Weight each population member so that the best population member
  // gets a weight of PopulationSize and the worst population member
  // gets a weight of 1 with a line inbetween.  This is arbitrary and could
  // be improved

  // To do this in the continuous case, the probability density can
  // be inverted so that # = PopulationSize * (1-sqrt(1-ran(0,1)))
  // for our case # = int( PopulationSize * (1-sqrt(1-ran(0,1))) );

  return int( PopulationSize * ( 1-sqrt( 1-ran1(&iseed) ) ) );
}

void CKGeneticAlgorithm1::generateNewPopulation()
{
  // Generate the new parameter vectors
  // Keep the best member of the previous population

  Array1D < Array1D<double> > NewParameters(PopulationSize);

  NewParameters(0) = *Population.get( 0 ).getParameters();

  for(int i=1; i<PopulationSize; i++)
    {
      ParameterScorePair Parent1 = Population.get( selectParent() );
      ParameterScorePair Parent2 = Population.get( selectParent() );

      // Form the child with corssing over
      Array1D<double> Child = crossover( *Parent1.getParameters(), 
					 *Parent2.getParameters());

      // Mutate the child
      mutate( Child );

      NewParameters(i) = Child;
    }

  // Evaluate the scores for the population
  Array1D <QMCObjectiveFunctionResult> OFresults = 
    OF->evaluate( NewParameters );

  // Add the parameters and scores to the population
  Population.clear();
  for(int i=0; i<PopulationSize; i++)
    {
      ParameterScorePair PSP(OFresults(i).getScore(),
			     NewParameters(i) );

      Population.add( PSP );
    }
}

Array1D<double> CKGeneticAlgorithm1::optimize(Array1D<double> & InitialGuess)
{
  initializePopulation(InitialGuess);

  // Initialize the best solution
  BestPopulationMember = Population.get(0);

  bool done = false;
  int step = 0;
  while( !done )
    {
      step++;

      // Generate a new population
      generateNewPopulation();

      // Update the best population member
      if( Population.get(0) < BestPopulationMember )
	{
	  BestPopulationMember = Population.get(0);
	}

      // Output the results from this iteration
      cout << "Step: " << step << "\tBest Score: " 
	   << BestPopulationMember.getScore() << endl;

      // see if done
      if( step > 30 )
	{
	  done = true;
	}
    }

  return *BestPopulationMember.getParameters();
}















