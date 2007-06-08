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
}

void CKGeneticAlgorithm1::initializePopulation(Array1D<double> &initialguess)
{
  Array1D < Array1D<double> > ParamArray(PopulationSize);

  //make sure that we have our best known result included
  ParamArray(0) = initialguess;

  // Generate the rest of the population randomly
  for(int i=1; i<ParamArray.dim1(); i++)
    {
      Array1D<double> Parameters = initialguess;
      mutate(Parameters, InitialDistributionWidth);      
      ParamArray(i) = Parameters;
    }

  // Evaluate the scores for the population

  sw.reset();
  sw.start();
  Array1D <QMCObjectiveFunctionResult> OFresults = OF->evaluate( ParamArray );
  sw.stop();

  // Add the parameters and scores to the population
  Population.clear();
  for(int i=0; i<ParamArray.dim1(); i++)
    {
      ParameterScorePair PSP(OFresults(i),ParamArray(i) );

      Population.add( PSP );
    }
}

void CKGeneticAlgorithm1::mutate(Array1D<double>& parameters, double width)
{
  /*
    If the parameter ends up being too small... then it won't be any good
    anyway, so we can just choose to ignore them here.
   */
  for(int i=0; i<parameters.dim1(); i++)
    {
      double change;
      do {
	change = width*ran.gasdev();
      } while(parameters(i) + change <= 5e-2);
      parameters(i) += change;
    }
}

Array1D<double> CKGeneticAlgorithm1::crossover(Array1D<double> &Parent1, 
					       Array1D<double> &Parent2)
{
  Array1D<double> result(Parent1.dim1());

  for(int i=0; i<result.dim1(); i++)
    {
      double alpha = ran.unidev();
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

  return int( Population.size() * ( 1-sqrt( 1-ran.unidev() ) ) );
}

void CKGeneticAlgorithm1::generateNewPopulation()
{
  // Generate the new parameter vectors
  // Keep the best several from the previous population

  Array1D < Array1D<double> > NewParameters(PopulationSize);

  for(int i=0; i<PopulationSize; i++)
    {
      ParameterScorePair Parent1 = Population.get( selectParent() );
      ParameterScorePair Parent2 = Population.get( selectParent() );

      // Form the child with crossing over
      Array1D<double> Child = crossover( *Parent1.getParameters(), 
					 *Parent2.getParameters());

      // Mutate the child
      /*
	Just to speed up any parameter drift, let's throw in
	a couple with higher distribution widths.
       */
      int numHigh = (int)(0.1 * PopulationSize);
      if(i < numHigh) mutate( Child, InitialDistributionWidth );
      else            mutate( Child, MutationRate );

      NewParameters(i) = Child;
    }

  // Evaluate the scores for the population
  sw.reset();
  sw.start();
  Array1D <QMCObjectiveFunctionResult> OFresults = 
    OF->evaluate( NewParameters );
  sw.stop();

  // Add the parameters and scores to the population
  for(int i=0; i<PopulationSize; i++)
    {
      ParameterScorePair PSP(OFresults(i),NewParameters(i) );

      Population.add( PSP );
    }

  /*
    Since the list is automatically sorted, this will
    save the best PopulationSize known results.
    
    If we save too many, efficiency is lost
    when selecting parents.
  */
  int numToKeep = min(100,PopulationSize);
  Population.resize(numToKeep);
}

Array1D<double> CKGeneticAlgorithm1::optimize(Array1D<double> & InitialGuess,
					      double value,
					      Array1D<double> & gradient,
					      Array2D<double> & hessian,
					      double misc)
{
  initializePopulation(InitialGuess);

  bool done = false;
  int step = 0;

  while( !done )
    {
      step++;

      // Generate a new population
      // after we've shown the results from the initialization
      if(step > 1)
	generateNewPopulation();
      
      cout << "\nIteration: " << step << " took (" << sw << ")" << endl;

      /*
	Print out some of the best so far, so we can see
	how much agreement there is.
      */
      int numToShow = 10;
      if(numToShow > min(100,PopulationSize))
	numToShow = min(100,PopulationSize);

      for(int i=0; i<numToShow; i++)
	{
	  printf("%3i) ",i);
	  cout << Population.get(i);
	}

      // see if done
      if( step > 3 )
	{
	  done = true;
	}
    }

  /*
    The Population array is automatically sorted, so the 0th
    element is always the best available.
   */
  return *Population.get(0).getParameters();
}

