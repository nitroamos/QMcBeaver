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

#ifndef CKGENETICALGORITHM1_H
#define CKGENETICALGORITHM1_H

#include "SortedParameterScorePairList.h"
#include "QMCObjectiveFunction.h"
#include "QMCObjectiveFunctionResult.h"
#include "QMCOptimizationAlgorithm.h"
#include "Random.h"
#include "Stopwatch.h"

/**
  A moderately greedy genetic algorithm for trying to globally optimize
  a function dreamed up by David Randall (Chip) Kent IV.  As is standard
  in the field, optimization means minimization.

  Mutation is accomplished by adding a N-dimensional gaussian random variable
  to the population member.

  The amount of each parent contributed to a child is determined by a uniform
  random variable.

  A linear probability distribution is used to select which population member
  will be a parent.  The best members have better probabilities of being
  selected.
  */

class CKGeneticAlgorithm1 : public QMCOptimizationAlgorithm
{
public:
  /**
    Constructs and inializes this optimization algorithm.

    @param function function to optimize.
    @param populationsize number of members in the population used to 
    optimize the function.  This is a positive number.
    @param mutationrate a positive number describing how much mutation is 
    introduced into the population.  Larger numbers correspond to more
    mutation.
    @param distributionwidth a positive number describing how far the initial
    population members spread from the initial guess.
    */
  CKGeneticAlgorithm1(QMCObjectiveFunction * function, 
		      int populationsize,
		      double mutationrate,
		      double distributionwidth);

  Array1D<double> optimize(Array1D<double> & initialGuess,
			   double value,
			   Array1D<double> & gradient,
			   Array2D<double> & hessian,
			   double);

private:
  Stopwatch sw;

  /**
    Objective function to optimize.
    */
  QMCObjectiveFunction* OF;

  /**
    Rate that there are mutations in the population.
    */
  double MutationRate;

  /**
    Number of members in the population.
    */
  int PopulationSize;

  /**
    Spread of the initial population around the initial guess.
    */
  double InitialDistributionWidth;

  /**
    Population ordered so that the better parameter sets have smaller indices.
    */
  SortedParameterScorePairList Population;

  /**
    Best parameter set found.
    */
  ParameterScorePair BestPopulationMember;

  /**
    Mutates the input parameters by adding an N-dimensional gaussian 
    random variable with a standard deviation of <code>MutationRate</code>.
    
    It will make sure that no parameter goes below zero.

    @param parameters parameters to be mutated and the returned mutated 
    parameters.
    @param width the width of the gaussian drift applied to each parameter
    */
  void mutate(Array1D<double>& parameters, double width);

  /**
    Generates a child by breeding two parents.  The fraction of each
    parent is chosen with a uniform random distribution.

    @param parent1 first partent of the child.
    @param parent2 second parent of the child.
    */
  Array1D<double> crossover(Array1D<double> &parent1, 
			    Array1D<double> &parent2);

  /**
    Selects a parent for breeding.  The probability distribution of 
    parent <code>i</code> being selected is a linear distribution with the
    first member having a weight of <code>PopulationSize</code> and the
    last member having a weight of <code>1</code>. 
    \f[
    \rho(member_{i}) = 2\frac{PopulationSize-i}{PopulationSize^{2}}
    \f]
    The details of this distriubution are worked out for the continuous
    case and used without modification here.

    @return index of the selected parent.
    */
  int selectParent();

  /**
    Breeds and mutates the current population to generate a new population.
    */
  void generateNewPopulation();

  /**
    Create the initial population of points surrounding the initial guess.
    The points are distributed with respect to an N-dimensional gaussian
    distribution with a standard deviation of 
    <code>InitialDistributionWidth</code>.

    @param initialguess starting point for the optimization.
    */
  void initializePopulation(Array1D<double> &initialguess);
};

#endif
