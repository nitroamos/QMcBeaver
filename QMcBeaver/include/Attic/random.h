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

#ifndef RANDOM_H
#define RANDOM_H

/**
  @file random.h
  Library of functions for generating random numbers.
  */

/**
  Generates a uniform random number on \f$[0,1]\f$ using the ran1 algorithm 
  from numerical recipes.

  @param idum random number seed
  @return uniform random number on \f$[0,1]\f$.
  */
double ran1(long *idum);

/**
  Generates a gaussian distributed random number with unit variance using the
  gasdev algorithm from numerical recipes.  To generate a gaussian random
  number with another variance, multiply this result by the desired standard
  deviation.

  @param idum random number seed
  @return gaussian random number with unit variance
  */

double gasdev(long *idum);


/**
  Generates an exponentially distributed random number whith the probability
  distribution function
  \f[
  p(x)dx = e^{-x}dx.
  \f]
  To generate points according to 
  \f[
  p(x)dx = a e^{-a*x}dx 
  \f]
  divide the result of this function by \f$a\f$.
  The random variable is a number on the positive real axis.

  This method was developed by David Randall ``Chip'' Kent IV and works by
  directly inverting the distribution.

  @param idum random number seed
  @return exponential random number
  */

double expdev(long *idum);

/**
  Generates a random number distributed with the probability
  distribution function
  \f[
  p(x)dx = \frac{1}{2}\sin(x)dx.
  \f]
  The random variable is a number in \f$[0,\Pi]\f$

  This method was developed by David Randall ``Chip'' Kent IV and works by
  directly inverting the distribution.

  @param idum random number seed
  @return random number distributed with respect to \f$\sin(x)/2\f$.
  */

double sindev(long *idum);

/**
  Generates a random number distributed whith the probability distribution 
  function
  \f[
  p(x)dx = \frac{1}{2}x^{2}e^{-x}dx.
  \f]
  The random variable is a number on the positive real axis.

  This method was developed by David Randall ``Chip'' Kent IV and works by
  using a rejection type algorithm.

  @param idum random number seed
  @return random number distributed with respect to \f$x^{2}e^{-x}/2\f$.
  */

double randomDistribution1(long *idum);



#endif
