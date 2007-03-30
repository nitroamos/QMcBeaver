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

#include <iostream>
#include <string>
#include <math.h>
#include <limits>

#ifdef USESPRNG
#include "sprng_cpp.h"
#endif

#ifdef PARALLEL
#include "mpi.h"
#endif

#define NTAB 32

using namespace std;
/**
   @file random.h
   Library of functions for generating random numbers.

   SPRNG is available at
   http://sprng.cs.fsu.edu/

   You do not need to compile the MPI version of SPRNG. The
   only difference in the libraries is that it makes sure
   make_sprng_seed produces the same seed on all processes;
   which is unnecessary for our purposes.
*/

class Random
{
 private:
  /**
     Store the initial value for the seed
  */
  long start;

  /**
     When using ran1, this represents the current
     state of the initial seed.
  */
  long current;

#ifdef USESPRNG
  /**
     The identifier for our SPRNG stream.
  */
  Sprng * stream;

#endif

  //variables used by ran1
  long iy;
  long iv[NTAB];

  /**
     Generates a uniform random number on \f$[0,1]\f$ using the ran1 algorithm 
     from numerical recipes.

     "Minimal" random number generator of Park and Miller with Bays-Durham shuffle
     and added safeguards. Returns a uniform random deviate between 0.0 and 1.0, exclusive
     of the end point values. Call with *idum a negative number to initialize, and thereafter
     do not alter *idum between sucessive deviates in a sequence. RNMX should approximate
     the largest floating value that is less than 1.
     
     @param idum random number seed
     @return uniform random number on \f$[0,1]\f$.
  */
  double ran1(long *idum);

 public:

  Random();

  ~Random();

  Random(long seed);

  Random(const Random & rhs);

  /**
     This is called to initialize a Random object.
     
     @param seed is the seed used for the stream of random numbers
     @param rank is used to generate a different seed for each processor,
                 ignored if SPRNG is available.
  */
  void initialize(long seed, int rank);

  /**
     Print the current status of the random number generator.

     @param strm A strm to print to.
  */
  void printStream(ostream & strm);

  /**
     Checkpoint the Random object
     
     @param strm where to write the checkpoint data
  */
  void writeXML(ostream & strm);

  /**
     Recover a checkpointed Random object.

     @param strm where the stream info is available to be read.
  */
  void readXML(istream & strm);

  /**
     This will generate a random integer with a uniform
     deviation. If the SPRNG library is available, it will
     be used.
     
     @param idum random number seed
     @return uniform random number on \f$[0,INT_MAX]\f$.
  */
  int intdev();

  /**
     This will generate a random number with a uniform
     deviation. If the SPRNG library is available, it will
     be used.
     
     @param idum random number seed
     @return uniform random number on \f$[0,1]\f$.
  */
  double unidev();
    
  /**
     Generates a gaussian distributed random number with unit variance using the
     gasdev algorithm from numerical recipes.  To generate a gaussian random
     number with another variance, multiply this result by the desired standard
     deviation.
     
     @param idum random number seed
     @return gaussian random number with unit variance
  */
  double gasdev();
  
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
  double expdev();

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
  double sindev();
  
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
  
  double randomDistribution1();
};  

/**
   We'll maintain one Random object for all other files
   to use. The object lives in QMCManager.
*/
extern Random ran;

#endif
