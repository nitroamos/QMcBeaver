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

#ifndef QMCSTATISTIC_H
#define QMCSTATISTIC_H

#include <iostream>
#include <string>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

using namespace std;


/**
  Statistical information on a set of data.
  */

class QMCStatistic
{
private:
  long double sum;
  long double sum2;
  long double weights;
  unsigned long nsamples;

public:
  /**
    Creates a zeroed out instance of the class and generates the MPI type if 
    it has not been done.
    */

  QMCStatistic();


  /**
    Sets all of the data in the object to zero.
    */

  void zeroOut();


  /**
    Gets the number of data samples entered into the object.
    */

  unsigned long getNumberSamples() const;


  /**
    Gets the average of the data entered into the object.
    */

  double getAverage() const;


  /**
    Gets the variance of the data entered into the object.
    */

  double getVariance() const;


  /**
     Gets the standard deviation of the data entered into the object.
  */

  double getStandardDeviation() const;


  /**
    Adds a new data sample to the object.

    @param s new sample data
    @param weight statistical weight of the sample
    */

  void newSample(long double s, long double weight);

  /**
     Sets two object equal.
    */

  void operator = ( const QMCStatistic &rhs);

  /**
    Returns the sum of two QMCStatistics.
    */

  QMCStatistic operator + (const QMCStatistic &rhs);

  /**
		This will multiply sum and sum2 by factor w.
		
		@param w the reweighting factor
    */
  void reWeight(double w);

  /**
    Writes the state of this object to an XML stream.

    @param strm XML stream
    */

  void toXML(ostream& strm);


  /**
    Loads the state of this object from an XML stream.

    @param strm XML stream
    @return whether the read was successful
    */
  bool readXML(istream& strm);


  /**
    Formats and prints the statistic to a stream.
    */
  friend ostream& operator << (ostream& strm, const QMCStatistic & rhs);


#ifdef PARALLEL

private:
  /**
    A flag which tells if MPI_TYPE has been generated.
    */

  static bool mpiTypeCreated;

  
  /**
    Build MPI_TYPE.
    */

  static void buildMpiType();

  /**
    Build MPI_REDUCE.
    */

  static void buildMpiReduce();


  /**
    An MPI function which allows MPI_Reduce to be used in adding QMCStatistics.
    */

  static void Reduce_Function(QMCStatistic *in, QMCStatistic *inout, 
			      int *len, MPI_Datatype *dptr);

public:

  /**
    The MPI data type for a QMCStatistic.
    */
  static MPI_Datatype MPI_TYPE;


  /** 
    The MPI operation for performing MPI_Reduce on QMCStatistics.
    */
  
  static MPI_Op MPI_REDUCE;

#endif

};

#endif
