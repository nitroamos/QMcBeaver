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

//DCL=15 will decorrelate us to around 30 thousand steps
//DCL=20 will decorrelate us to around 1 million steps
//DCL = Decorrelation Length in log2(n)

#ifndef QMCPROPERTY_H
#define QMCPROPERTY_H

#define DCL 30

#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "QMCStatistic.h"
#include "Array1D.h"
#include "Array2D.h"
#include "LU.h"

using namespace std;

/**
  All of the statistical information used in calculating a quantity or 
  property during a calculation.
*/

class QMCProperty
{
private:
  QMCStatistic DeCorr[DCL];
  int DeCorr_flags[DCL];
  double DeCorr_sample[DCL];
  double DeCorr_weight[DCL];

  //static const int numberParams = 4;
  ///**
  //   Parameters used in fitting the decorrelation plot to a function to
  //   estimate the standard deviation.
  //*/
  //double stdevFittingParameters[numberParams];

  /**
    Gets the decorrelation depth- the index of the decorrelation array where
    the variance reaches a plateau.
  */
  int getDecorrDepth();

public:

  /**
    Creates a zeroed out instance of the class and generates the MPI types if 
    they have not been done.
  */
  QMCProperty();

  /**
    Sets all of the data in the object to zero.
  */
  void zeroOut();

  /**
    Adds a new data sample to the object.
    @param s new sample data
    @param weight statistical weight of the sample
  */
  void newSample(double s, double weight);

  /**
    Gets the number of data samples entered into the object.
    @return number of samples in the object.
  */
  unsigned long getNumberSamples();

  /**
    Gets the average of the data entered into the object.
    @return average of the data in the object.
  */
  double getAverage();

  /**
    Gets the variance of the data entered into the object.
    @return variance of the data in the object.
  */
  double getVariance();

  /**
    Gets the serially correlated variance of the data entered into the object.
    @return serially correlated variance of the data in the object.
  */
  double getSeriallyCorrelatedVariance();

  /**
    Gets the standard deviation of the data entered into the object.
    @return standard deviation of the data in the object.
  */
  double getStandardDeviation();

  /**
    Gets the serially correlated standard deviation of the data entered 
    into the object.
    @return serially correlated standard deviation of the data in the object.
  */
  double getSeriallyCorrelatedStandardDeviation();

  /**
    Gets the standard deviation of the standard deviation.
    @return standard deviation of the standard deviation.
  */
  double getStandardDeviationStandardDeviation();

  /**
     Sets two objects equal.
  */
  void operator = ( const QMCProperty &rhs);

	/**
		This will add the data from rhs to our data
	*/
  void operator += ( QMCProperty &rhs);

  /**
    Returns the sum of two QMCproperties.
    @param rhs QMCProperty to add to this one.
    @return sum of these two QMCProperties.
  */
  QMCProperty operator + ( QMCProperty &rhs);

	/**
		This will change the overall importance
		of our data. This function basically just
		calls reWeight on the QMCStatistic DeCorr
		objects.
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
    Formats and prints the all the decorrelation data 
    for the property to a stream.
  */
  void printAll(ostream & strm);

  /**
    Formats and prints the property to a stream as a single
    line.
  */
  friend ostream& operator <<(ostream& strm, QMCProperty &rhs);

  /**
    Gets the correct variance for the block.  Where the variance is 
    calculated using
    \f[
    \sigma^2 = \frac{\left(\frac{1}{n}\sum^{n}_{i=1}x^{2}_{i}\right) -
    \left(\frac{1}{n}\sum^{n}_{i=1}x_{i}\right)^{2}}{n-1}
    \f]
    as is described in the Dynamic Distributable Decorrelation Algorithm 
    (DDDA) paper by Kent and Feldmann.
    @param i block number.
    @return variance for block.
  */
  double getBlockVariance(int i);

  /**
    Gets the correct standard deviation for the block.  Where the standard
    deviation is calculated using
    \f[
    \sigma = \sqrt{\frac{\left(\frac{1}{n}\sum^{n}_{i=1}x^{2}_{i}\right) -
    \left(\frac{1}{n}\sum^{n}_{i=1}x_{i}\right)^{2}}{n-1}}
    \f]
    as is described in the Dynamic Distributable Decorrelation Algorithm 
    (DDDA) paper by Kent and Feldmann.
    @param i block number.
    @return standard deviation for block.
  */
  double getBlockStandardDeviation(int i);

  /**
    Gets the correct standard deviation of the standard deviation for the 
    block.  Where the standard deviation of the standard deviation is 
    calculated using
    \f[
    \sigma = \sqrt{\frac{\left(\frac{1}{n}\sum^{n}_{i=1}x^{2}_{i}\right) -
    \left(\frac{1}{n}\sum^{n}_{i=1}x_{i}\right)^{2}}{n-1}}
    \left(\frac{1}{\sqrt{2(n-1)}}\right)
    \f]
    as is described in the Flyvbjerg-Petersen data blocking paper 
    (JCP 1989 p461).
    @param i block number.
    @return standard deviation of the standard deviation for block.
  */
  double getBlockStandardDeviationStandardDeviation(int i);

  /**
    Gets the correct standard deviation of the variance for the 
    block.  Where the standard deviation of the variance is 
    calculated using
    \f[
    \sigma = \frac{\left(\frac{1}{n}\sum^{n}_{i=1}x^{2}_{i}\right) -
    \left(\frac{1}{n}\sum^{n}_{i=1}x_{i}\right)^{2}}{n-1}
    \sqrt{\frac{2}{n-1}})
    \f]
    as is described in the Flyvbjerg-Petersen data blocking paper 
    (JCP 1989 p461).
    @param i block number.
    @return standard deviation of the standard deviation for block.
  */
  double getBlockVarianceStandardDeviation(int i);

 private:

  /**
    Calculates the value and gradient of the function used to fit the
    standard deviation.  This function is
    \f[
    r^2(p) = \sum\left[ y_{i} - f(x_{i};p) \right]^2
    \f]
    where
    \f[
    f(x_{i};p) = \frac{p(0)^{2}e^{p(2)^{2}x+p(3)^{2}x^{2}}}
    {1+p(1)^{2}e^{p(2)^{2}x+p(3)^{2}x^{2}}}
    \f]
    and
    \f[
    x_{i} = log_{2}(blockSize_{i}).
    \f]
    @param params set of parameters to evaluate the function with
    @param standardDeviations the standard deviation calculated for 
    each block size.  These could be calculated on the fly but have 
    been added here to improve performance.
    @param standardDeviationsErrors one standard deviation error in the 
    standard deviation calculated for each block size.  These could be
    calculated on the fly but have been added here to improve performance.
    @param functionValue calculated function value is returned here
    @param gradientValue calculated gradient value is returned here
  */
  static void calculateObjectiveFunction(Array1D<double> & params, 
                                    Array1D<double> & standardDeviations,
				    Array1D<double> & standardDeviationsErrors,
				    double & functionValue, 
				    Array1D<double> & gradientValue);

  /**
    Calculates the value of the function used to fit the
    standard deviation.  This function is
    \f[
    r^2(p) = \sum\left[ y_{i} - f(x_{i};p) \right]^2
    \f]
    where
    \f[
    f(x_{i};p) = \frac{p(0)^{2}e^{p(2)^{2}x+p(3)^{2}x^{2}}}
    {1+p(1)^{2}e^{p(2)^{2}x+p(3)^{2}x^{2}}}
    \f]
    and
    \f[
    x_{i} = log_{2}(blockSize_{i}).
    \f]
    @param params set of parameters to evaluate the function with
    @param standardDeviations the standard deviation calculated for 
    each block size.  These could be calculated on the fly but have 
    been added here to improve performance.
    @param standardDeviationsErrors one standard deviation error in the 
    standard deviation calculated for each block size.  These could be
    calculated on the fly but have been added here to improve performance.
    @param functionValue calculated function value is returned here
  */
  static void calculateObjectiveFunction(Array1D<double> & params, 
				    Array1D<double> & standardDeviations, 
				    Array1D<double> & standardDeviationsErrors,
				    double & functionValue);

  /**
    Calculates the value of the 1-D objective function from the line search.
    This function is 
    \f[
    \phi(\alpha) = r^2(p+\alpha d)
    \f]
    where \f$p\f$ is the current set of parameters, \f$d\f$ is the search
    direction, and \f$\alpha\f$ is the step length parameter. 
    @param params set of parameters to evaluate the function with
    @param searchDirection direction the line search is being performed along.
    @param stepLength parameter determining how long of a step to take in the
    line search.
    @param standardDeviations the standard deviation calculated for 
    each block size.  These could be calculated on the fly but have 
    been added here to improve performance.
    @param standardDeviationsErrors one standard deviation error in the 
    standard deviation calculated for each block size.  These could be
    calculated on the fly but have been added here to improve performance.
    @return value of the 1-D objective function for the given step length 
    parameter.
    */
  static double calculateLineSearchObjectiveFunction(Array1D<double> & params,
				  Array1D<double> & searchDirection,
				  double stepLength,
				  Array1D<double> & standardDeviations,
				  Array1D<double> & standardDeviationsErrors); 

  /**
    Calculates the derivative of the 1-D objective function from the line 
    search.  This function is 
    \f[
    \phi(\alpha) = r^2(p+\alpha d)
    \f]
    where \f$p\f$ is the current set of parameters, \f$d\f$ is the search
    direction, and \f$\alpha\f$ is the step length parameter.
    @param params set of parameters to evaluate the function with
    @param searchDirection direction the line search is being performed along.
    @param stepLength parameter determining how long of a step to take in the
    line search.
    @param standardDeviations the standard deviation calculated for 
    each block size.  These could be calculated on the fly but have 
    been added here to improve performance.
    @param standardDeviationsErrors one standard deviation error in the 
    standard deviation calculated for each block size.  These could be
    calculated on the fly but have been added here to improve performance.
    @return derivate of the 1-D objective function for the given step length 
    parameter.
   */
  static double calculateLineSearchObjectiveFunctionDerivative(
				  Array1D<double> & params,
				  Array1D<double> & searchDirection,
				  double stepLength,
				  Array1D<double> & standardDeviations,
				  Array1D<double> & standardDeviationsErrors); 

  /**
    Perform a cubic interpolation to choose a step length in the 
    interval [a_lo,a_hi]. See Nocedal and Wright p 56-58.
    @param alphaLo lower bound on the step length.
    @param alphaHi upper bound on the step length.
    @param phi_0 objective function value at the current parameters
    @param phi_alphaLo 1-D objective function value at the smallest step 
    length.
    @param phi_alphaHi 1-D objective function value at the largest step length.
    @param phiPrime_0 derivative of the 1-D objective function at the
    current parameters.
    @return calculated step length 
   */
  static double cubicInterpolateStep(double alphaLo, double alphaHi, 
				     double phi_0, 
				     double phi_alphaLo, double phi_alphaHi, 
				     double phiPrime_0);

  /**
    Finds a point in the interval [a_lo, a_hi] that satisfy the strong Wolfe
    conditions.  See Nocedal and Wright p 60.
    @param alphaLo lower bound on the step length.
    @param alphaHi upper bound on the step length.
    @param phi_0 objective function value at the current parameters
    @param phi_alphaLo 1-D objective function value at the smallest step 
    length.
    @param phi_alphaHi 1-D objective function value at the largest step length.
    @param phiPrime_0 derivative of the 1-D objective function at the
    current parameters.
    @param params set of parameters to evaluate the function with
    @param searchDirection direction the line search is being performed along.
    @param standardDeviations the standard deviation calculated for 
    each block size.  These could be calculated on the fly but have 
    been added here to improve performance.
    @param standardDeviationsErrors one standard deviation error in the 
    standard deviation calculated for each block size.  These could be
    calculated on the fly but have been added here to improve performance.
    @return step length satisfying the strong Wolfe conditions.
   */
  static double zoom(double alphaLo, double alphaHi, double phi_0, 
		     double phi_alphaLo, double phi_alphaHi, 
		     double phiPrime_0, Array1D<double> & params,
		     Array1D<double> & searchDirection,
		     Array1D<double> & standardDeviations,
		     Array1D<double> & standardDeviationsErrors);

  /**
    Finds a point between 0 and the maximum step length which satisfies the
    Wolfe conditions.  See Nocedal and Wright p 59.
    @param alphaGuess guess at a step length satisfying the Wolfe conditions.
    @param params set of parameters to evaluate the function with
    @param searchDirection direction the line search is being performed along.
    @param gradient gradient of the objective function.
    @param functionValue value of the objective function at the current 
    set of parameters.
    @param standardDeviations the standard deviation calculated for 
    each block size.  These could be calculated on the fly but have 
    been added here to improve performance.
    @param standardDeviationsErrors one standard deviation error in the 
    standard deviation calculated for each block size.  These could be
    calculated on the fly but have been added here to improve performance.
    @return step length satisfying the strong Wolfe conditions.
   */
  static double wolfeStepLength(double alphaGuess, Array1D<double> & params,
				Array1D<double> & searchDirection,
				Array1D<double> & gradient,
				double functionValue,
				Array1D<double> & standardDeviations,
				Array1D<double> & standardDeviationsErrors);
 
  /**
    Generates an initial guess set of parameters used in fitting the 
    standard deviation decorrelation plot.
  */
  void generateInitialGuessFittingParameters();

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
    An MPI function which allows MPI_Reduce to be used in adding QMCproperties.
  */
  static void Reduce_Function(QMCProperty *in, QMCProperty *inout, 
                              int *len, MPI_Datatype *dptr);

public:

  /**
    The MPI data type for a QMCProperty.
    */
  static MPI_Datatype MPI_TYPE;

  /** 
    The MPI operation for performing MPI_Reduce on QMCproperties.
    */
  static MPI_Op MPI_REDUCE;

#endif

};

#endif
