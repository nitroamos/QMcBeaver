#ifndef QMCGRCOMP_H
#define QMCGRCOMP_H

#include <iostream>
#include <math.h>

using namespace std;

/**
  Either the numerator or the denominator of a Greens function ratio.  The 
  purpose of this is to make the evaluation of ratios of very large or very 
  small numbers more stable.  
  This object has the form k*a^b*exp(c)  
*/

class QMCGreensRatioComponent
{
 public:
  
  /**
    Creates a new uninitialized instance of this class.
  */

  QMCGreensRatioComponent();

  /**
    Creates a new instance of this class with a total value of 'k.'
  */

  QMCGreensRatioComponent(double k);

  /** 
    Creates a new instance of this class.
  */

  QMCGreensRatioComponent(double k, double a, double b, double c);

  /** 
    Creates  new instance of this class and makes it equivalent to another 
    instance of this class.
    @ param rhs object to set this equal to.
  */

  QMCGreensRatioComponent( const QMCGreensRatioComponent & rhs );

  /** 
    Deallocates the memory allocated by this object.
  */

  ~QMCGreensRatioComponent();

  /** 
    Sets two QMCGreensRatioComponent objects equal.
    @ param rhs object to set this object equal to.
  */

  void operator=( const QMCGreensRatioComponent & rhs );

  /**
    @ param X object to add to this object.
    @ return sum of these two GreensRatioComponents 
  */

  QMCGreensRatioComponent operator + ( const QMCGreensRatioComponent & rhs );

  /** 
    @ param denom object to divide this object by
    @ return Greens function ratio of these two GreensRatioComponents.
  */

  double DivideBy(QMCGreensRatioComponent &denom);

  /**
    @ param X object to multiply this object by
    @ return product of these two GreensRatioComponents.
  */

  void MultiplyBy(QMCGreensRatioComponent &X);

  /**
    Writes the state of this object to an XML stream.
    @ param strm XML stream.
  */

  void toXML(ostream & strm);

  /**
    Returns the overall value of the object.
  */

  double getValue();

 private:

  double k,a,b,c;

  /**
    Initializes all of the data members.
  */

  void initialize();

  /**
    na^nb/da^db = ra^rb
  */

  void SimplifyRatioPowers(double na, double nb, double da, double db, 
			   double &ra, double &rb);

};

#endif
