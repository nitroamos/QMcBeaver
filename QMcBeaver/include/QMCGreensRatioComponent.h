#ifndef QMCGRCOMP_H
#define QMCGRCOMP_H

#include <iostream>
#include <math.h>
#include "IeeeMath.h"

using namespace std;

/**
  Either the numerator or the denominator of a Greens function ratio.  The 
  purpose of this is to make the evaluation of ratios of very large or very 
  small numbers more stable.  
  This object has the form k*a^b*exp(c)

  Notes from Amos:
  This class has been extended a bit from its Green's function only usage.
  This is now intended to be a general purpose class to hold and handle super
  large or super small numbers. In particular, it is now used in 
  QMCFunctions to handle super small psi multiplied by super huge jastrow term.
  
  The full set of math operators have been defined including a cast to double
  operator. One word of warning, if GRC is a QMCGreensRatioComponent object,
  then -1.0*GRC is not necessarily equal to GRC*-1.0. This is because the 
  latter calls the QMCGreensRatioComponent operator * (double) function, while 
  the former casts GRC to a double (which *might* be a problem), then 
  multiplies by -1.0. 

  Perhaps in the future we will want to make this class even better. 
  Specifically, include a "balancing" function to shift large exponents out of
  the k degree of freedom and into the c degree of freedom. Also maybe add 
  checks in the divideBy function.

  Some class philosophy:
  I sort of intend this class to be a data type... an ExtendedDouble type, with
  all the appropriate overloaded operators.
  It now has 3 different classes of math operators:
  1) Ones that return by reference the object operated on (divideBy, 
  multiplyBy, add)
  2) Ones that return nothing (the +=, *=, /=, etc)
  3) Ones that return a new QMCGreensRatioComponent, leaving the operands 
  unaffected (*, +, /, etc)
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
     @param k total value for the new object.
  */
  QMCGreensRatioComponent(double k);
  
  /** 
      Creates a new instance of this class of the form k*a^b*c.
      @param k coefficient
      @param a base
      @param b power
      @param c exponent
  */
  QMCGreensRatioComponent(double k, double a, double b, double c);
  
  /** 
      Creates  new instance of this class and makes it equivalent to another 
      instance of this class.
      @param rhs object to set this equal to.
  */
  QMCGreensRatioComponent( const QMCGreensRatioComponent & rhs );
  
  /** 
      Deallocates the memory allocated by this object.
  */
  ~QMCGreensRatioComponent();
  
  /**
     Divides this QMCGreensRatioComponent by another one and returns a double.
     @param denom object to divide this object by
     @return *this
  */
  QMCGreensRatioComponent & divideBy(const QMCGreensRatioComponent &denom);
  
  /**
     Multiplies two QMCGreensRatioComponents together.
     @param *this
  */
  QMCGreensRatioComponent & multiplyBy(const QMCGreensRatioComponent &rhs);
  
  /**
     Adds two QMCGreensRatioComponents together.
     @param *this
  */
  QMCGreensRatioComponent & add(const QMCGreensRatioComponent &rhs);
  
  
  QMCGreensRatioComponent operator + ( const QMCGreensRatioComponent & rhs ) const;
  QMCGreensRatioComponent operator - ( const QMCGreensRatioComponent & rhs ) const;
  QMCGreensRatioComponent operator * ( const QMCGreensRatioComponent & rhs ) const;
  QMCGreensRatioComponent operator / ( const QMCGreensRatioComponent & rhs ) const;
  void operator += ( const QMCGreensRatioComponent & rhs );
  void operator -= ( const QMCGreensRatioComponent & rhs );  
  void operator *= ( const QMCGreensRatioComponent & rhs );
  void operator /= ( const QMCGreensRatioComponent & rhs );
  void operator =  ( const QMCGreensRatioComponent & rhs );
  
  QMCGreensRatioComponent operator + ( const double & rhs ) const;
  QMCGreensRatioComponent operator - ( const double & rhs ) const;
  QMCGreensRatioComponent operator * ( const double & rhs ) const;
  QMCGreensRatioComponent operator / ( const double & rhs ) const;
  void operator += ( const double & rhs );
  void operator -= ( const double & rhs );  
  void operator *= ( const double & rhs );
  void operator /= ( const double & rhs );
  void operator =  ( double rhs );
  
  /**
     A typecasting function. Obviously, care must be taken for unintended
     casting (see comments for the class, above), but this is quite a 
     convenience I think. It just calls the getValue() function.
  */
  operator double() const;
  
  /**
     Writes the state of this object to an XML stream.
     @param strm XML stream.
  */
  void toXML(ostream & strm);
  
  /*
    Formats the object
  */
  friend ostream& operator << (ostream& strm, const QMCGreensRatioComponent & rhs);

  /**
     Gets the overall value of the object.
     @return total value of this object.
  */
  double getValue() const;

  /**
     If the QMCGreensRatioComponent is bad (e.g. inf, nan, etc)
     then this will return false.
  */
  bool isNotValid() const;

  /**
     To test whether we can divide by this number.
  */
  bool isZero() const;

 private:
  
  /** 
      This class has the form k*a^b*c.
  */
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

