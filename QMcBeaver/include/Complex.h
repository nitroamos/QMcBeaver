//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2000-2.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

#ifndef Complex_H
#define Complex_H

#include <iostream>
#include <math.h>

using namespace std;

/**
  An implementation of a complex number with the associated basic functions.
  */

class Complex
{
private:
  double re;
  double im;

public:
  /**
    Creates an object and initializes it to (0,0).
    */
  Complex();

  /**
    Creates and initializes this object.

    @param re real part of this number.
    @param im imaginary part of this number.
    */
  Complex(double re, double im);

  /**
    Creates an new instance of this object which is equal to another instance.

    @param rhs object this new object will be set equal to.
    */
  Complex(const Complex & rhs);

  /**
    Real part of this number.
    
    @return real part of this number.
    */
  double real();

  /**
    Imaginary part of this number.

    @return imaginary part of this number.
    */
  double imaginary();

  /**
    Sets two complex numbers equal.

    @rhs number to set this one equal to.
    */
  void operator=(const Complex & rhs);

  /**
    Sets a complex number and a real number equal.

    @rhs number to set this one equal to.
    */
  void operator=(const double & rhs);


  /**
    Adds two complex numbers.
    
    @return sum of the arguments.
    */
  Complex operator+(const Complex & rhs);

  /**
    Adds a complex and a real number.
    
    @return sum of the arguments.
    */
  Complex operator+(const double & rhs);


  /**
    Subtracts two complex numbers.

    @return difference of the arguments.
    */
  Complex operator-(const Complex & rhs);

  /**
    Subtracts a complex and a real number.

    @return difference of the arguments.
    */
  Complex operator-(const double & rhs);

  /**
    Multiplies two complex number.

    @return product of the arguments.
    */
  Complex operator*(const Complex & rhs);

  /**
    Multiplies a complex and a real number.

    @return product of the arguments.
    */
  Complex operator*(const double & rhs);

  /**
    Divides two complex numbers.

    @return result of the division.
    */
  Complex operator/(const Complex & rhs);

  /**
    Calculates the complex conjugate of this number.

    @return complex conjugate of this number.
    */
  Complex conjugate();

  /**
    Calculates the magniutde of this complex number.
    \f[
    c.abs() = \sqrt{\left( c.re()^{2} + c.im()^{2} \right)}
    \f]

    @return magnitude of this complex number.
    */
  double abs();

  /**
    Calculates the square root of this complex number.

    @return square root of this complex number.
    */
  Complex squareroot();

  /**
    Write the number to an output stream.
    */
  friend ostream & operator<<(ostream & strm, Complex & c);
};

#endif
