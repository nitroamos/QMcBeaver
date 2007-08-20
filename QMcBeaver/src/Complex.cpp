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

#include "Complex.h"

Complex::Complex()
{
  re = 0.0;
  im = 0.0;
}

Complex::Complex(double re, double im)
{
  this->re = re;
  this->im = im;
}

Complex::Complex(const Complex & rhs)
{
  *this = rhs;
}

double Complex::real()
{
  return re;
}

double Complex::imaginary()
{
  return im;
}

void Complex::operator=(const Complex & rhs)
{
  this->re = rhs.re;
  this->im = rhs.im;
}

void Complex::operator=(const double & rhs)
{
  this->re = rhs;
  this->im = 0;
}

Complex Complex::operator+(const Complex & rhs)
{
  Complex result;

  result.re = this->re + rhs.re;
  result.im = this->im + rhs.im;

  return result;
}

Complex Complex::operator+(const double & rhs)
{
  Complex result;

  result.re = this->re + rhs;
  result.im = this->im;

  return result;
}

Complex Complex::operator-(const Complex & rhs)
{
  Complex result;

  result.re = this->re - rhs.re;
  result.im = this->im - rhs.im;

  return result;
}

Complex Complex::operator-(const double & rhs)
{
  Complex result;

  result.re = this->re - rhs;
  result.im = this->im;

  return result;
}

Complex Complex::operator*(const Complex & rhs)
{
  Complex result;

  result.re = this->re*rhs.re - this->im*rhs.im;
  result.im = this->im*rhs.re + this->re*rhs.im;

  return result;
}

Complex Complex::operator*(const double & rhs)
{
  Complex result;

  result.re = this->re*rhs;
  result.im = this->im*rhs;

  return result;
}

Complex Complex::operator/(const Complex & rhs)
{
  Complex result;

  if( fabs(rhs.re) >= fabs(rhs.im) )
    {
      double r = rhs.im/rhs.re;
      double den = rhs.re+r*rhs.im;
      result.re = (this->re+r*this->im)/den;
      result.im = (this->im-r*this->re)/den;
    }
  else
    {
      double r = rhs.re/rhs.im;
      double den = rhs.im+r*rhs.re;
      result.re = (this->re*r+this->im)/den;
      result.im = (this->im*r-this->re)/den;
    }

  return result;
}

Complex Complex::conjugate()
{
  Complex result;

  result.re = this->re;
  result.im = -this->im;

  return result;
}

double Complex::abs()
{
  double x = fabs(re);
  double y = fabs(im);
  double ans;

  if(x == 0.0) 
    {
      ans = y;
    }
  else if(y == 0.0) 
    {
      ans = x;
    }
  else if( x > y)
    {
      double temp = y/x;
      ans = x*sqrt(1.0+temp*temp);
    }
  else
    {
      double temp = x/y;
      ans = y*sqrt(1.0+temp*temp);
    }

  return ans;
}

Complex Complex::squareroot()
{
  Complex c;

  if ((this->re == 0.0) && (this->im == 0.0)) 
    {
      c.re = 0.0;
      c.im = 0.0;
      return c;
    } 
  else 
    {
      double x = fabs(this->re);
      double y = fabs(this->im);
      double w,r;

      if (x >= y) 
	{
	  r = y/x;
	  w = sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
	} 
      else 
	{
	  r = x/y;
	  w = sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
	}
      
      if (this->re >= 0.0) 
	{
	  c.re = w;
	  c.im = this->im/(2.0*w);
	} 
      else 
	{
	  c.im = (this->im >= 0) ? w : -w;
	  c.re = this->im/(2.0*c.im);
	}
      return c;
    }
}

ostream & operator<<(ostream & strm, Complex & c)
{
  if(c.real() >= 0.0) strm << " ";
  strm << c.real();
  strm << (c.imaginary() > 0? "+i":"-i") << fabs(c.imaginary());

  return strm;
}









