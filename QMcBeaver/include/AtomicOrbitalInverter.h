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

#ifndef ATOMIC_ORBITAL_INVERTER_H
#define ATOMIC_ORBITAL_INVERTER_H

#define TINY (1e-12)
#define PI 3.14159265359

#include <iostream>
#include <math.h>

#include "Array1D.h"
#include "LinearSpline.h"
#include "random.h"
#include "DistributionInverter.h"
//#include "QMCBasisFunctionCoefficients.h" //NOT DONE

using namespace std;
/**
   This class can invert any atomic orbital.
*/

class AtomicOrbitalInverter
{
public:

  /**
     Creates a new uninitialized instance of this class.
  */
  AtomicOrbitalInverter();
  
  /**
     Copy the rhs oject into this current object.

     @param rhs object to set this equal to.
  */
  void operator=( const AtomicOrbitalInverter & rhs );

  /**
     Set the number of grid points in each cooridnant.
  */
  void set_npoints(int points_r,int points_theta,int points_phi);

  /**
     Initialize the total atomic orbital.
     coeff array c and exponents b in: 
     Sum(ci*exp(-bi*x^2)exp(-bi*y^2)exp(-bi*z^2),i)
     @param Xn is the power on X
     @param Yn is the power on Y
     @param Zn is the power on Z
     @param B is the array of the exponents
     @param C is the array of the prefactors 
  */
  void initialize(int Xn, int Yn, int Zn,Array1D<double> B,Array1D<double> C);


  //NOT DONE
  //  /**
  // Initialize the total atomic orbital.
  // @param BSC is the basis set coefficients
  //*/
  //void initialize(QMCBasisFunctionCoefficients &BSC);
  //END NOT DONE


  /**
     given a uniform random number this will return the xyz coords
  */
  void get_xyz(long & iseed, double &x, double &y, double &z);

  /**
     Use the current b and c arrays to evaluate the gaussian part
  */
  double eval_gaussians(double r);

  /**
     Invert the current gaussian part
  */
  double invert_gaussians(long & iseed);

private:

  Array1D<double> b;
  Array1D<double> c;
  Array1D<double> x_inp;
  Array1D<double> y_inp;
  int xn,yn,zn;          //the polarization powers in front of the AO
  double cutoff1;         //tells us how far out to do the linear spline
  double cutoff2;         //tells us how far out to do the linear spline
  double cutoff3;         //tells us how far out to do the linear spline
  int npoints_r;       //tells us how many points do do in each linear spline
  int npoints_theta;   //tells us how many points do do in each linear spline
  int npoints_phi;     //tells us how many points do do in each linear spline
  DistributionInverter r_form;      //r part of atomic orbital
  DistributionInverter theta_form;  //theta part of atomic orbital
  DistributionInverter phi_form;    //phi part of atomic orbital

};


#endif





