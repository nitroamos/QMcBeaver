//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2002.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

#ifndef CubicSplineWithGeometricProgressionGrid_H
#define CubicSplineWithGeometricProgressionGrid_H

#include "CubicSpline.h"

/**
   A 1-dimensional (\f$\mathbf{R}^{1} \rightarrow \mathbf{R}^{1}\f$) cubic
   spline interpolation with a grid that is assumed to be spaced according
   to a geometric relationship for faster evaluation.  
   \f[
   x_{i+1} = \beta x_{i}
   \f]
   \f$x_{0}\f$ is set equal to the first datum used to initialize the spline.
*/
  

class CubicSplineWithGeometricProgressionGrid : public CubicSpline
{
public:
  /**
    Constructs an uninitialized spline.
    */
  CubicSplineWithGeometricProgressionGrid();

  /**
    Sets the value for \beta and \f$x_{0}\f$ used in generating this grid.   
    \f[
    x_{i+1} = e x_{i}
    \f]

    @param x0 the first point in the grid \f$x_{i+1} = \beta x_{i}\f$.
    */
  void setGridParameters(double x0, double beta);

  void evaluate(double x);

  /**
     Sets two CubicSplineWithGeometricProgressionGrid objects equal.

    @param rhs object to set this object equal to
  */

  void operator=( const CubicSplineWithGeometricProgressionGrid & rhs );


private:
  double x0;

  /**
     \frac{1}{ln(\beta)}
  */
  double inverselnbeta;
};

#endif
