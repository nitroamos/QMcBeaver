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

#include "CubicSplineWithGeometricProgressionGrid.h"

CubicSplineWithGeometricProgressionGrid::
CubicSplineWithGeometricProgressionGrid():CubicSpline()
{
  x0   = 0;
  inverselnbeta = 0;
}

void CubicSplineWithGeometricProgressionGrid::setGridParameters(double x0,
								double beta)
{
  this->x0   = x0;
  this->inverselnbeta = 1.0/log(beta);
}

void CubicSplineWithGeometricProgressionGrid::evaluate(double x)
{
  // which box of the domain is x in?
  int box = (int) ( log(x/x0) * inverselnbeta );

  CubicSpline::evaluate(x,box);
}

void CubicSplineWithGeometricProgressionGrid::
operator=( const CubicSplineWithGeometricProgressionGrid & rhs )
{
  CubicSpline::operator=(rhs);
  this->x0   = rhs.x0;
  this->inverselnbeta = rhs.inverselnbeta;
}
