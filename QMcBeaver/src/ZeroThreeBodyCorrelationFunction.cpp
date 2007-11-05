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

#include "ZeroThreeBodyCorrelationFunction.h"

void ZeroThreeBodyCorrelationFunction::initializeParameters(int en, int ee,
                            Array1D<double> &Parameters, int power, double max)
{
  grad1.allocate(3);
  grad1 = 0.0;
  
  grad2.allocate(3);
  grad2 = 0.0;
}

bool ZeroThreeBodyCorrelationFunction::setElectron(bool first, Array1D<double> &xyz, double dist)
{
  return false;
}

void ZeroThreeBodyCorrelationFunction::evaluate(Array1D<double> &xyz12, double r12)
{
  // Do nothing
}

double ZeroThreeBodyCorrelationFunction::getFunctionValue()
{
  return 0.0;
}

double ZeroThreeBodyCorrelationFunction::get_p_a(int ai)
{
  return 0.0;
}

Array1D<double>* ZeroThreeBodyCorrelationFunction::getElectron1Gradient()
{
  return &grad1;
}

Array1D<double>* ZeroThreeBodyCorrelationFunction::getElectron2Gradient()
{
  return &grad2;
}

double ZeroThreeBodyCorrelationFunction::get_p2_xa(bool e1, int xyz, int ai)
{
  return 0.0;
}

double ZeroThreeBodyCorrelationFunction::getLaplacianValue()
{
  return 0.0;
}

double ZeroThreeBodyCorrelationFunction::get_p3_xxa(int ai)
{
  return 0.0;
}

double ZeroThreeBodyCorrelationFunction::getCutoffDist()
{
  return 0.0;
}
