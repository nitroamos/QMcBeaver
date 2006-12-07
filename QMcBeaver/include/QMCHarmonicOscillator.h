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

#ifndef QMCHarmonicOscillator_H
#define QMCHarmonicOscillator_H

#include <iostream>
#include <sstream>

#include "QMCInput.h"
#include "QMCGreensRatioComponent.h"
#include "QMCConfigIO.h"
#include "IeeeMath.h"
#include "QMCWalkerData.h"
#include "QMCFunctions.h"

using namespace std;

class QMCHarmonicOscillator : public QMCFunctions
{
public:
  /**
    Creates a new instance of the class.
  */
  QMCHarmonicOscillator() {};

  /**
    Creates a new instance of the class and initializes it with the data 
    controlling the QMC calculation.

    @param input input data for the calculation
  */
  QMCHarmonicOscillator(QMCInput *input);

  /**
    Creates a new instance of the class that is identical to another 
    instance of QMCHarmonicOscillator.

    @param rhs object to make a copy of
  */
  QMCHarmonicOscillator(const QMCHarmonicOscillator & rhs )
    {
      *this = rhs;
    }

  /**
    Deallocates all memory used by the object.
  */
  ~QMCHarmonicOscillator() {};

  /**
    Evaluates all of the calculated properties at X and places the calculated
    data into the QMCWalkerData struct provided. Two overloaded functions are
    provided, one of them processes a array of parameters, the other processes
    just one (useful during a QMCWalker's initialization)

    @param X \f$3N\f$ dimensional configuration of electrons represented by 
    a \f$N \times 3\f$ matrix
    @param data all the data that a QMCWalker should require
    @param writeConfig if the program is writing configs, we need to know here.
    if true, the walkerData.configOutput will be given its info
  */

  void evaluate(Array1D<QMCWalkerData *> &walkerData, Array1D<Array2D<double> * > &xData, int num, bool writeConfig);

  void evaluate(Array2D<double> &X, QMCWalkerData & data)
    {
      Array1D< Array2D<double>* > tempX;
      tempX.allocate(1);
      tempX(0) = &X;
      Array1D< QMCWalkerData * > tempWD;
      tempWD.allocate(1);
      tempWD(0) = &data;
      evaluate(tempWD,tempX,1,false);
    }

  /**
    Sets two QMCHarmonicOscillator objects equal.

    @param rhs object to set this object equal to
  */
  void operator=(const QMCHarmonicOscillator & rhs )
    {
      Input = rhs.Input;
    }

 private:
  double a, w, w2div2, a2;

};

#endif





