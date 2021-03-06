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

#ifndef QMCFUNCTIONS_H
#define QMCFUNCTIONS_H

#include <iostream>
#include <sstream>

#include "QMCInput.h"
#include "QMCHartreeFock.h"
#include "QMCWalkerData.h"
#include "Stopwatch.h"

using namespace std;

class QMCFunctions
{
public:
  /**
    Creates a new instance of the class.
  */
  QMCFunctions() {};

  /**
    Creates a new instance of the class and initializes it with the data 
    controlling the QMC calculation.

    @param input input data for the calculation
  */
  QMCFunctions(QMCInput *input) : Input(input) {};
 
  /**
    Creates a new instance of the class that is identical to another 
    instance of QMCFunctions.

    @param rhs object to make a copy of
  */
  QMCFunctions(const QMCFunctions & rhs ) {};

  /**
    Deallocates all memory used by the object.
  */
  virtual ~QMCFunctions() {};

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
  virtual void evaluate(Array2D<double> &X, QMCWalkerData & data) = 0;

  virtual void evaluate(Array1D<QMCWalkerData *> &walkerData, 
		Array1D<Array2D<double> * > &xData, int num) = 0;

  /**
     This must be called after calculate_Psi_quantities.
     It assumes that we're only interested Jastrow parameter modifications,
     recalculate the energy for each set of parameters.
  */
  virtual void calculate_CorrelatedSampling(Array1D<QMCWalkerData *> &walkerData,
				    Array1D<Array2D<double> * > &xData,
				    int num){};

  virtual int getNumTimers()
    {
      return 0;
    }

  virtual void aggregateTimers(Array1D<Stopwatch> & timers,
			       int & idx){};

  /**
    Sets two QMCFunctions objects equal.

    @param rhs object to set this object equal to
  */
  void operator=(const QMCFunctions & rhs )
    {
      Input = rhs.Input;
      nalpha = rhs.nalpha;
      nbeta = rhs.nbeta;
    }

 protected:
  QMCInput *Input; 
 
  int nalpha, nbeta;
};

#endif





