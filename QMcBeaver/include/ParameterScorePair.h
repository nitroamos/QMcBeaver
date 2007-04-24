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

#ifndef PARAMETERSCOREPAIR_H
#define PARAMETERSCOREPAIR_H

#include "Array1D.h"
#include "QMCObjectiveFunctionResult.h"
#include <iomanip>

/**
  A container which holds a set of parameters and an associated scalar
  score value.
  */

class ParameterScorePair
{
public:
  /**
    Creates an uninitialized instance of this class with no allocated memory.
    */
  ParameterScorePair();

  /**
    Creates an uninitialized instance of this class and sets the score and
    parameter values.

    @param score Score.
    @param parameters Parameters.
    */
  ParameterScorePair(QMCObjectiveFunctionResult score, Array1D<double> & parameters);

  /**
    Creates an instance of this class which is equal to another instance.
    */
  ParameterScorePair(const ParameterScorePair &PSP);

  /**
    Gets the score.

    @return score.
    */
  double getScore() const;

  /**
    Gets the parameters.

    @return paramters.
    */
  Array1D<double> * getParameters();

  /**
    Set two ParameterScorePair objects equal.

    @param rhs object to set this object equal to.
    */
  void operator=(const ParameterScorePair & rhs);

  /**
    An operator which orders ParameterScorePair objects based on their scores.
    */
  bool operator<(const ParameterScorePair &PSP) const;

  /**
     Prints the contents of this object in a human readable format.
  */
  friend ostream& operator<<(ostream & strm, 
			     const ParameterScorePair & rhs);

private:
  QMCObjectiveFunctionResult Score;
  Array1D<double> Parameters;
};

#endif
