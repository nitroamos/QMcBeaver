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

#include "ParameterScorePair.h"

ParameterScorePair::ParameterScorePair()
{
}

ParameterScorePair::ParameterScorePair(double score, 
				       Array1D<double> & parameters)
{
  this->Score = score;
  this->Parameters = parameters;
}

ParameterScorePair::ParameterScorePair(const ParameterScorePair & PSP)
{
  *this = PSP;
}

void ParameterScorePair::operator=(const ParameterScorePair &PSP)
{
  Score      = PSP.Score;
  Parameters = PSP.Parameters;
}

bool ParameterScorePair::operator<(ParameterScorePair &PSP)
{
  bool returnvalue = true;
  if( Score < PSP.Score )
    {
      returnvalue = true;
    }
  else
    {
      returnvalue = false;
    }

  return returnvalue;
}

double ParameterScorePair::getScore()
{
  return Score;
}

Array1D<double> * ParameterScorePair::getParameters()
{
  return &Parameters;
}
