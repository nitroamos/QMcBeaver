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

ParameterScorePair::ParameterScorePair(QMCObjectiveFunctionResult score, 
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

bool ParameterScorePair::operator<(const ParameterScorePair &PSP) const
{
  bool returnvalue = true;
  if( getScore() < PSP.getScore() )
    {
      returnvalue = true;
    }
  else
    {
      returnvalue = false;
    }

  return returnvalue;
}

ostream& operator<<(ostream & strm, const ParameterScorePair & rhs)
{
  strm << "Score " << setw(20) << rhs.getScore() << endl;
  strm << "     ";
  strm << "Parameters    " << rhs.Parameters;
  strm << "     ";
  strm << "Energy      " << setw(20) << rhs.Score.getEnergyAve() << " +/- " << rhs.Score.getEnergyVar();
  strm << " (# samples = " << rhs.Score.getNumberSamples() << ")" << endl;
  return strm;
}

double ParameterScorePair::getScore() const
{
  double score = Score.getScore();
  return score;
}

Array1D<double> * ParameterScorePair::getParameters()
{
  return &Parameters;
}
