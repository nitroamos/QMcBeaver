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

#include "QMCWavefunction.h"

QMCWavefunction::QMCWavefunction()
{
  Norbitals = 0; 
  Nbasisfunc = 0; 
  Nalpha = 0; 
  Nbeta = 0; 
  Nelectrons = 0;
}

int QMCWavefunction::getNumberOrbitals()
{
  return Norbitals;
}

int QMCWavefunction::getNumberBasisFunctions()
{
  return Nbasisfunc;
}

int QMCWavefunction::getNumberAlphaElectrons()
{
  return Nalpha;
}

int QMCWavefunction::getNumberBetaElectrons()
{
  return Nbeta;
}

int QMCWavefunction::getNumberElectrons()
{
  return Nelectrons;
}
  

QMCWavefunction QMCWavefunction::operator=( const QMCWavefunction & rhs )
{
  Norbitals       = rhs.Norbitals;
  Nbasisfunc      = rhs.Nbasisfunc;
  Nalpha          = rhs.Nalpha;
  Nbeta           = rhs.Nbeta;
  Nelectrons      = rhs.Nelectrons;
  Coeffs          = rhs.Coeffs;
  AlphaOccupation = rhs.AlphaOccupation;
  BetaOccupation  = rhs.BetaOccupation;
  return *this;
}

istream& operator >>(istream &strm,QMCWavefunction &rhs)
{
  rhs.Coeffs.allocate(rhs.Nbasisfunc,rhs.Norbitals);
  rhs.AlphaOccupation.allocate(rhs.Norbitals);
  rhs.BetaOccupation.allocate(rhs.Norbitals);

  string temp_string;
  for(int i=0; i<rhs.Norbitals; i++)
    {
      strm >> rhs.AlphaOccupation(i) 
	   >> rhs.BetaOccupation(i);
      
      for(int j=0; j<rhs.Nbasisfunc;j++)
	{
	  strm >> rhs.Coeffs(j,i);
	}
    }

  rhs.Nalpha = 0;
  rhs.Nbeta = 0;
  for(int i=0; i<rhs.Norbitals; i++)
    {
      rhs.Nalpha += rhs.AlphaOccupation(i);
      rhs.Nbeta += rhs.BetaOccupation(i);
    } 

  rhs.Nelectrons = rhs.Nalpha + rhs.Nbeta;
  return strm;
}

void QMCWavefunction::read(int numberOrbitals, int numberBasisFunctions,
			   string runfile)
{
  Norbitals  = numberOrbitals;
  Nbasisfunc = numberBasisFunctions;

  ifstream input_file(runfile.c_str());

  if(!input_file)
    {
      cerr << "ERROR: Could not open " << runfile << "!" << endl;
      exit(1);
    }

  string temp_string;
  input_file >> temp_string;
  while((temp_string != "&wavefunction") && (input_file.eof() != 1))
    {
      input_file >> temp_string;
    }

  input_file >> *this;
}

ostream& operator <<(ostream& strm, QMCWavefunction& rhs)
{
 strm << "&wavefunction" << endl;
 for(int i=0; i<rhs.Norbitals; i++)
 {
  strm << rhs.AlphaOccupation(i) << "\t"
       << rhs.BetaOccupation(i) << endl;
  for(int j=0; j<rhs.Nbasisfunc; j++)
    {
      if( j%3 == 0 && j > 0 )
        {
          strm << endl;
        }

      strm << rhs.Coeffs(j,i) << "\t";
    }
  strm << endl;
 }
 strm << "&" << endl;
 return strm;
}



