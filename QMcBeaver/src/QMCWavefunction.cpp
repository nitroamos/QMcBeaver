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
    Ndeterminants = 0;
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

int QMCWavefunction::getNumberDeterminants()
{
    return Ndeterminants;
}

QMCWavefunction QMCWavefunction::operator=( const QMCWavefunction & rhs )
{
    Norbitals       = rhs.Norbitals;
    Nbasisfunc      = rhs.Nbasisfunc;
    Nalpha          = rhs.Nalpha;
    Nbeta           = rhs.Nbeta;
    Nelectrons      = rhs.Nelectrons;
    Ndeterminants   = rhs.Ndeterminants;
    Coeffs          = rhs.Coeffs;
    CI_coeffs       = rhs.CI_coeffs;
    AlphaOccupation = rhs.AlphaOccupation;
    BetaOccupation  = rhs.BetaOccupation;
    return *this;
}

/** The QSC didn't like my little commenting trick, so I just got rid
    of the old streaming operator -- the old files are still in CVS of course...

    This operator >> will read input files in the multi configuration style.
*/
istream& operator >>(istream &strm,QMCWavefunction &rhs)
{
  rhs.Coeffs.allocate(rhs.Norbitals,rhs.Nbasisfunc);
  rhs.AlphaOccupation.allocate(rhs.Ndeterminants,rhs.Norbitals);
  rhs.BetaOccupation.allocate(rhs.Ndeterminants,rhs.Norbitals);
  rhs.CI_coeffs.allocate(rhs.Ndeterminants);
    
  string temp_string;
  for(int i=0; i<rhs.Norbitals; i++)
    for(int j=0; j<rhs.Nbasisfunc;j++)
      {
	strm >> rhs.Coeffs(i,j);
      }
            
  strm >> temp_string;
  strm >> temp_string;
    
  for (int i=0; i<rhs.Ndeterminants; i++)
    for (int j=0; j<rhs.Norbitals; j++)
      {
	strm >> rhs.AlphaOccupation(i,j);
      }
            
  strm >> temp_string;
  strm >> temp_string;
    
  for (int i=0; i<rhs.Ndeterminants; i++)
    for (int j=0; j<rhs.Norbitals; j++)
      {
	strm >> rhs.BetaOccupation(i,j);
      }  
            
  strm >> temp_string;
  strm >> temp_string;
    
  for (int i=0; i<rhs.Ndeterminants; i++)
    {
      strm >> rhs.CI_coeffs(i);
    }
    
  rhs.Nalpha = 0;
  rhs.Nbeta = 0;
  for(int i=0; i<rhs.Norbitals; i++)
    {
      rhs.Nalpha += rhs.AlphaOccupation(0,i);
      rhs.Nbeta += rhs.BetaOccupation(0,i);
    } 
    
  rhs.Nelectrons = rhs.Nalpha + rhs.Nbeta;
  return strm;
}

void QMCWavefunction::read(int numberOrbitals, int numberBasisFunctions,
                           int numberDeterminants, string runfile)
{
    Norbitals  = numberOrbitals;
    Nbasisfunc = numberBasisFunctions;
    Ndeterminants = numberDeterminants;
    
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
        for(int j=0; j<rhs.Nbasisfunc; j++)
        {
            if( j%3 == 0 && j > 0 )
            {
                strm << endl;
            }
            strm << rhs.Coeffs(i,j) << "\t";
        }
        strm << endl << endl;
    }
    
    strm << "Alpha Occupation" << endl;
    for (int i=0; i<rhs.Ndeterminants; i++)
    {
        for (int j=0; j<rhs.Norbitals; j++)
        {
            strm << rhs.AlphaOccupation(i,j) << "\t";
        }
        strm << endl;
    }
    
    strm << endl << "Beta Occupation" << endl;
    for (int i=0; i<rhs.Ndeterminants; i++)
    {
        for (int j=0; j<rhs.Norbitals; j++)
        {
            strm << rhs.BetaOccupation(i,j) << "\t";
        }
        strm << endl;
    }
    
    strm << endl << "CI Coeffs" << endl;
    for (int i=0; i<rhs.Ndeterminants; i++)
    {
        strm << rhs.CI_coeffs(i) << endl;
    }
    strm << endl << endl;
    
    strm << "&" << endl;
    return strm;
}



