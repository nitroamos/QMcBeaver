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

#include "QMCBasisFunctionCoefficients.h"
#include "StringManipulation.h"
#include <iomanip>

QMCBasisFunctionCoefficients::QMCBasisFunctionCoefficients() 
{
  N_Orbitals = 0;
}

int QMCBasisFunctionCoefficients::getNumberBasisFunctions()
{
  return N_Orbitals;
}

void QMCBasisFunctionCoefficients::operator=(
				   const QMCBasisFunctionCoefficients & rhs)
{
  this->Label         = rhs.Label;
  this->N_Orbitals    = rhs.N_Orbitals;
  this->Max_Gaussians = rhs.Max_Gaussians;
  this->Coeffs        = rhs.Coeffs;
  this->xyz_powers    = rhs.xyz_powers;
  this->N_Gauss       = rhs.N_Gauss;
  this->Type          = rhs.Type;
}

void type_to_xyz(string type, int &x, int &y, int &z)
{
  type = StringManipulation::toAllLower(type);
  const char * bf = type.c_str();
  const int maxl = 7;
  char lnames[maxl] = {'s','p','d','f','g','h','i'};

  int idx = -1;
  for(int i=0; i<maxl; i++)
    if(bf[0] == lnames[i]){
      idx = i;
      break;
    }
  
  if(idx == -1){
    clog << "ERROR: Unknown Basis Function (" << type
	 << ")" << endl;
    clog << "You requested type = " << bf[0] << " but we  only have ";
    for(int i=0; i<maxl; i++)
      clog << lnames[i];
    clog << " available." << endl;
    exit(1);
  }

  if((int)type.size()-1 != idx){
    clog << "ERROR: Basis Function (" << type << ")"
	 << " doesn't have angular momentum = " << idx << "." << endl;
    exit(1);
  }

  x = y = z = 0;
  for(unsigned int i=1; i<type.size(); i++)
    {
      switch (bf[i]){
      case 'x': x++; break;
      case 'y': y++; break;
      case 'z': z++; break;
      default:
	clog << "Error: unknown angular component " << bf[i]
	     << " in basis function " << type << endl;
	exit(0);
      }
    }
}

istream& operator >>(istream &strm, QMCBasisFunctionCoefficients &rhs)
{
  strm >> rhs.Label;
  rhs.Label = StringManipulation::toFirstUpperRestLower( rhs.Label );

  strm >> rhs.N_Orbitals;
  strm >> rhs.Max_Gaussians;

  rhs.Coeffs.allocate(rhs.N_Orbitals,rhs.Max_Gaussians,2);
  rhs.xyz_powers.allocate(rhs.N_Orbitals,3);
  rhs.N_Gauss.allocate(rhs.N_Orbitals);
  rhs.Type.allocate(rhs.N_Orbitals);

  int x,y,z;
  rhs.lmax = 0;
  for(int i=0; i<rhs.N_Orbitals; i++)
    {
      strm >> rhs.N_Gauss(i);
      strm >> rhs.Type(i);
      
      type_to_xyz(rhs.Type(i),x,y,z);
      rhs.xyz_powers(i,0) = x;
      rhs.xyz_powers(i,1) = y;
      rhs.xyz_powers(i,2) = z;
      rhs.lmax = max(rhs.lmax,x+y+z);
      
      for(int j=0; j<rhs.N_Gauss(i); j++)
	{
	  strm >> rhs.Coeffs(i,j,0);
	  strm >> rhs.Coeffs(i,j,1);
	}
    }
  return strm;
}

void QMCBasisFunctionCoefficients::read(string runfile)
{
  ifstream input_file(runfile.c_str());

  if(!input_file)
    {
      cerr << "ERROR: Could not open " << runfile << "!" << endl;
      exit(1);
    }

  string temp_string;
  input_file >> temp_string;
  while((temp_string != "&basis") && (input_file.eof() != 1))
    {
      input_file >> temp_string;
    }

  input_file >> *this;
}

ostream& operator <<(ostream& strm, QMCBasisFunctionCoefficients& rhs)
{
  int width = 20;
  strm << setw(3) << left << rhs.Label
       << setw(5) << right << rhs.N_Orbitals
       << setw(5) << rhs.Max_Gaussians << endl;

  strm.unsetf(ios::fixed);
  strm.unsetf(ios::scientific);

  for(int i=0; i<rhs.N_Orbitals; i++)
    {
      strm << right;
      strm << setw(3) << rhs.N_Gauss(i) << " "
	   << setw(5) << left << rhs.Type(i) << endl;
      for(int j=0; j<rhs.N_Gauss(i); j++)
	{
	  strm << setw(width) << right << rhs.Coeffs(i,j,0) << "  ";

	  if(rhs.Coeffs(i,j,1) < 0.0) strm << " -";
	  else                        strm << "  ";
	  strm << setw(width) << left << fabs(rhs.Coeffs(i,j,1)) << endl;
	}
    }
  strm << right;
  return strm;
}

