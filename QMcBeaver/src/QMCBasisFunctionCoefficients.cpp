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
  if(type == "s" )         {x=0; y=0; z=0;}
  else if(type == "px" )   {x=1; y=0; z=0;}
  else if(type == "py" )   {x=0; y=1; z=0;}
  else if(type == "pz" )   {x=0; y=0; z=1;}
  else if(type == "dxx" )  {x=2; y=0; z=0;}
  else if(type == "dxy" )  {x=1; y=1; z=0;}
  else if(type == "dxz" )  {x=1; y=0; z=1;}
  else if(type == "dyy" )  {x=0; y=2; z=0;}
  else if(type == "dyz" )  {x=0; y=1; z=1;}
  else if(type == "dzz" )  {x=0; y=0; z=2;}
  else if(type == "fxxx" ) {x=3; y=0; z=0;}
  else if(type == "fyyy" ) {x=0; y=3; z=0;}
  else if(type == "fzzz" ) {x=0; y=0; z=3;}
  else if(type == "fyyx" ) {x=1; y=2; z=0;}
  else if(type == "fxxy" ) {x=2; y=1; z=0;}
  else if(type == "fxxz" ) {x=2; y=0; z=1;}
  else if(type == "fzzx" ) {x=1; y=0; z=2;}
  else if(type == "fzzy" ) {x=0; y=1; z=2;}
  else if(type == "fyyz" ) {x=0; y=2; z=1;}
  else if(type == "fxyz" ) {x=1; y=1; z=1;}

  else
    {
      cerr << "ERROR: Unknown Basis Function Type (" << type
           << ")" << endl;
      exit(1);
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

