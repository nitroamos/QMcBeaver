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

/**************************************************************************
This SOFTWARE has been authored or contributed to by an employee or 
employees of the University of California, operator of the Los Alamos 
National Laboratory under Contract No. W-7405-ENG-36 with the U.S. 
Department of Energy.  The U.S. Government has rights to use, reproduce, 
and distribute this SOFTWARE.  Neither the Government nor the University 
makes any warranty, express or implied, or assumes any liability or 
responsibility for the use of this SOFTWARE.  If SOFTWARE is modified 
to produce derivative works, such modified SOFTWARE should be clearly 
marked, so as not to confuse it with the version available from LANL.   

Additionally, this program is free software; you can distribute it and/or 
modify it under the terms of the GNU General Public License. Accordingly, 
this program is  distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY;  without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A  PARTICULAR PURPOSE.  See the GNU General Public License 
for more details. 
**************************************************************************/



#include "QMCBasisFunction.h"

QMCBasisFunction::QMCBasisFunction()
{
}

void QMCBasisFunction::initialize(QMCFlags * FLAGS, QMCMolecule * MOL)
{
  flags = FLAGS; 
  Molecule = MOL; 
  Xcalc.allocate(3);
}

void QMCBasisFunction::initializeInterpolations()
{
  use_radial_interpolation = false;

  if( flags->use_basis_function_interpolation == 0 )
    {
      // Don't need to do anything ...
    }
  else if( flags->use_basis_function_interpolation == 1 )
    {
      // Determine the maximum number of orbitals 
      int norb = 0;
      
      for(int i=0; i<BFCoeffs.dim1(); i++)
	{
	  if( BFCoeffs(i).getNumberBasisFunctions() > norb )
	    {
	      norb = BFCoeffs(i).getNumberBasisFunctions();
	    }
	}
      
      // Allocate the splines
      RadialFunctionInterpolation.allocate(BFCoeffs.dim1(),norb);
      RadialFunctionFirstDerivativeInterpolation.allocate(BFCoeffs.dim1(),
							  norb);
      RadialFunctionSecondDerivativeInterpolation.allocate(BFCoeffs.dim1(),
							   norb);

      for(int bfc_number=0; bfc_number<BFCoeffs.dim1(); bfc_number++)
	{
	  for( int orbital=0; 
	       orbital<BFCoeffs(bfc_number).getNumberBasisFunctions(); 
	       orbital++)
	  {
	    initializeInterpolation(bfc_number,orbital,
			   RadialFunctionInterpolation(bfc_number,orbital),0);
	    initializeInterpolation(bfc_number,orbital,
	     RadialFunctionFirstDerivativeInterpolation(bfc_number,orbital),1);
	    initializeInterpolation(bfc_number,orbital,
	    RadialFunctionSecondDerivativeInterpolation(bfc_number,orbital),2);
	  }
       }

     // This must be set after the splines are initialized or else
     // it will try to evaluate the spline to initialize the spline
     use_radial_interpolation = true;
    }
  else
    {
      cerr << "ERROR: Incorrect value for use_basis_function_splines input!" 
	   << endl;
      exit(0);
    }
}

void QMCBasisFunction::initializeInterpolation(int bfc_number,int orbital, 
				CubicSplineWithGeometricProgressionGrid & S,
					       int whichDerivative)
{
  Array1D<double> x(flags->number_basis_function_interpolation_grid_points);
  Array1D<double> y(flags->number_basis_function_interpolation_grid_points);

  // Determine the grid parameter
  // The maximum distance is determined by the parameter
  const double maxDistance = 60.0;
  const double beta = pow(maxDistance/
			  flags->basis_function_interpolation_first_point,
		  2.0/flags->number_basis_function_interpolation_grid_points);

  // Determine the BF Coefficients
  QMCBasisFunctionCoefficients * BFC = &BFCoeffs(bfc_number);

  // generate x and y

  for(int i=0; i<x.dim1(); i++)
    {
      // geometric type grid
      x(i) = pow(beta,i) * (flags->basis_function_interpolation_first_point * 
			    flags->basis_function_interpolation_first_point);

      switch( whichDerivative )
	{
	case 0:
	  y(i) = radialFunction(*BFC,orbital,x(i));
	  break;
	case 1:
	  y(i) = radialFunctionFirstDerivative(*BFC,orbital,x(i));
	  break;
	case 2:
	  y(i) = radialFunctionSecondDerivative(*BFC,orbital,x(i));
	  break;
	default:
	  cerr << "ERROR: Trying to calculate unavailable derivatives in "
	       << "QMCBasisFunction::initializeInterpolation(...)" << endl;
	  exit(0);
	}
    }

  S.setGridParameters(x(0),beta);
  S.initializeWithFunctionValues(x,y,0,0);
}

void QMCBasisFunction::operator=( const QMCBasisFunction & rhs )
{
  N_BasisFunctions = rhs.N_BasisFunctions;
  flags            = rhs.flags;
  Molecule         = rhs.Molecule;
  Xcalc            = rhs.Xcalc;
  BFCoeffs         = rhs.BFCoeffs;
  BFLookupTable    = rhs.BFLookupTable;

  RadialFunctionInterpolation = rhs. RadialFunctionInterpolation;
  RadialFunctionFirstDerivativeInterpolation = 
    rhs.RadialFunctionFirstDerivativeInterpolation;
  RadialFunctionSecondDerivativeInterpolation = 
    rhs.RadialFunctionSecondDerivativeInterpolation;

  use_radial_interpolation = rhs.use_radial_interpolation;
}

int QMCBasisFunction::getNumberBasisFunctions(int i)
{
  return BFCoeffs(i).getNumberBasisFunctions();
}

istream& operator >>(istream &strm, QMCBasisFunction &rhs)
{
  rhs.BFCoeffs.allocate(rhs.flags->Natoms);

  for(int i=0; i<rhs.flags->Natoms; i++) strm >> rhs.BFCoeffs(i);

  rhs.N_BasisFunctions = 0;
  for(int i=0; i<rhs.flags->Natoms; i++)
    rhs.N_BasisFunctions += rhs.BFCoeffs(i).getNumberBasisFunctions();

  rhs.BFLookupTable.allocate(rhs.N_BasisFunctions,2);

  int ii = 0;
  int Atom = 0;
  while(ii < rhs.N_BasisFunctions)
    {
      for(int Orbital=0; Orbital < 
	    rhs.BFCoeffs(Atom).getNumberBasisFunctions(); Orbital++)
	   {
	    rhs.BFLookupTable(ii,0) = Atom;
	    rhs.BFLookupTable(ii,1) = Orbital;
       ii++;
      }
      Atom++;
    }

  rhs.initializeInterpolations();

  return strm;
}

ostream& operator <<(ostream& strm,QMCBasisFunction& rhs)
{
 strm << "&basis" << endl;
 for(int i=0; i<rhs.flags->Natoms; i++)
         strm << rhs.BFCoeffs(i);
 strm << "&" << endl;
 return strm;
}

void QMCBasisFunction::read(string runfile)
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

  if(temp_string == "&basis") input_file >> *this;
}

double QMCBasisFunction::radialFunction
       (QMCBasisFunctionCoefficients& BFC, int orbital, double r_sq)
{
  double temp = 0.0;

  int nGaussians = BFC.N_Gauss(orbital);

  double **coeffs = BFC.Coeffs.array()[orbital];

  for(int i=0; i<nGaussians; i++)
    {
      temp += coeffs[i][1]*exp(-coeffs[i][0]*r_sq);
    }

  return temp;
}

double QMCBasisFunction::radialFunctionFirstDerivative
       (QMCBasisFunctionCoefficients& BFC, int orbital, double r_sq)
{
  double temp = 0.0;

  for(int i=0; i<BFC.N_Gauss(orbital); i++)
    {
      temp += BFC.Coeffs(orbital,i,1)*BFC.Coeffs(orbital,i,0)
        *exp(-BFC.Coeffs(orbital,i,0)*r_sq);
    }

  temp *= (-2*sqrt(r_sq));

  return temp;
}

double QMCBasisFunction::radialFunctionSecondDerivative
              (QMCBasisFunctionCoefficients& BFC, int orbital, double r_sq)
{
  double temp = 0.0;

  for(int i=0; i<BFC.N_Gauss(orbital); i++)
    {
      temp += ( -2*BFC.Coeffs(orbital,i,0) + 4 * BFC.Coeffs(orbital,i,0) *
                BFC.Coeffs(orbital,i,0) * r_sq ) * BFC.Coeffs(orbital,i,1) *
        exp(-BFC.Coeffs(orbital,i,0)*r_sq);
    }

  return temp;
}

void QMCBasisFunction::evaluateBasisFunction(int which_BF, Array2D<double>& X,
	  int el_number, double& Psi, Array1D<double>& Grad, double& Laplacian)
{
  int atom_index = BFLookupTable(which_BF,0);
  int orbital_index = BFLookupTable(which_BF,1);
  
  for(int i=0; i<3; i++)
    {
      Xcalc(i) = X(el_number,i) - Molecule->Atom_Positions(atom_index,i);
    }

  double x2 = Xcalc(0)*Xcalc(0);
  double y2 = Xcalc(1)*Xcalc(1);
  double z2 = Xcalc(2)*Xcalc(2);

  double r_sq = x2+y2+z2;

  int a = BFCoeffs(atom_index).xyz_powers(orbital_index,0);
  int b = BFCoeffs(atom_index).xyz_powers(orbital_index,1);
  int c = BFCoeffs(atom_index).xyz_powers(orbital_index,2);
 
  double xyz_term = pow(Xcalc(0),a)*pow(Xcalc(1),b)*pow(Xcalc(2),c);

  double** coeffs = BFCoeffs(atom_index).Coeffs.array()[orbital_index];

  int nGaussians = BFCoeffs(atom_index).N_Gauss(orbital_index);

  Psi = 0.0;
  Grad(0) = 0.0;
  Grad(1) = 0.0;
  Grad(2) = 0.0;
  Laplacian = 0.0;

  for (int i=0; i<nGaussians; i++)
    {
      double p0 = coeffs[i][0];
      double p1 = coeffs[i][1];
      double exp_term = p1*exp(-p0*r_sq);
      double temp = -2.0*p0;
      
      Psi += exp_term;

      Grad(0) += (a/Xcalc(0)+temp*Xcalc(0))*exp_term;
      Grad(1) += (b/Xcalc(1)+temp*Xcalc(1))*exp_term;
      Grad(2) += (c/Xcalc(2)+temp*Xcalc(2))*exp_term;
      
      Laplacian += (a*(a-1.0)/x2 + b*(b-1.0)/y2 + c*(c-1.0)/z2 - \
		                    (4.0*(a+b+c)+6.0-4.0*r_sq*p0)*p0)*exp_term;
    }

  Psi *= xyz_term;
  Grad(0) *= xyz_term;
  Grad(1) *= xyz_term;
  Grad(2) *= xyz_term;
  Laplacian *= xyz_term;
}
  


