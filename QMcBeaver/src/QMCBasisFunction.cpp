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

  for(int i=0; i<BFC.N_Gauss(orbital); i++)
    {
      temp += BFC.Coeffs(orbital,i,1)
	*exp(-BFC.Coeffs(orbital,i,0)*r_sq);
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


double QMCBasisFunction::basis_function(int which_BFC, 
					QMCBasisFunctionCoefficients& BFC, 
					int orbital, Array1D <double>& X)
{
  double psi = 0;
  double r_sq = X(0)*X(0) + X(1)*X(1) + 
                X(2)*X(2);

  double xyz_term = fastPower(X(0),BFC.xyz_powers(orbital,0))
    *fastPower(X(1),BFC.xyz_powers(orbital,1))
    *fastPower(X(2),BFC.xyz_powers(orbital,2));

  if( use_radial_interpolation )
    {
      RadialFunctionInterpolation(which_BFC,orbital).evaluate(r_sq);
      psi = xyz_term * 
	RadialFunctionInterpolation(which_BFC,orbital).getFunctionValue();
    }
  else
    {
      psi = xyz_term * radialFunction(BFC, orbital, r_sq);
    }

  return psi;
}

Array1D <double> QMCBasisFunction::grad_basis_function(int which_BFC,
	  QMCBasisFunctionCoefficients& BFC, int orbital, Array1D <double>& X)
{
  int a = BFC.xyz_powers(orbital,0);
  int b = BFC.xyz_powers(orbital,1);
  int c = BFC.xyz_powers(orbital,2);

  Array1D <double> Grad(3);
  double r_sq = X(0)*X(0) + X(1)*X(1) + X(2)*X(2);

  Grad(0) = 0.0;
  Grad(1) = 0.0;
  Grad(2) = 0.0;

  double xyz_term = fastPower(X(0),a)*fastPower(X(1),b)
    *fastPower(X(2),c);

  if( use_radial_interpolation )
    {
      RadialFunctionInterpolation(which_BFC,orbital).evaluate(r_sq);
      RadialFunctionFirstDerivativeInterpolation(which_BFC,orbital).
	evaluate(r_sq);

      double radialFunctionValue = RadialFunctionInterpolation(which_BFC,
					 orbital).getFunctionValue();
      double radialFunctionFirstDerivative = 
	RadialFunctionFirstDerivativeInterpolation(which_BFC,orbital).
	getFunctionValue();

      double r = sqrt(r_sq);
      
      Grad(0) = xyz_term*(a/X(0) * radialFunctionValue + X(0)/r * 
			  radialFunctionFirstDerivative);

      Grad(1) = xyz_term*(b/X(1) * radialFunctionValue + X(1)/r *
			  radialFunctionFirstDerivative);

      Grad(2) = xyz_term*(c/X(2) * radialFunctionValue + X(2)/r * 
			  radialFunctionFirstDerivative);
    }
  else
    {
      double exp_term = 0;

      for(int i=0; i<BFC.N_Gauss(orbital); i++)
	{
	  exp_term = BFC.Coeffs(orbital,i,1)
	    *exp(-BFC.Coeffs(orbital,i,0)*r_sq);
	  
	  Grad(0) += (a/X(0)-2.0*
			     BFC.Coeffs(orbital,i,0)*X(0))
	                     *exp_term*xyz_term;
	  
	  Grad(1) += (b/X(1)-2.0*
			     BFC.Coeffs(orbital,i,0)*X(1))
	                     *exp_term*xyz_term;
	  
	  Grad(2) += (c/X(2)-2.0*
			     BFC.Coeffs(orbital,i,0)*X(2))
	                     *exp_term*xyz_term;
	}
    }

  return Grad;
}

double QMCBasisFunction::laplacian_basis_function(int which_BFC,
	 QMCBasisFunctionCoefficients& BFC, int orbital, Array1D <double>& X)
{

  double x2 = X(0)*X(0);
  double y2 = X(1)*X(1);
  double z2 = X(2)*X(2);
  double r_sq = x2 + y2 + z2;

  int a = BFC.xyz_powers(orbital,0);
  int b = BFC.xyz_powers(orbital,1);
  int c = BFC.xyz_powers(orbital,2);

  double xyz_term = fastPower(X(0),a)*fastPower(X(1),b)
                    *fastPower(X(2),c);

  double temp = 0;

  if( use_radial_interpolation )
    {
      RadialFunctionInterpolation(which_BFC,orbital).evaluate(r_sq);
      RadialFunctionFirstDerivativeInterpolation(which_BFC,orbital).
	evaluate(r_sq);
      RadialFunctionSecondDerivativeInterpolation(which_BFC,orbital).
	evaluate(r_sq);


      double radialFunctionValue = RadialFunctionInterpolation(which_BFC,
					 orbital).getFunctionValue();
      double radialFunctionFirstDerivative = 
	RadialFunctionFirstDerivativeInterpolation(which_BFC,orbital).
	getFunctionValue();
      double radialFunctionSecondDerivative = 
	RadialFunctionSecondDerivativeInterpolation(which_BFC,orbital).
	getFunctionValue();


      double r = sqrt(r_sq);

      temp = (a*(a-1)/x2 + b*(b-1)/y2 + c*(c-1)/z2) * radialFunctionValue +
	2.0/r*(a+b+c+1) * radialFunctionFirstDerivative + 
	radialFunctionSecondDerivative;
    }
  else
    {
      for(int i=0; i<BFC.N_Gauss(orbital); i++)
	{
	  temp += (a*(a-1)/x2 + b*(b-1)/y2 + c*(c-1)/z2
		   -(4*(a+b+c)+6)*BFC.Coeffs(orbital,i,0)
		   +4*(r_sq)*BFC.Coeffs(orbital,i,0)
		   *BFC.Coeffs(orbital,i,0))
	           *BFC.Coeffs(orbital,i,1)
	           *exp(-BFC.Coeffs(orbital,i,0)*r_sq);
	}
    }

  return temp*xyz_term;
}


double QMCBasisFunction::getPsi(int which_BF, Array2D <double>& X, 
				int el_number)
{
  for(int i=0; i<3; i++)
    {
      Xcalc(i) = X(el_number,i) -
	Molecule->Atom_Positions(
	  BFLookupTable(which_BF,0),i);
    }

  return basis_function(BFLookupTable(which_BF,0),
			BFCoeffs(BFLookupTable(which_BF,0)),
			BFLookupTable(which_BF,1), Xcalc);
}

Array1D <double> QMCBasisFunction::getGradPsi(int which_BF, 
					      Array2D <double>& X,
					      int el_number)
{
  for(int i=0; i<3; i++)
    {
      Xcalc(i) = X(el_number,i) -
	Molecule->Atom_Positions(
	  BFLookupTable(which_BF,0),i);
    }

  return grad_basis_function(BFLookupTable(which_BF,0),
			    BFCoeffs(BFLookupTable(which_BF,0)),
			    BFLookupTable(which_BF,1), Xcalc);
}

double QMCBasisFunction::getLaplacianPsi(int which_BF, Array2D <double>& X,
					 int el_number)
{
  for(int i=0; i<3; i++)
    {
      Xcalc(i) = X(el_number,i) -
	Molecule->Atom_Positions(
	  BFLookupTable(which_BF,0),i);
    }

  return laplacian_basis_function(BFLookupTable(which_BF,0),
			   BFCoeffs(BFLookupTable(which_BF,0)),
			   BFLookupTable(which_BF,1), Xcalc);
}








