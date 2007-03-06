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
{}

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

void QMCBasisFunction::
initializeInterpolation(int bfc_number,int orbital,
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

QMCBasisFunctionCoefficients* QMCBasisFunction::getBFCoeffs(int i)
{
  return &BFCoeffs(i);
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

double QMCBasisFunction::
radialFunction(QMCBasisFunctionCoefficients& BFC,
               int orbital, double r_sq)
{
  double temp = 0.0;

  int nGaussians = BFC.N_Gauss(orbital);

  qmcfloat **coeffs = BFC.Coeffs.array()[orbital];

  for(int i=0; i<nGaussians; i++)
    {
      temp += coeffs[i][1]*exp(-coeffs[i][0]*r_sq);
    }

  return temp;
}

double QMCBasisFunction::
radialFunctionFirstDerivative(QMCBasisFunctionCoefficients& BFC,
                              int orbital, double r_sq)
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

void QMCBasisFunction::
evaluateBasisFunctions(Array2D<double>& X, int start, int stop,
                       Array2D<qmcfloat>& chi_value,
                       Array2D<qmcfloat>& chi_grx,
                       Array2D<qmcfloat>& chi_gry,
                       Array2D<qmcfloat>& chi_grz,
                       Array2D<qmcfloat>& chi_laplacian)
{
#if defined SINGLEPRECISION || defined QMC_GPU
  const float TOOSMALL = 1e-35;
#else
  const double TOOSMALL = 1e-306;
#endif

  int el = 0, bf;
  int a, b, c, nGaussians;
  int numBF;

  qmcfloat x, y, z, x2, y2, z2, r_sq, xyz;
  qmcfloat p0, p1, exp_term, temp;
  qmcfloat chi, chi_gradx, chi_grady, chi_gradz, chi_lap;

  for (int index=start; index<=stop; index++)
    {
      bf = 0;
      for (int atom=0; atom<flags->Natoms; atom++)
        {
          for (int i=0; i<3; i++)
            {
              Xcalc(i) = X(index,i) - Molecule->Atom_Positions(atom,i);
            }
          x = Xcalc(0);
          y = Xcalc(1);
          z = Xcalc(2);
          x2 = x*x;
          y2 = y*y;
          z2 = z*z;
          r_sq = x2+y2+z2;
          numBF = BFCoeffs(atom).getNumberBasisFunctions();
          for (int j=0; j<numBF; j++)
            {
              a = BFCoeffs(atom).xyz_powers(j,0);
              b = BFCoeffs(atom).xyz_powers(j,1);
              c = BFCoeffs(atom).xyz_powers(j,2);

              xyz = pow(x,a)*pow(y,b)*pow(z,c);

              qmcfloat** coeffs = BFCoeffs(atom).Coeffs.array()[j];

              nGaussians = BFCoeffs(atom).N_Gauss(j);
              chi       = 0;
              chi_gradx = 0;
              chi_grady = 0;
              chi_gradz = 0;
              chi_lap   = 0;
              for (int i=0; i<nGaussians; i++)
                {
                  p0 = coeffs[i][0];
                  p1 = coeffs[i][1];
                  exp_term = p1*exp(-p0*r_sq)*xyz;
                  temp = -2.0*p0;

                  chi       += exp_term;
                  chi_gradx += (a/x + temp*x) * exp_term;
                  chi_grady += (b/y + temp*y) * exp_term;
                  chi_gradz += (c/z + temp*z) * exp_term;
                  chi_lap   += (a*(a-1.0)/x2 + b*(b-1.0)/y2 + c*(c-1.0)/z2
                                - (4.0*(a+b+c)+6.0-4.0*r_sq*p0)*p0)*exp_term;
                }

              /*
		This next block is unneeded if we're flushing underflow to zero
		or if the computer doesn't slow down too much with denormal
		values.

		Allowing denormal floating point values can really slow down the machine since
		they are likely handled in software (as opposed to hardware).

		Further, I was finding that some machines were crashing when transcendental
		math functions were evaluated with denormals later in the code...
	      */	      
	      /*
              if(chi > 0 && chi < TOOSMALL) chi = TOOSMALL;
              if(chi < 0 && chi > -1.0*TOOSMALL) chi = -1.0*TOOSMALL;
              if(chi_gradx > 0 && chi_gradx < TOOSMALL) chi_gradx = TOOSMALL;
              if(chi_gradx < 0 && chi_gradx > -1.0*TOOSMALL) chi_gradx = -1.0*TOOSMALL;
              if(chi_grady > 0 && chi_grady < TOOSMALL) chi_grady = TOOSMALL;
              if(chi_grady < 0 && chi_grady > -1.0*TOOSMALL) chi_grady = -1.0*TOOSMALL;
              if(chi_gradz > 0 && chi_gradz < TOOSMALL) chi_gradz = TOOSMALL;
              if(chi_gradz < 0 && chi_gradz > -1.0*TOOSMALL) chi_gradz = -1.0*TOOSMALL;
              if(chi_lap > 0 && chi_lap < TOOSMALL) chi_lap = TOOSMALL;
              if(chi_lap < 0 && chi_lap > -1.0*TOOSMALL) chi_lap = -1.0*TOOSMALL;
	      */

              chi_value(el,bf)     = (qmcfloat)chi;
              chi_grx(el,bf)       = (qmcfloat)chi_gradx;
              chi_gry(el,bf)       = (qmcfloat)chi_grady;
              chi_grz(el,bf)       = (qmcfloat)chi_gradz;
              chi_laplacian(el,bf) = (qmcfloat)chi_lap;
              bf++;
            }
        }
      el++;
    }
}

void QMCBasisFunction::
evaluateBasisFunctions(Array2D<double>& X, Array2D<qmcfloat>& chi_value)
{
  //This line helps prevent some floating point errors
  const double TOOSMALL = 1e-306;
  int el = 0, bf;
  int a, b, c, nGaussians;
  int numBF;
  double x, y, z, x2, y2, z2, r_sq, xyz;
  double p0, p1, exp_term, temp;
  double chi;
  for (int index=0; index<X.dim1(); index++)
    {
      bf = 0;
      for (int atom=0; atom<flags->Natoms; atom++)
        {
          for (int i=0; i<3; i++)
            {
              Xcalc(i) = X(index,i) - Molecule->Atom_Positions(atom,i);
            }
          x = Xcalc(0);
          y = Xcalc(1);
          z = Xcalc(2);

          x2 = x*x;
          y2 = y*y;
          z2 = z*z;
          r_sq = x2+y2+z2;
          numBF = BFCoeffs(atom).getNumberBasisFunctions();
          for (int j=0; j<numBF; j++)
            {
              a = BFCoeffs(atom).xyz_powers(j,0);
              b = BFCoeffs(atom).xyz_powers(j,1);
              c = BFCoeffs(atom).xyz_powers(j,2);

              xyz = pow(x,a)*pow(y,b)*pow(z,c);

              qmcfloat** coeffs = BFCoeffs(atom).Coeffs.array()[j];

              nGaussians = BFCoeffs(atom).N_Gauss(j);
              chi       = 0;
              for (int i=0; i<nGaussians; i++)
                {
                  p0 = coeffs[i][0];
                  p1 = coeffs[i][1];
                  exp_term = p1*exp(-p0*r_sq)*xyz;
                  temp = -2.0*p0;
                  chi       += exp_term;
                }

              //This block could safely be commented out for most compilers
              if(chi > 0 && chi < TOOSMALL) chi = TOOSMALL;
              if(chi < 0 && chi > -1.0*TOOSMALL) chi = -1.0*TOOSMALL;
              chi_value(el,bf)     = (qmcfloat)chi;
              //printf("chi: %3i %3i %10.7f\n",el,bf,chi);
              bf++;
            }//sum over basis functions
        }//sum over atoms
      el++;
    }//sum over index
}


