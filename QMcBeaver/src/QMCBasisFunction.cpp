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
#include "IeeeMath.h"
#include "MathFunctions.h"

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

void QMCBasisFunction::evaluateBasisFunctions(Array2D<double>& X,
					      int start, int stop,
					      Array2D<qmcfloat>& chi_val,
					      Array2D<qmcfloat>& chi_grx,
					      Array2D<qmcfloat>& chi_gry,
					      Array2D<qmcfloat>& chi_grz,
					      Array2D<qmcfloat>& chi_lap)
{
  int el = 0, bf;
  int a, b, c;

  qmcfloat r_sq;
  qmcfloat p0, exp_term, temp1, temp2, temp3x, temp3y, temp3z;
  qmcfloat x_val, x_grx, x_gry, x_grz, x_lap;

  temp1 = 0;
  temp2 = 0;
  temp3x = 0;
  temp3y = 0;
  temp3z = 0;
  for (int idx=start; idx<=stop; idx++)
    {
      bf = 0;
      p0 = 0.0;
      for (int atom=0; atom<flags->Natoms; atom++)
        {
          for (int xyz=0; xyz<3; xyz++)
	    Xcalc(xyz) = X(idx,xyz) - Molecule->Atom_Positions(atom,xyz);
	  
	  const int numl = BFCoeffs(atom).lmax + 1;
	  const int numn = max(3,numl);
	  /*
	    xyz_abc[x,y, or z][l][n];
	    if n = 0: then x^n
	    if n = 1: then l / x
	    if n = 2: then l(l-1) / x^2
	  */
	  double xyz_abc[3][numn][3];
	  for(int xyz=0;xyz<3;xyz++)
	    {
	      xyz_abc[xyz][0][0] = 1.0;
	      xyz_abc[xyz][1][0] = Xcalc(xyz);
	      for(int abc=2; abc<numn; abc++)
		xyz_abc[xyz][abc][0] = xyz_abc[xyz][1][0] * xyz_abc[xyz][abc-1][0];

	      xyz_abc[xyz][0][1] = 0.0;
	      xyz_abc[xyz][0][2] = 0.0;
	      for(int abc=1; abc<numn; abc++)
		{
		  xyz_abc[xyz][abc][1] = abc / xyz_abc[xyz][1][0];
		  xyz_abc[xyz][abc][2] = abc*(abc-1.0)/xyz_abc[xyz][2][0];
		}
	    }
          r_sq = xyz_abc[0][2][0]+xyz_abc[1][2][0]+xyz_abc[2][2][0];
	  
	  double xyzA[numl][numl][numl][2];
	  for(int ai=0; ai<numl; ai++)
	    for(int bi=0; bi<numl; bi++)
	      for(int ci=0; ci<numl; ci++)
		{
		  int sum = ai+bi+ci;
		  if(sum>numl) continue;
		  xyzA[ai][bi][ci][0] =
		    xyz_abc[0][ai][0]*
		    xyz_abc[1][bi][0]*
		    xyz_abc[2][ci][0];
		  xyzA[ai][bi][ci][1] = 4.0*sum+6.0;
		}

          const int numBF = BFCoeffs(atom).getNumberBasisFunctions();
          for (int j=0; j<numBF; j++)
            {
              a = BFCoeffs(atom).xyz_powers(j,0);
              b = BFCoeffs(atom).xyz_powers(j,1);
              c = BFCoeffs(atom).xyz_powers(j,2);
              qmcfloat** coeffs = BFCoeffs(atom).Coeffs.array()[j];

              x_val = 0;
              x_grx = 0;
              x_gry = 0;
              x_grz = 0;
              x_lap = 0;
              const int nGaussians = BFCoeffs(atom).N_Gauss(j);
              for (int i=0; i<nGaussians; i++)
                {
		  if(fabs(p0 - coeffs[i][0]) > 1e-10){
		    p0     = coeffs[i][0];
		    temp1  = exp(-p0*r_sq);		    
		    temp2  = 4.0*r_sq*p0*p0;
		    temp3x = -2.0*p0*xyz_abc[0][1][0];
		    temp3y = -2.0*p0*xyz_abc[1][1][0];
		    temp3z = -2.0*p0*xyz_abc[2][1][0];
		  }

                  exp_term = temp1*coeffs[i][1]*xyzA[a][b][c][0];

                  x_val += exp_term;
                  x_grx += exp_term * (xyz_abc[0][a][1] + temp3x);
                  x_gry += exp_term * (xyz_abc[1][b][1] + temp3y);
                  x_grz += exp_term * (xyz_abc[2][c][1] + temp3z);
                  x_lap += exp_term * (xyz_abc[0][a][2] +
				       xyz_abc[1][b][2] +
				       xyz_abc[2][c][2] -
				       p0*xyzA[a][b][c][1]+temp2);
                }
	      
              chi_val(el,bf) = (qmcfloat)x_val;
              chi_grx(el,bf) = (qmcfloat)x_grx;
              chi_gry(el,bf) = (qmcfloat)x_gry;
              chi_grz(el,bf) = (qmcfloat)x_grz;
              chi_lap(el,bf) = (qmcfloat)x_lap;
              bf++;
            }
        }
      el++;
    }
}

void QMCBasisFunction::evaluateBasisFunctions(Array2D<double>& X,
					      Array2D<qmcfloat>& chi_val)
{
  int bf;
  int a, b, c;

  qmcfloat r_sq;
  qmcfloat p0, temp1;
  qmcfloat x_val;

  temp1 = 0;
  for (int idx=0; idx<X.dim1(); idx++)
    {
      bf = 0;
      p0 = 0.0;
      for (int atom=0; atom<flags->Natoms; atom++)
        {
          for (int xyz=0; xyz<3; xyz++)
	    Xcalc(xyz) = X(idx,xyz) - Molecule->Atom_Positions(atom,xyz);
	  
	  const int numl = BFCoeffs(atom).lmax + 1;

	  /*
	    xyz_abc[x,y, or z][l][n];
	    if n = 0: then x^n
	    if n = 1: then l / x
	    if n = 2: then l(l-1) / x^2
	  */
	  double xyz_abc[3][numl];
	  for(int xyz=0;xyz<3;xyz++)
	    {
	      xyz_abc[xyz][0] = 1.0;
	      xyz_abc[xyz][1] = Xcalc(xyz);
	      for(int abc=2; abc<numl; abc++)
		xyz_abc[xyz][abc] = xyz_abc[xyz][1] * xyz_abc[xyz][abc-1];
	    }
          r_sq = xyz_abc[0][2]+xyz_abc[1][2]+xyz_abc[2][2];
	  
	  double xyzA[numl][numl][numl];
	  for(int ai=0; ai<numl; ai++)
	    for(int bi=0; bi<numl; bi++)
	      for(int ci=0; ci<numl; ci++)
		{
		  int sum = ai+bi+ci;
		  if(sum>numl) continue;
		  xyzA[ai][bi][ci] =
		    xyz_abc[0][ai]*
		    xyz_abc[1][bi]*
		    xyz_abc[2][ci];
		}

          const int numBF = BFCoeffs(atom).getNumberBasisFunctions();
          for (int j=0; j<numBF; j++)
            {
              a = BFCoeffs(atom).xyz_powers(j,0);
              b = BFCoeffs(atom).xyz_powers(j,1);
              c = BFCoeffs(atom).xyz_powers(j,2);
              qmcfloat** coeffs = BFCoeffs(atom).Coeffs.array()[j];

              x_val = 0;
              const int nGaussians = BFCoeffs(atom).N_Gauss(j);
              for (int i=0; i<nGaussians; i++)
                {
		  if(fabs(p0 - coeffs[i][0]) > 1e-10){
		    p0     = coeffs[i][0];
		    temp1  = exp(-p0*r_sq);		    
		    //temp1  = MathFunctions::fast_exp(-p0*r_sq);
		  }

                  x_val += temp1*coeffs[i][1];
                }

              chi_val(idx,bf) = (qmcfloat)x_val*xyzA[a][b][c];
              bf++;
            }
        }
    }
}

void QMCBasisFunction::basisFunctionsOnGrid(Array2D<double>& grid,
					    int nuc,
					    Array2D< Array1D<double> > & angularCi,
					    Array2D<qmcfloat>& chi_val)
{
  int bf;

  qmcfloat p0, temp1;
  temp1 = 0;
  //First just handle the radial gaussians: X = \sum pi1 exp(-pi0 r)
  for (int idx=0; idx<grid.dim1(); idx++)
    {
      bf = 0;
      p0 = 0.0;
      for (int atom=0; atom<flags->Natoms; atom++)
        {
	  if(atom == nuc && idx > 0){

	    const int numBF = BFCoeffs(atom).getNumberBasisFunctions();
	    for (int j=0; j<numBF; j++)
	      {
		// All of the grid points are the same distance from this nucleus, so the radial
		// parts are identical
		chi_val(idx,bf) = chi_val(0,bf);		
		bf++;
	      }

	  } else {

	    double r_sq = 0.0;
	    for (int xyz=0; xyz<3; xyz++){
	      Xcalc(xyz) = grid(idx,xyz) - Molecule->Atom_Positions(atom,xyz);
	      r_sq += Xcalc(xyz)*Xcalc(xyz);
	    }
	    
	    const int numBF = BFCoeffs(atom).getNumberBasisFunctions();	   
	    for (int j=0; j<numBF; j++)
	      {
		qmcfloat** coeffs = BFCoeffs(atom).Coeffs.array()[j];
		
		double x_val = 0;
		const int nGaussians = BFCoeffs(atom).N_Gauss(j);
		for (int i=0; i<nGaussians; i++)
		  {
		    if(fabs(p0 - coeffs[i][0]) > 1e-10){
		      p0       = coeffs[i][0];
		      double x = -p0*r_sq;
		      //temp1    = MathFunctions::fast_exp(x); 
		      temp1  = exp(x);		      
		    }
		    
		    x_val += temp1*coeffs[i][1];
		  }
		
		chi_val(idx,bf) = (qmcfloat)x_val;
		bf++;
	      }
	  }
        }
    }

  //Now add in the angular components x^a y^b z^c
  double r = 0.0;
  for (int xyz=0; xyz<3; xyz++)
    {
      Xcalc(xyz) = grid(0,xyz) - Molecule->Atom_Positions(nuc,xyz);
      r += Xcalc(xyz)*Xcalc(xyz);
    }
  r = sqrt(r);

  const int maxrn = 4;
  double rn[maxrn];
  rn[0] = 1.0;
  for(int i=1; i<maxrn; i++)
    rn[i] = rn[i-1]*r;


  Array2D<double> ang(angularCi.dim1(),angularCi.dim2());
  ang = 0.0;
  double * chi_p = chi_val.array();
  Array1D<double> * angci_p = angularCi.array();
  for(int i=0; i<angularCi.dim1()*angularCi.dim2(); i++)
    {
      double * ang_p = angci_p[i].array();
      register double sum = ang_p[0];
      for(int l=1; l<angci_p[i].dim1(); l++)
	sum += ang_p[l] * rn[l];
      chi_p[i] *= sum;
    }
}


void QMCBasisFunction::angularGrid(Array2D<double>& grid,
				   int nuc,
				   Array2D< Array1D<double> > & angularCi)
{
  int bf;
  int a, b, c;

  angularCi.allocate(grid.dim1(),N_BasisFunctions);

  for (int idx=0; idx<grid.dim1(); idx++)
    {
      bf = 0;

      double x = grid(idx,0);
      double y = grid(idx,1);
      double z = grid(idx,2);
      //double r = sqrt(x*x+y*y+z*z);//this should be 1.0

      for (int atom=0; atom<flags->Natoms; atom++)
        {
	  double x0 = Molecule->Atom_Positions(atom,0) - Molecule->Atom_Positions(nuc,0);
	  double y0 = Molecule->Atom_Positions(atom,1) - Molecule->Atom_Positions(nuc,1);
	  double z0 = Molecule->Atom_Positions(atom,2) - Molecule->Atom_Positions(nuc,2);

          const int numBF = BFCoeffs(atom).getNumberBasisFunctions();
          for (int j=0; j<numBF; j++)
            {
              a = BFCoeffs(atom).xyz_powers(j,0);
              b = BFCoeffs(atom).xyz_powers(j,1);
              c = BFCoeffs(atom).xyz_powers(j,2);
	      const int numl = a+b+c+1;	      
	      Array1D<double> poly(numl);
	      poly = 0.0;
	      poly(0) = 1.0;

	      bool print = !true;
	      if(idx != 0) print = false;
	      if(print){
		if(a > 0)	printf("(%8.5f r-%8.5f)^%i ",x,x0,a);
		else	printf("%24s","");
		if(b > 0)	printf("(%8.5f r-%8.5f)^%i ",y,y0,b);
		else	printf("%24s","");
		if(c > 0)	printf("(%8.5f r-%8.5f)^%i ",z,z0,c);
		else	printf("%24s","");
		printf(" %i%i%i = ",a,b,c);
	      }

	      for(int ai=0; ai<a; ai++)
		{
		  Array1D<double> temp = poly;
		  for(int l=1; l<numl; l++)
		    poly(l) = x*temp(l-1) - x0*temp(l);
		  poly(0) = - x0*temp(0);
		}

	      for(int bi=0; bi<b; bi++)
		{
		  Array1D<double> temp = poly;
		  for(int l=1; l<numl; l++)
		    poly(l) = y*temp(l-1) - y0*temp(l);
		  poly(0) = - y0*temp(0);
		}

	      for(int ci=0; ci<c; ci++)
		{
		  Array1D<double> temp = poly;
		  for(int l=1; l<numl; l++)
		    poly(l) = z*temp(l-1) - z0*temp(l);
		  poly(0) = - z0*temp(0);
		}

	      int firstNonZero = 0;
	      for(int i=poly.dim1()-1; i>=0; i--)
		if(fabs(poly(i)) > 1e-15) {
		  firstNonZero = i;
		  break;
		}

	      Array1D<double> temp(firstNonZero+1);
	      for(int i=0; i<temp.dim1(); i++)
		temp(i) = poly(i);
	      poly = temp;
	      
	      if(print) cout << poly << endl;
	      angularCi(idx,bf) = poly;
              bf++;
            }//bf
        }//atom
    }//grid
}


