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

#include "QMCFunctions.h"

QMCFunctions::QMCFunctions()
{
}

QMCFunctions::QMCFunctions(QMCInput *INPUT)
{
  initialize(INPUT);
}

QMCFunctions::QMCFunctions(const QMCFunctions & rhs )
{
  *this = rhs;
}

void QMCFunctions::operator=(const QMCFunctions & rhs )
{
  Input = rhs.Input; 

  Alpha = rhs.Alpha; 
  Beta  = rhs.Beta; 

  Jastrow = rhs.Jastrow;  
  PE      = rhs.PE;
 
  Psi                    = rhs.Psi;
  Laplacian_PsiRatio     = rhs.Laplacian_PsiRatio;
  Grad_PsiRatio          = rhs.Grad_PsiRatio;
  Modified_Grad_PsiRatio = rhs.Modified_Grad_PsiRatio;
  E_Local                = rhs.E_Local;
}

void QMCFunctions::initialize(QMCInput *INPUT)
{
  Input = INPUT;

  Alpha.initialize(Input,0,Input->WF.getNumberAlphaElectrons()-1);
  Beta.initialize(Input,Input->WF.getNumberAlphaElectrons(),
		  Input->WF.getNumberElectrons()-1);

  PE.initialize(Input);
  Jastrow.initialize(Input);
  
  Grad_PsiRatio.allocate(Input->WF.getNumberElectrons(),3);
}

void QMCFunctions::evaluate(Array2D<double> &X)
{
  Alpha.evaluate(X);
  Beta.evaluate(X);

  Jastrow.evaluate(X);
  PE.evaluate(X);

  calculate_Psi();
  calculate_Grad_PsiRatio();
  calculate_Modified_Grad_PsiRatio(X);
  calculate_Laplacian_PsiRatio();
  calculate_E_Local();
}

void QMCFunctions::calculate_Psi()
{
  Psi = Alpha.getPsi() * Beta.getPsi() * Jastrow.getJastrow();
}

void QMCFunctions::calculate_Grad_PsiRatio()
{
  for(int i=0; i<Input->WF.getNumberAlphaElectrons(); i++)
    {
      Array2D<double> * AlphaGrad = Alpha.getGradPsiRatio();
      for(int j=0; j<3; j++)
	{
	  Grad_PsiRatio(i,j) = (*AlphaGrad)(i,j);
	}
    }

  for(int i=0; i<Input->WF.getNumberBetaElectrons(); i++)
    {
      Array2D<double> * BetaGrad = Beta.getGradPsiRatio();
      for(int j=0; j<3; j++)
	{
	  Grad_PsiRatio(i+Input->WF.getNumberAlphaElectrons(),j) = 
	    (*BetaGrad)(i,j);
	}
    }

  Array2D<double> * grad_jastrow_ratio = Jastrow.getGradientLnJastrow();

  for(int i=0; i<Input->WF.getNumberElectrons(); i++)
    {
      for(int j=0; j<3; j++)
	{
	  Grad_PsiRatio(i,j) += (*grad_jastrow_ratio)(i,j);
	}
    }
}

void QMCFunctions::calculate_Modified_Grad_PsiRatio(Array2D<double> &X)
{
  // Call this after calculate_Grad_PsiRatio() is called

  if( Input->flags.QF_modification_type == "none" )
    {
      // do not modify the quantum force
      Modified_Grad_PsiRatio = Grad_PsiRatio;
    }
  else if( Input->flags.QF_modification_type == "umrigar93_equalelectrons" )
    {
      // from Umrigar, Nightingale, and Runge JCP 99(4) 2865; 1993 eq 34.
      // The QF for each electron is changed by the same factor

      double a = Input->flags.umrigar93_equalelectrons_parameter;

      if( a>1 || a<=0 )
	{
	  cerr << "ERROR: Improper value for " 
	       << "umrigar93_equalelectrons_parameter!\n The value should be"
	       << " between 0 and 1 inclusive of 1 and exclusive of 0." 
	       << endl;
	  exit(0);
	}

      // Calculate the magnitude squared of the Grad_PsiRatio
      double magsqQF = 0.0;
      for(int i=0; i<Grad_PsiRatio.dim1(); i++)
	{
	  for(int j=0; j<Grad_PsiRatio.dim2(); j++)
	    {
	      magsqQF += Grad_PsiRatio(i,j) * 
		Grad_PsiRatio(i,j);
	    }
	}

      double factor = ( -1.0 + sqrt( 1.0 + 2*a*magsqQF*Input->flags.dt )) /
	( a*magsqQF*Input->flags.dt);

      Modified_Grad_PsiRatio = Grad_PsiRatio;
      Modified_Grad_PsiRatio *= factor;
    }
  else if( Input->flags.QF_modification_type == "umrigar93_unequalelectrons" )
    {
      // from Umrigar, Nightingale, and Runge JCP 99(4) 2865; 1993 eq 35-36.
      // The QF for each electron is changed by a different factor

      Modified_Grad_PsiRatio = Grad_PsiRatio;

      for(int i=0; i<Modified_Grad_PsiRatio.dim1(); i++)
	{
	  // Find the closest nucleus to electron i
	  int closest_nucleus_index  = 
	    Input->Molecule.findClosestNucleusIndex(X,i);
	  
	  double closest_nucleus_distance_squared = 0.0;
	  
	  for(int xyz=0; xyz<3; xyz++)
	    {
	      double temp = X(i,xyz) - 
		Input->Molecule.Atom_Positions(closest_nucleus_index,xyz);
	      
	      closest_nucleus_distance_squared += temp * temp;
	    }

	  // now we have the closest nucleus identified

	  // Charge of this nucleus
	  int closest_nucleus_Z = 
	    Input->Molecule.Z(closest_nucleus_index);

	  // Unit vector from nearest nucleus to the electron
	  Array1D<double> closest_nucleus_to_electron_norm_vector(3);
	  for(int xyz=0; xyz<3; xyz++)
	    {
	      closest_nucleus_to_electron_norm_vector(xyz) = 
		X(i,xyz) - 
		Input->Molecule.Atom_Positions(closest_nucleus_index,xyz);
	    }
	  closest_nucleus_to_electron_norm_vector *= 1.0/sqrt(
	       closest_nucleus_to_electron_norm_vector * 
	       closest_nucleus_to_electron_norm_vector);

	  // Unit vector in the direction of electron velocity
	  Array1D<double> electron_velocity_norm_vector(3);
	  for(int xyz=0; xyz<3; xyz++)
	    {
	      electron_velocity_norm_vector(xyz) = 
		Grad_PsiRatio(i,xyz);
	    }
	  double electron_velocity_norm_squared =
	    electron_velocity_norm_vector * electron_velocity_norm_vector;
	  electron_velocity_norm_vector *= 1.0/sqrt( 
					    electron_velocity_norm_squared );

	  // Now calculate the a(r) for the equation
	  double temp = closest_nucleus_Z * closest_nucleus_Z * 
	    closest_nucleus_distance_squared;

	  double a = 0.5*(1.0 + closest_nucleus_to_electron_norm_vector * 
			  electron_velocity_norm_vector) + 
	    temp/(10.0*(4.0+temp));

	  // Now calculate the factor the QF of this electron is modified by
	  double factor = (-1.0 + sqrt(1.0+
	     2.0*a*electron_velocity_norm_squared*Input->flags.dt))/
	    (a*electron_velocity_norm_squared*Input->flags.dt);

	  // Now adjust the QF of electron I
	  for(int xyz=0; xyz<3; xyz++)
	    {
	      Modified_Grad_PsiRatio(i,xyz) *= factor;
	    }
	}
    }
  else
    {
      cerr << "ERROR: Unknown value for QF_modification_type set!" << endl;
      exit(0);
    }
}

void QMCFunctions::calculate_Laplacian_PsiRatio()
{
  // must calculate Grad_PsiRatio before calling this

  Laplacian_PsiRatio = Alpha.getLaplacianPsiRatio() +
    Beta.getLaplacianPsiRatio() + Jastrow.getLaplacianLnJastrow();

  Array2D<double> * grad_jastrow_ratio = Jastrow.getGradientLnJastrow();

  for(int i=0; i<Input->WF.getNumberElectrons(); i++)
    {
      for(int j=0; j<3; j++)
	{
	  Laplacian_PsiRatio += (*grad_jastrow_ratio)(i,j) *
	    (2*Grad_PsiRatio(i,j) - (*grad_jastrow_ratio)(i,j));
	}
    }
}

void QMCFunctions::calculate_E_Local()
{
  // must calculate Laplacian_PsiRatio before calling this

  E_Local = -0.5 * Laplacian_PsiRatio + PE.getEnergy();
}

double QMCFunctions::getPsi()
{
  return Psi;
}

double QMCFunctions::getLocalEnergy()
{
  return E_Local;
}

double QMCFunctions::getKineticEnergy()
{
  return -0.5 * Laplacian_PsiRatio;
}

double QMCFunctions::getPotentialEnergy()
{
  return PE.getEnergy();
}

Array2D<double> * QMCFunctions::getGradPsiRatio()
{
  return &Grad_PsiRatio;
}

Array2D<double> * QMCFunctions::getModifiedGradPsiRatio()
{
  return &Modified_Grad_PsiRatio;
}

void QMCFunctions::writeCorrelatedSamplingConfiguration(ostream& strm)
{
  // Print the sum of the laplacians from the alpha and beta electrons
  strm << "D1" << endl;
  strm << "\t" << (Alpha.getLaplacianPsiRatio()+ Beta.getLaplacianPsiRatio()) 
       << endl;

  // print the GradPsiRatio excluding the Jastrow
  strm << "D2" << endl;

  Array2D<double> * TempGrad = Alpha.getGradPsiRatio();

  for(int i=0; i<Input->WF.getNumberAlphaElectrons(); i++)
    {
      for(int j=0; j<3; j++)
	{
	  strm << "\t" << (*TempGrad)(i,j);
	}
      strm << endl;
    }

  TempGrad = Beta.getGradPsiRatio();
  for(int i=0; i<Input->WF.getNumberBetaElectrons(); i++)
    {
      for(int j=0; j<3; j++)
 	strm << "\t" << (*TempGrad)(i,j);
      strm << endl;
    }
  
  strm << "J\t" << endl;
  strm << "\t" << Jastrow.getJastrow() << endl;
  strm << "PE\t" << endl;
  strm << "\t" << PE.getEnergy() << endl;
}

bool QMCFunctions::isSingular()
{
  if( Alpha.isSingular() || Beta.isSingular() )
    {
      return true;
    }
  else
    {
      return false;
    }
}





