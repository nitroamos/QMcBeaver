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


#include "QMCSlater.h"

void QMCSlater::operator=(const QMCSlater & rhs )
{
  Psi                = rhs.Psi;
  Laplacian_PsiRatio = rhs.Laplacian_PsiRatio;
  Grad_PsiRatio      = rhs.Grad_PsiRatio;

  Input = rhs.Input;
  BF    = rhs.BF;
  WF    = rhs.WF;

  Start          = rhs.Start;
  Stop           = rhs.Stop;
  occupation     = rhs.occupation;

  Singular = rhs.Singular;

  // This was done in this round about way to avoid some weird compiler errors
  // dealing with when const was used.  In the simplest format this would be
  // allocate( rhs.D.dim1() );
  D = rhs.D;
  int size = D(0).dim1();
  allocate( size );
}

void QMCSlater::initialize
                   (QMCInput *INPUT, int startEl, int stopEl, Array2D<int> occ)
{
  Input = INPUT;
  BF = &Input->BF;
  WF = &Input->WF;

  Singular = false;

  setStartAndStopElectronPositions(startEl, stopEl); 
  occupation = occ;
}

void QMCSlater::allocate(int N)
{
  D.allocate(WF->getNumberDeterminants());
  D_inv.allocate(WF->getNumberDeterminants());
  for (int i=0; i<WF->getNumberDeterminants(); i++)
    {
      D(i).allocate(N,N);
      D_inv(i).allocate(N,N);
    }

  Laplacian_D.allocate(WF->getNumberDeterminants(),N,N);
  Grad_D.allocate(WF->getNumberDeterminants(),N,N,3);
  Singular.allocate(WF->getNumberDeterminants());

  occupation.allocate(WF->getNumberDeterminants(),WF->getNumberOrbitals());

  Psi.allocate(WF->getNumberDeterminants());
  Laplacian_PsiRatio.allocate(WF->getNumberDeterminants());
  Grad_PsiRatio.allocate(WF->getNumberDeterminants(),N,3);

  Chi1D.allocate(WF->getNumberBasisFunctions());
  Chi1D_laplacian.allocate(WF->getNumberBasisFunctions());
  Chi2D.allocate(WF->getNumberBasisFunctions(),3);
  Grad1e.allocate(3);
}

void QMCSlater::setStartAndStopElectronPositions(int StartEl, int StopEl)
{
  Start = StartEl;
  Stop  = StopEl;

  allocate(StopEl-StartEl+1);
}

void QMCSlater::evaluate(Array2D<double> &X)
{
  for (int electron=Start; electron<=Stop; electron++)
    {
      update_Ds(electron,X);
    }

  for (int i=0; i<WF->getNumberDeterminants(); i++)
    {
      update_D_inverse_and_Psi(i);
    }
  
  if( !isSingular() )
    {
      calculate_DerivativeRatios();
    }
}

void QMCSlater::update_D_inverse_and_Psi(int i)
{
  // The LU Decomposition inverse used here is O(1*N^3)
  // Updating one electron at a time is O(2*N^3-N^2)

  bool calcOK = true;

  determinant_and_inverse(D(i),D_inv(i),Psi(i),&calcOK);

  Singular(i) = !calcOK;
}

void QMCSlater::update_Ds(int electron, Array2D<double> &X)
{
  int nBasisFunc = WF->getNumberBasisFunctions();

  double *chi1DArray = Chi1D.array();
  double *chi1D_laplacianArray = Chi1D_laplacian.array();
  double **chi2DArray = Chi2D.array();

  for (int i=0; i<nBasisFunc; i++)
    {
      double value;
      double laplacian;

      BF->evaluateBasisFunction(i,X,electron,value,Grad1e,laplacian);

      chi1DArray[i] = value;
      chi1D_laplacianArray[i] = laplacian;
      for (int j=0; j<3; j++)
	{
	  chi2DArray[i][j] = Grad1e(j);
	}
    }

  int offset = electron-Start;
  double **wfCoeffsArray = WF->Coeffs.array();

  for (int i=0; i<WF->getNumberDeterminants(); i++)
    {
      double *dArray = D(i).array()[offset];
      double *laplacianDarray = Laplacian_D.array()[i][offset];
      double **gradDArray = Grad_D.array()[i][offset];

      for (int j=0; j<D(i).dim2(); j++)
	{
	  dArray[j] = 0.0;
	  laplacianDarray[j] = 0.0;

	  for (int k=0; k<3; k++)
	    {
	      gradDArray[j][k] = 0.0;
	    }
	}

      int orbital_index = 0;
      for (int j=0; j<WF->getNumberOrbitals(); j++)
	{
	  if (occupation(i,j) == 1)
	    {
	      for (int k=0; k<nBasisFunc; k++)
		{
		  dArray[orbital_index] += wfCoeffsArray[k][j] * chi1DArray[k];
		  laplacianDarray[orbital_index] += wfCoeffsArray[k][j] * \
                                                       chi1D_laplacianArray[k];

		  for (int l=0; l<3; l++)
		    {
		      gradDArray[orbital_index][l] += wfCoeffsArray[k][j] * \
                                                              chi2DArray[k][l];
		    }
		}
	      if (orbital_index == D(i).dim2()) break;
	      else orbital_index++;
	    }
	}
    }
}

void QMCSlater::calculate_DerivativeRatios()
{
  for (int i=0; i<WF->getNumberDeterminants(); i++)
    {
      double** d_invArray = D_inv(i).array();

      Laplacian_PsiRatio(i) = 0.0;
      double** laplacianDArray = Laplacian_D.array()[i];

      double** grad_psiratioArray = Grad_PsiRatio.array()[i];
      double*** grad_dArray = Grad_D.array()[i];

      for (int el=0; el<D(i).dim1(); el++)
	{
	  for (int j=0; j<D(i).dim1(); j++)
	    {
	      Laplacian_PsiRatio(i) += laplacianDArray[el][j] * \
		                                             d_invArray[j][el];
	    }

	  for (int k=0; k<3; k++)
	    {
	      grad_psiratioArray[el][k] = 0.0;
	      for (int l=0; l<D(i).dim1(); l++)
		{
		  grad_psiratioArray[el][k] += grad_dArray[el][l][k] * \
                                                             d_invArray[l][el];
		}
	    }
	}
    }
}

Array1D<double>* QMCSlater::getPsi()
{
  return &Psi;
}

Array1D<double>* QMCSlater::getLaplacianPsiRatio()
{
  return &Laplacian_PsiRatio;
}

Array3D<double>* QMCSlater::getGradPsiRatio()
{
  return &Grad_PsiRatio;
}

bool QMCSlater::isSingular()
{
  for (int i=0; i<WF->getNumberDeterminants(); i++)
    {
      if (Singular(i)) return true;
    }
  return false;
}






