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

  Singular = rhs.Singular;

  // This was done in this round about way to avoid some weird compiler errors
  // dealing with when const was used.  In the simplest formt this would be
  // allocate( rhs.D.dim1() );
  D = rhs.D;
  int size = D.dim1();
  allocate( size );
}

void QMCSlater::initialize(QMCInput *INPUT, int startEl, int stopEl)
{
  Input = INPUT;
  BF = &Input->BF;
  WF = &Input->WF;

  Singular = false;

  setStartAndStopElectronPositions(startEl, stopEl); 
}

void QMCSlater::allocate(int N)
{
 D.allocate(N,N);
 D_inv.allocate(N,N);
 Laplacian_D.allocate(N,N);
 Grad_D.allocate(N,N,3);
 Grad_PsiRatio.allocate(N,3);
}

void QMCSlater::setStartAndStopElectronPositions(int StartEl, int StopEl)
{
  Start = StartEl;
  Stop  = StopEl;

  allocate(StopEl-StartEl+1);
}

void QMCSlater::evaluate(Array2D<double> &X)
{
  for(int electron=Start; electron<=Stop; electron++)
    {
      update_Ds(electron,X);
    }
  
  update_D_inverse_and_Psi();

  if( !isSingular() )
    {
      calculate_DerivativeRatios();
    }
}

void QMCSlater::update_D_inverse_and_Psi()
{
  // The LU Decomposition inverse used here is O(1*N^3)
  // Updating one electron at a time is O(2*N^3-N^2)

  bool calcOK = true;

  determinant_and_inverse(D,D_inv,Psi,&calcOK);

  Singular = !calcOK;
}

void QMCSlater::update_Ds(int electron, Array2D<double> &X)
{
  update_D(electron,X);
  update_Laplacian_D(electron,X);
  update_Grad_D(electron,X);
}

void QMCSlater::update_D(int electron, Array2D<double> &X)
{
  Array1D<double> Chi(WF->getNumberBasisFunctions());

  // calculate the value of the basis functions evaluated at the
  // current electron position
  for(int i=0; i<WF->getNumberBasisFunctions(); i++)
    {
      Chi(i) = BF->getPsi(i,X,electron);
    }

  //Calculate the new orbital values at the trial position
  for(int i=0; i<D.dim2(); i++)
    {
      D(electron-Start,i) = 0;
      for(int k=0; k<WF->getNumberBasisFunctions(); k++)
	{
	  D(electron-Start,i) += WF->Coeffs(k,i)
	    *Chi(k);
	}
    }
}  

void QMCSlater::update_Laplacian_D(int electron, Array2D<double> &X)
{
  Array1D<double> Chi(WF->getNumberBasisFunctions());

  // calculate the value of the basis functions evaluated at the
  // current electron position
  for(int i=0; i<WF->getNumberBasisFunctions(); i++)
    {
      Chi(i) = BF->getLaplacianPsi(i,X,electron);
    }

  //Calculate the new orbital values at the trial position
  for(int i=0; i<Laplacian_D.dim2(); i++)
    {
      Laplacian_D(electron-Start,i) = 0;
      for(int k=0; k<WF->getNumberBasisFunctions(); k++)
	{
	  Laplacian_D(electron-Start,i) += WF->Coeffs(k,i)
	    *Chi(k);
	}
    }
}  

void QMCSlater::update_Grad_D(int electron, Array2D<double> &X)
{
  Array2D<double> Chi(WF->getNumberBasisFunctions(),3);
  Array1D<double> Temp(3);

  //calculate the value of the basis functions evaluated at the
  //trial electron position
  for(int i=0; i<WF->getNumberBasisFunctions(); i++)
    {
      Temp = BF->getGradPsi(i,X,electron);
      for(int j=0; j<3; j++)
	Chi(i,j) = Temp(j);
    }

  //Calculate the new orbital values at the trial position
  for(int i=0; i<Grad_D.dim2(); i++)
    {
      for(int j=0; j<3; j++)
	{
	  Grad_D(electron-Start,i,j) = 0;
	  for(int k=0; k<WF->getNumberBasisFunctions(); k++)
	    {
	      Grad_D(electron-Start,i,j) += WF->Coeffs(k,i)
		*Chi(k,j);
	    }
	}
    }
}  

void QMCSlater::calculate_DerivativeRatios()
{
  calculate_Grad_PsiRatio();
  calculate_Laplacian_PsiRatio();
}

void QMCSlater::calculate_Laplacian_PsiRatio()
{
  Laplacian_PsiRatio = 0;

  for(int el=0; el<D.dim1(); el++)
    {
      for(int i=0;i<D.dim1();i++)
	{
	  Laplacian_PsiRatio += 
	    Laplacian_D(el,i)*
	    D_inv(i,el);
	}
    }
}

void QMCSlater::calculate_Grad_PsiRatio()
{
  for(int el=0; el<D.dim1(); el++)
    {
      for(int j=0;j<3;j++)
	{
	  Grad_PsiRatio(el,j) = 0.0;
	  for(int i=0;i<D.dim1();i++)
	    {
	      Grad_PsiRatio(el,j) += 
		Grad_D(el,i,j)*
		D_inv(i,el);
	    }
	}
    }
}

double QMCSlater::getPsi()
{
  return Psi;
}

double QMCSlater::getLaplacianPsiRatio()
{
  return Laplacian_PsiRatio;
}

Array2D<double> * QMCSlater::getGradPsiRatio()
{
  return &Grad_PsiRatio;
}

bool QMCSlater::isSingular()
{
  return Singular;
}






