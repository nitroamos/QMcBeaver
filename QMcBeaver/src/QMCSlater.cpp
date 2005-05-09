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
  Input              = rhs.Input;
  Psi                = rhs.Psi;
  Laplacian_PsiRatio = rhs.Laplacian_PsiRatio;
  Grad_PsiRatio      = rhs.Grad_PsiRatio;
  if (Input->flags.calculate_bf_density == 1)
    Chi_Density      = rhs.Chi_Density;
  
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
  allocate( Stop-Start+1 );
}

void QMCSlater::initialize
(QMCInput *INPUT, int startEl, int stopEl, Array2D<int> occ)
{
  Input = INPUT;
  BF = &Input->BF;
  WF = &Input->WF;
    
  setStartAndStopElectronPositions(startEl, stopEl); 
  occupation = occ;

  for(int i=0; i<Singular.dim1(); i++)
    Singular(i) = false;

  int orbital_index;
  int numRows;
  int numOrbs = WF->getNumberOrbitals();
  for(int i=0; i<WF->getNumberDeterminants(); i++)
    {
      numRows = 0;
      orbital_index = 0;
      for(int j=0; j<=numOrbs; j++)
	{
	  if (j<numOrbs && occupation(i,j) == 1)
	    {
	      numRows++;
	      orbital_index++;
	    } else {
	      if(numRows != 0)
		{
		  WF_coeffs(i).setRows(orbital_index-numRows,j-numRows,numRows,WF->Coeffs);
		  numRows = 0;
		}
	      if (orbital_index == D(0,i).dim2()) break;
	    }
	}
    }
}

void QMCSlater::allocate(int N)
{
  int ndet = WF->getNumberDeterminants();  
  int nbasisfunc = WF->getNumberBasisFunctions();

  D.allocate(WALKERS_PER_PASS,ndet);
  D_inv.allocate(WALKERS_PER_PASS,ndet);
  Laplacian_D.allocate(WALKERS_PER_PASS,ndet);
  Grad_D.allocate(WALKERS_PER_PASS,ndet,3);

  Singular.allocate(WALKERS_PER_PASS);      
  Psi.allocate(WALKERS_PER_PASS);
  Laplacian_PsiRatio.allocate(WALKERS_PER_PASS);
  Grad_PsiRatio.allocate(WALKERS_PER_PASS);
  
  if (Input->flags.calculate_bf_density == 1)
    Chi_Density.allocate(WALKERS_PER_PASS);

  for(int j=0; j<WALKERS_PER_PASS; j++)
    {      
      for (int i=0; i<ndet; i++)
	{
	  D(j,i).allocate(N,N);
	  D_inv(j,i).allocate(N,N);
	  Laplacian_D(j,i).allocate(N,N);
	  for(int k=0; k<3; k++)
	    Grad_D(j,i,k).allocate(N,N);
	}
      
      Singular(j).allocate(ndet);
      Psi(j).allocate(ndet);
      Laplacian_PsiRatio(j).allocate(ndet);
      Grad_PsiRatio(j).allocate(ndet,N,3);
      if (Input->flags.calculate_bf_density == 1)
	Chi_Density(j).allocate(nbasisfunc);      
    }

  WF_coeffs.allocate(ndet);
  for(int i=0; i<ndet; i++)
    WF_coeffs(i).allocate(N,nbasisfunc);

  occupation.allocate(ndet,WF->getNumberOrbitals());
  
  Chi.allocate(N,nbasisfunc);
  Chi_laplacian.allocate(N,nbasisfunc);
  Chi_gradient.allocate(3);
  for(int j=0; j<3; j++)
    Chi_gradient(j).allocate(N,nbasisfunc);
}

QMCSlater::~QMCSlater()
{
  for(int j=0; j < WALKERS_PER_PASS; j++){
    for (int i=0; i<WF->getNumberDeterminants(); i++)
      {
	D(j,i).deallocate();
	D_inv(j,i).deallocate();
	Laplacian_D(j,i).deallocate();
	for (int k=0; k<3; k++)
	  Grad_D(j,i,k).deallocate();
      }
    Singular(j).deallocate();
    Psi(j).deallocate();
    Laplacian_PsiRatio(j).deallocate();
    Grad_PsiRatio(j).deallocate();    
    if (Input->flags.calculate_bf_density == 1)
      Chi_Density(j).deallocate();
  }

  D.deallocate();
  D_inv.deallocate();
  Laplacian_D.deallocate();
  Grad_D.deallocate();
  Singular.deallocate();
  Psi.deallocate();
  Laplacian_PsiRatio.deallocate();
  Grad_PsiRatio.deallocate();    
  
  occupation.deallocate();

  if (Input->flags.calculate_bf_density == 1)
    Chi_Density.deallocate();

  Chi.deallocate();
  Chi_laplacian.deallocate();
  for (int j=0; j<3; j++)
    Chi_gradient(j).deallocate();
  Chi_gradient.deallocate();

  WF_coeffs.deallocate();
}

void QMCSlater::setStartAndStopElectronPositions(int StartEl, int StopEl)
{
  Start = StartEl;
  Stop  = StopEl;
  
  allocate(StopEl-StartEl+1);
}

void QMCSlater::evaluate(Array1D<Array2D<double>*> &X, int num)
{
  update_Ds(X, num);
  
  update_D_inverse_and_Psi(num);
  
  for(int i=0; i<num; i++){
    if( !isSingular(i) )
      {
	calculate_DerivativeRatios(i);
      }
  }
}

void QMCSlater::update_D_inverse_and_Psi(int num)
{
  // The LU Decomposition inverse used here is O(1*N^3)
  // Updating one electron at a time is O(2*N^3-N^2)
  
  bool calcOK = true;
  
  for(int j=0; j<num; j++)
    for (int i=0; i<WF->getNumberDeterminants(); i++)
      {
	determinant_and_inverse(D(j,i),D_inv(j,i),(Psi(j))(i),&calcOK);
	(Singular(j))(i) = !calcOK;
      }
}

/**
   This function contains the meat of the QMC calculation.
   1) basis functions are calculated for each electron position
   for each determinant:
   2) a wavefunction coefficient matrix is composed of the input basis function coefficients
   some monkey business is required here because not all consecutive orbitals in the
   input coefficients are occupied.
   3) the Slater matrix is calculated along with the associated gradient and laplacian
   matricies:
   D(numElec x numOrb) = Chi(numElec x numBasisFunction) * WF_coeffs(numBasisFunction x numOrb)
   
   Note: the coefficient matricies are transposed. hand-coded matrix multiplication is faster
   with a transposed matrix, and it enables the use of memcpy for Array2D's row-major data.
*/
void QMCSlater::update_Ds(Array1D<Array2D<double>*> &X, int num)
{
  for(int walker = 0; walker < num; walker++)
    {
      BF->evaluateBasisFunctions(*X(walker),Start,Stop,
				 Chi,
				 Chi_gradient(0),
				 Chi_gradient(1),
				 Chi_gradient(2),
				 Chi_laplacian);
      for(int i=0; i<WF->getNumberDeterminants(); i++)
	{      
	  D(walker,i)           = Chi * WF_coeffs(i);
	  Laplacian_D(walker,i) = Chi_laplacian * WF_coeffs(i);
	  Grad_D(walker,i,0)    = Chi_gradient(0) * WF_coeffs(i);
	  Grad_D(walker,i,1)    = Chi_gradient(1) * WF_coeffs(i);
	  Grad_D(walker,i,2)    = Chi_gradient(2) * WF_coeffs(i);
	}

      if (Input->flags.calculate_bf_density == 1)
	{
	  Chi_Density(walker) = 0.0;
	  for (int i=0; i<WF->getNumberBasisFunctions(); i++)
	    for (int j=0; j<D(0,0).dim1(); j++)
	      {
		(Chi_Density(walker))(i) += Chi(j,i);
	      }
	} 
    }
}

/**
   This function calculates the derivative ratios. that is, del psi/psi, and lap psi/psi.
   D_inv is the transpose of the true inverse, but it turns out that if we use the transpose
   here, then the calculation of the ratios turns into a sort of dot product, which not only does ATLAS
   know how to do it, it makes the code look neater. the hand-coded dot product probably doesn't
   save much time.
   At this time, D_inv and Grad_D and Laplacian_D are all qmcfloat type. The explicit typecast
   when creating the final result (double) should emphasize this.
*/
void QMCSlater::calculate_DerivativeRatios(int k)
{
  int numElectrons = (D(0,0)).dim1();
  for(int i=0; i<WF->getNumberDeterminants(); i++)
    {
      double** grad_psiratioArray = Grad_PsiRatio(k).array()[i];
      
      (Laplacian_PsiRatio(k))(i) = (double)((Laplacian_D(k,i)).dotAllElectrons(D_inv(k,i)));
      
      for(int j=0; j<numElectrons; j++)
	{
	  grad_psiratioArray[j][0] = (double)((Grad_D(k,i,0)).dotOneElectron(D_inv(k,i),j));
	  grad_psiratioArray[j][1] = (double)((Grad_D(k,i,1)).dotOneElectron(D_inv(k,i),j));
	  grad_psiratioArray[j][2] = (double)((Grad_D(k,i,2)).dotOneElectron(D_inv(k,i),j));
	}
    }
}

Array1D<double>* QMCSlater::getPsi(int i)
{
  return &Psi(i);
}

Array1D<double>* QMCSlater::getLaplacianPsiRatio(int i)
{
  return &Laplacian_PsiRatio(i);
}

Array3D<double>* QMCSlater::getGradPsiRatio(int i)
{
  return &Grad_PsiRatio(i);
}

Array1D<double>* QMCSlater::getChiDensity(int i)
{
  return &Chi_Density(i);
}

bool QMCSlater::isSingular(int j)
{
  for (int i=0; i<WF->getNumberDeterminants(); i++)
    {
      if (Singular(j)(i)) return true;
    }
  return false;
}
