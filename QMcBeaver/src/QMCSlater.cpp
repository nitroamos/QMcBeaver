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
  int ndet = WF->getNumberDeterminants();
  
  D.allocate(ndet);
  D_inv.allocate(ndet);
  Laplacian_D.allocate(ndet);
  Grad_D.allocate(ndet,3);
  
  for (int i=0; i<ndet; i++)
    {
      D(i).allocate(N,N);
      D_inv(i).allocate(N,N);
      Laplacian_D(i).allocate(N,N);
      for(int j=0; j<3; j++)
	Grad_D(i,j).allocate(N,N);
    }
  
  Singular.allocate(ndet);
  
  occupation.allocate(ndet,WF->getNumberOrbitals());
  
  Psi.allocate(ndet);
  Laplacian_PsiRatio.allocate(ndet);
  Grad_PsiRatio.allocate(ndet,N,3);
  
  int nbasisfunc = WF->getNumberBasisFunctions();
  
  Chi.allocate(N,nbasisfunc);
  Chi_laplacian.allocate(N,nbasisfunc);
  Chi_gradient.allocate(3);
  for(int j=0; j<3; j++)
    Chi_gradient(j).allocate(N,nbasisfunc);
  
  WF_coeffs.allocate(N,nbasisfunc);
}

void QMCSlater::setStartAndStopElectronPositions(int StartEl, int StopEl)
{
  Start = StartEl;
  Stop  = StopEl;
  
  allocate(StopEl-StartEl+1);
}

void QMCSlater::evaluate(Array2D<double> &X)
{
  update_Ds(X);
  
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
  
  for (int i=0; i<WF->getNumberDeterminants(); i++)
    {
      determinant_and_inverse(D(i),D_inv(i),Psi(i),&calcOK);
      Singular(i) = !calcOK;
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
void QMCSlater::update_Ds(Array2D<double> &X){
  int numOrbs;
  BF->evaluateBasisFunctions(X,Start,Stop,
			     Chi,
			     Chi_gradient(0),
			     Chi_gradient(1),
			     Chi_gradient(2),
			     Chi_laplacian);
  
  for(int i=0; i<WF->getNumberDeterminants(); i++){
    
    int orbital_index = 0;
    int numRows = 0;
    numOrbs = WF->getNumberOrbitals();
    for(int j=0; j<=numOrbs; j++){
      if (j<numOrbs && occupation(i,j) == 1){
	numRows++;
	orbital_index++;
      } else {
	if(numRows != 0){
	  WF_coeffs.setRows(orbital_index-numRows,j-numRows,numRows,WF->Coeffs);
	  numRows = 0;
	}
	if (orbital_index == D(i).dim2()) break;
      }
    }
    
    D(i)           = Chi * WF_coeffs;
    Laplacian_D(i) = Chi_laplacian * WF_coeffs;
    Grad_D(i,0)    = Chi_gradient(0) * WF_coeffs;
    Grad_D(i,1)    = Chi_gradient(1) * WF_coeffs;
    Grad_D(i,2)    = Chi_gradient(2) * WF_coeffs;
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
void QMCSlater::calculate_DerivativeRatios()
{
  int numElectrons = (D(0)).dim1();
  for(int i=0; i<WF->getNumberDeterminants(); i++)
    {
      double** grad_psiratioArray = Grad_PsiRatio.array()[i];
      
      Laplacian_PsiRatio(i) = (double)((Laplacian_D(i)).dotAllElectrons(D_inv(i)));
      
      for(int j=0; j<numElectrons; j++){
	grad_psiratioArray[j][0] = (double)((Grad_D(i,0)).dotOneElectron(D_inv(i),j));
	grad_psiratioArray[j][1] = (double)((Grad_D(i,1)).dotOneElectron(D_inv(i),j));
	grad_psiratioArray[j][2] = (double)((Grad_D(i,2)).dotOneElectron(D_inv(i),j));
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
