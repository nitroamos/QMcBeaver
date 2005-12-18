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

QMCSlater::QMCSlater()
{
  Start = 0;
  Stop = 0;
}

void QMCSlater::operator=(const QMCSlater & rhs )
{
  Input              = rhs.Input;
  Psi                = rhs.Psi;
  Laplacian_PsiRatio = rhs.Laplacian_PsiRatio;
  Grad_PsiRatio      = rhs.Grad_PsiRatio;
  if (Input->flags.calculate_bf_density == 1)
    Chi_Density      = rhs.Chi_Density;

  BF         = rhs.BF;
  WF         = rhs.WF;
  Start      = rhs.Start;
  Stop       = rhs.Stop;
  occupation = rhs.occupation;

  Singular = rhs.Singular;

  // This was done in this round about way to avoid some weird compiler errors
  // dealing with when const was used.  In the simplest format this would be
  // allocate( rhs.D.dim1() );
  D = rhs.D;
  allocate( Stop-Start+1 );
}

void QMCSlater::initialize(QMCInput *INPUT, int startEl, int stopEl, 
			   Array2D<int> *occ, Array2D<qmcfloat> *coeffs)
{
  Input = INPUT;
  BF = &Input->BF;
  WF = &Input->WF;

  setStartAndStopElectronPositions(startEl, stopEl);
  occupation = *occ;

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
            }
          else
            {
              if(numRows != 0)
                {
                  WF_coeffs(i).setRows(orbital_index-numRows,j-numRows,numRows,
				       *coeffs);
                  numRows = 0;
                }
              if (orbital_index == D(0,i).dim2()) break;
            }
        }
    }
#ifdef QMC_GPU
  // I should put some more effort into this question of copy 
  // constructors/assignment
  GPUQMCBasisFunction temp(*BF, (int)(WF_coeffs(0).dim1()), Input->flags.walkers_per_pass);
  gpuBF = temp;
  gpuMatMult = GPUQMCMatrix(Input, WF_coeffs,Input->flags.walkers_per_pass);
#endif
}

void QMCSlater::allocate(int N)
{
  int ndet = WF->getNumberDeterminants();
  int nbasisfunc = WF->getNumberBasisFunctions();
  int wpp = Input->flags.walkers_per_pass;

  D.allocate(wpp,ndet);
  D_inv.allocate(wpp,ndet);
  Laplacian_D.allocate(wpp,ndet);
  Grad_D.allocate(wpp,ndet,3);

  Singular.allocate(wpp);
  Psi.allocate(wpp);
  Laplacian_PsiRatio.allocate(wpp);
  Grad_PsiRatio.allocate(wpp);

  if (Input->flags.calculate_bf_density == 1)
    Chi_Density.allocate(wpp);

  for(int j=0; j<wpp; j++)
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

#ifdef QMC_GPU
  bfData.allocate(wpp*5);
  for(int j=0; j<bfData.dim1(); j++)
    bfData(j).allocate(N,nbasisfunc);

  /* The dimensions of these data are numWalkers x numDeterminants then numElec
     x numOrb.  What I'm doing here is filling the resultsCollector database 
     with addresses of our actual data storage. This creates an alias which I
     can then send to the GPUQMC code without duplicating data storage or 
     ruining the previously existing data design
  */
  resultsCollector.allocate(ndet, wpp*5);
  for(int iDet=0; iDet<ndet; iDet++)
    {
      for(int iWalker=0; iWalker<wpp; iWalker++)
        {
          resultsCollector(iDet,iWalker*5    ) = & D(iWalker, iDet);
          resultsCollector(iDet,iWalker*5 + 1) = & Grad_D(iWalker, iDet, 0);
          resultsCollector(iDet,iWalker*5 + 2) = & Grad_D(iWalker, iDet, 1);
          resultsCollector(iDet,iWalker*5 + 3) = & Grad_D(iWalker, iDet, 2);
          resultsCollector(iDet,iWalker*5 + 4) = & Laplacian_D(iWalker, iDet);
        }
    }

#else
  Chi.allocate(N,nbasisfunc);
  Chi_laplacian.allocate(N,nbasisfunc);
  Chi_gradient.allocate(3);
  for(int j=0; j<3; j++)
    Chi_gradient(j).allocate(N,nbasisfunc);
#endif
}

QMCSlater::~QMCSlater()
{
  if (Start != Stop)
    {
#ifdef QMC_GPU
      resultsCollector.deallocate();
#endif

      for(int j=0; j < Input->flags.walkers_per_pass; j++)
	{
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

      for (int i=0; i<WF->getNumberDeterminants(); i++)
	WF_coeffs(i).deallocate();
      WF_coeffs.deallocate();

      D.deallocate();
      D_inv.deallocate();
      Laplacian_D.deallocate();
      Grad_D.deallocate();
      Singular.deallocate();
      Psi.deallocate();
      Laplacian_PsiRatio.deallocate();
      Grad_PsiRatio.deallocate();

      if (Input->flags.calculate_bf_density == 1)
	Chi_Density.deallocate();

      occupation.deallocate();

#ifdef QMC_GPU
      //gpuBF.destroy();
      gpuMatMult.destroy();
      for (int j=0; j<bfData.dim1(); j++)
	bfData(j).deallocate();
      bfData.deallocate();
#else
      Chi.deallocate();
      Chi_laplacian.deallocate();
      for (int j=0; j<3; j++)
	Chi_gradient(j).deallocate();
      Chi_gradient.deallocate();
#endif
    }
}

void QMCSlater::setStartAndStopElectronPositions(int StartEl, int StopEl)
{
  Start = StartEl;
  Stop  = StopEl;

  allocate(StopEl-StartEl+1);
}

void QMCSlater::evaluate(int num)
{
  // The new idea here is to split all the work that the evaluate function used
  // to do.  This enables QMCFunction to control when the multiplication 
  // happens relative to the inverse

  // The LU Decomposition inverse used here is O(1*N^3)
  // Updating one electron at a time is O(2*N^3-N^2)
  bool calcOK = true;
  bool printPsi = !true;

#ifdef QMC_GPU
  gpuMatMult.getResults(resultsCollector);
#endif

  for(int i=0; i<num; i++)
    {
      for (int j=0; j<WF->getNumberDeterminants(); j++)
        {
          determinant_and_inverse(D(i,j),D_inv(i,j),(Psi(i))(j),&calcOK);
          (Singular(i))(j) = !calcOK;
        }

      if( !isSingular(i) )
        {
          calculate_DerivativeRatios(i);

          if(printPsi)
            {
              printf("%4d: ",i);
              for (int j=0; j<WF->getNumberDeterminants(); j++)
                {
                  printf("psi # %s%18.15e ",(Psi(i))(j)<0?"":" ",(Psi(i))(j) );
                  printf("lap ratio # %s%18.15g\n",(Laplacian_PsiRatio(i))(j)<0?"":" ",(Laplacian_PsiRatio(i))(j) );
                }
            }
        }
    }
}

/**
This function contains the meat of the QMC calculation.
1) Basis functions are calculated for each electron position
   for each determinant:
2) A wavefunction coefficient matrix is composed of the input basis function 
   coefficients.  Some monkey business is required here because not all 
   consecutive orbitals in the input coefficients are occupied.
3) The Slater matrix is calculated along with the associated gradient and 
   laplacian matricies:
   D(numElec x numOrb) = Chi(numElec x numBasisFunction) * 
                                           WF_coeffs(numBasisFunction x numOrb)
 
  Note: the coefficient matricies are transposed.  Hand-coded matrix 
  multiplication is faster with a transposed matrix, and it enables the use of 
  memcpy for Array2D's row-major data.
*/
void QMCSlater::update_Ds(Array1D<Array2D<double>*> &X, int num)
{
  static double averageM = 0, timeM = 0;
  static double averageB = 0, timeB = 0;
  static double numT = -5;
  static int nBF = WF_coeffs(0).dim2();
  static int nOE = WF_coeffs(0).dim1();
  static const bool showTimings = false;
  static const bool printBF = false;//is there a good way to do this?
  // 1 add + 1 mul per nBF * nOE * nOE
  double numOps = 5*num*(nBF * nOE * nOE * 2.0 - nOE * nOE) / 1000.0;

#ifdef QMC_GPU
  GLuint texture;
#endif
  Stopwatch sw = Stopwatch();
  sw.reset();

#ifdef QMC_GPU
  static int mat_multiplier = gpuMatMult.getNumIterations();
  static int bas_multiplier = gpuBF.getNumIterations();

  if(showTimings)
    {
      sw.reset(); sw.start();
    }
  texture = gpuBF.runCalculation(X,num, Start, Stop);
  if(showTimings)
    {
      sw.stop();
      timeB = sw.timeMS();
      if(numT >= 0) averageB += sw.timeMS();
    }

  if(showTimings)
    {
      sw.reset(); sw.start();
    }
  gpuMatMult.runCalculation(num,texture);
  if(showTimings)
    {
      sw.stop();
      timeM = sw.timeMS();
      if(numT >= 0) averageM += sw.timeMS();
    }

#else
  static int mat_multiplier = 1;
  static int bas_multiplier = 1;
  timeB = 0; timeM = 0;
  for(int walker = 0; walker < num; walker++)
    {
      if(showTimings)
        {
          sw.reset(); sw.start();
        }
      BF->evaluateBasisFunctions(*X(walker),Start,Stop,
                                 Chi,
                                 Chi_gradient(0),
                                 Chi_gradient(1),
                                 Chi_gradient(2),
                                 Chi_laplacian);
      if(showTimings)
        {
          sw.stop();
          timeB += sw.timeMS();
          if(numT >= 0) averageB += sw.timeMS();
        }

      if(num > 1 && false)
        {
          cout << "\npsi_bf" << walker << ":\n" << Chi;
          if(!false)
            {
              cout << "grx_bf" << walker << ":\n" << Chi_gradient(0);
              cout << "gry_bf" << walker << ":\n" << Chi_gradient(1);
              cout << "grz_bf" << walker << ":\n" << Chi_gradient(2);
              cout << "lap_bf" << walker << ":\n" << Chi_laplacian;
            }
        }

      if(showTimings)
        {
          sw.reset(); sw.start();
        }
      for(int i=0; i<WF->getNumberDeterminants(); i++)
        {
          D(walker,i)           = Chi * WF_coeffs(i);
          Laplacian_D(walker,i) = Chi_laplacian * WF_coeffs(i);
          Grad_D(walker,i,0)    = Chi_gradient(0) * WF_coeffs(i);
          Grad_D(walker,i,1)    = Chi_gradient(1) * WF_coeffs(i);
          Grad_D(walker,i,2)    = Chi_gradient(2) * WF_coeffs(i);
        }
      if(showTimings)
        {
          sw.stop();
          timeM += sw.timeMS();
          if(numT >= 0) averageM += sw.timeMS();
        }

      if(num >= 1 && false)
        {
          //cout << "\npsi_mm" << walker << ":\n" << D(walker,0);
          //cout << "grx_mm" << walker << ":\n" << Grad_D(walker,0,0);
          //cout << "gry_mm" << walker << ":\n" << Grad_D(walker,0,1);
          //cout << "grz_mm" << walker << ":\n" << Grad_D(walker,0,2);
          cout << "lap_mm" << walker << ":\n" << Laplacian_D(walker,0);
        }

    }

  for(int walker = 0; walker < num; walker++)
    {
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
#endif

  if(printBF)
    {/*
       cout << i << ": ";
       for (int j=0; j<WF->getNumberDeterminants(); j++){
       printf("psi %s%8.5e ",(Psi(i))(j)<0?"":" ",(Psi(i))(j) );
       printf("lap ratio %s%8.5g\n",(Psi(i))(j)<0?"":" ",(Laplacian_PsiRatio(i))(j) );
       }
     */
    }

  if(showTimings)
    {
      if( num > 1) numT++;
      cout << "bf: " << (int)(timeB/bas_multiplier+0.5) << " (" << (int)(averageB/(numT*bas_multiplier)+0.5) << ") ";
      cout << "mm: " << (int)(timeM+0.5) << " (" << (int)(averageM/numT+0.5) << ") ";
      cout << "; mflops: " << (int)(mat_multiplier*numOps/timeM+0.5) <<
      " (" << (int)(mat_multiplier*numOps/(averageM/numT)+0.5) << ")\n";
    }
}

/**
This function calculates the derivative ratios. that is, del psi/psi, and lap 
psi/psi.  D_inv is the transpose of the true inverse, but it turns out that if 
we use the transpose here, then the calculation of the ratios turns into a sort
of dot product, which not only does ATLAS know how to do it, it makes the code 
look neater. the hand-coded dot product probably doesn't save much time.
At this time, D_inv and Grad_D and Laplacian_D are all qmcfloat type. The 
explicit typecast when creating the final result (double) should emphasize 
this.
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
