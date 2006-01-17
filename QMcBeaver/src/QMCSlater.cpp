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

static const bool showTimings = false;
static const bool printPsi    = false;

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
              if (orbital_index == Stop-Start+1) break;
            }
        }
    }
#ifdef QMC_GPU
  // I should put some more effort into this question of copy 
  // constructors/assignment
  GPUQMCBasisFunction temp(*BF, (int)(WF_coeffs(0).dim1()), Input->flags.getNumGPUWalkers());
  gpuBF = temp;
  gpuMatMult = GPUQMCMatrix(Input, WF_coeffs,Input->flags.getNumGPUWalkers());
#endif
}

void QMCSlater::allocate(int N)
{
  int ndet = WF->getNumberDeterminants();
  int nbasisfunc = WF->getNumberBasisFunctions();
  int wpp = Input->flags.walkers_per_pass;
  int gpp = Input->flags.getNumGPUWalkers();


  D.allocate(wpp-gpp,ndet);
  D_inv.allocate(wpp-gpp,ndet);
  Laplacian_D.allocate(wpp-gpp,ndet);
  Grad_D.allocate(wpp-gpp,ndet,3);
  for(int j=0; j<D.dim1(); j++)
    {
      for (int i=0; i<D.dim2(); i++)
        {
          D(j,i).allocate(N,N);
          D_inv(j,i).allocate(N,N);
          Laplacian_D(j,i).allocate(N,N);
          for(int k=0; k<3; k++)
            Grad_D(j,i,k).allocate(N,N);
        }
    }
#ifdef QMC_GPU
  gpu_D.allocate(gpp,ndet);
  gpu_D_inv.allocate(gpp,ndet);
  gpu_Laplacian_D.allocate(gpp,ndet);
  gpu_Grad_D.allocate(gpp,ndet,3);
  for(int j=0; j<gpu_D.dim1(); j++)
    {
      for (int i=0; i<gpu_D.dim2(); i++)
        {
          gpu_D(j,i).allocate(N,N);
          gpu_D_inv(j,i).allocate(N,N);
          gpu_Laplacian_D(j,i).allocate(N,N);
          for(int k=0; k<3; k++)
            gpu_Grad_D(j,i,k).allocate(N,N);
        }
    }
#endif

  Singular.allocate(wpp);
  Psi.allocate(wpp);
  Laplacian_PsiRatio.allocate(wpp);
  Grad_PsiRatio.allocate(wpp);

  if (Input->flags.calculate_bf_density == 1)
    Chi_Density.allocate(wpp);

  for(int j=0; j<wpp; j++)
    {
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
  /* The dimensions of these data are numWalkers x numDeterminants then numElec
     x numOrb.  What I'm doing here is filling the resultsCollector database 
     with addresses of our actual data storage. This creates an alias which I
     can then send to the GPUQMC code without duplicating data storage or 
     ruining the previously existing data design.

     resultsCollector will be assigned the first gpu_walkers_per_pass elements.
  */
  resultsCollector.allocate(ndet, gpp*5);
  for(int iDet=0; iDet<ndet; iDet++)
    {
      for(int iWalker=0; iWalker<gpp; iWalker++)
        {
          resultsCollector(iDet,iWalker*5    ) = & gpu_D(iWalker, iDet);
          resultsCollector(iDet,iWalker*5 + 1) = & gpu_Grad_D(iWalker, iDet, 0);
          resultsCollector(iDet,iWalker*5 + 2) = & gpu_Grad_D(iWalker, iDet, 1);
          resultsCollector(iDet,iWalker*5 + 3) = & gpu_Grad_D(iWalker, iDet, 2);
          resultsCollector(iDet,iWalker*5 + 4) = & gpu_Laplacian_D(iWalker, iDet);
        }
    }
#endif

  Chi.allocate(N,nbasisfunc);
  Chi_laplacian.allocate(N,nbasisfunc);
  Chi_gradient.allocate(3);
  for(int j=0; j<3; j++)
    Chi_gradient(j).allocate(N,nbasisfunc);
}

QMCSlater::~QMCSlater()
{
  if (Start != Stop)
  {
    for(int j=0; j<D.dim1(); j++)
    {
      for (int i=0; i<D.dim2(); i++)
      {
        D(j,i).deallocate();
        D_inv(j,i).deallocate();
        Laplacian_D(j,i).deallocate();
        for(int k=0; k<3; k++)
          Grad_D(j,i,k).deallocate();
      }
    }
    D.deallocate();
    D_inv.deallocate();
    Laplacian_D.deallocate();
    Grad_D.deallocate();

#ifdef QMC_GPU
    for(int j=0; j<gpu_D.dim1(); j++)
    {
      for (int i=0; i<gpu_D.dim2(); i++)
      {
        gpu_D(j,i).deallocate();
        gpu_D_inv(j,i).deallocate();
        gpu_Laplacian_D(j,i).deallocate();
        for(int k=0; k<3; k++)
          gpu_Grad_D(j,i,k).deallocate();
      }
    }
    gpu_D.deallocate();
    gpu_D_inv.deallocate();
    gpu_Laplacian_D.deallocate();
    gpu_Grad_D.deallocate();
#endif

    for(int j=0; j < Input->flags.walkers_per_pass; j++)
    {
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

    Singular.deallocate();
    Psi.deallocate();
    Laplacian_PsiRatio.deallocate();
    Grad_PsiRatio.deallocate();

    if (Input->flags.calculate_bf_density == 1)
      Chi_Density.deallocate();

    occupation.deallocate();

#ifdef QMC_GPU
    //We do not need to deallocate the contents of
    //resultsCollector because the contents were
    //references to gpu_D, gpu_Laplacian, etc
    //and those were already deallocated.
    resultsCollector.deallocate();

    gpuMatMult.destroy();
#endif

    Chi.deallocate();
    Chi_laplacian.deallocate();
    for (int j=0; j<3; j++)
      Chi_gradient(j).deallocate();
    Chi_gradient.deallocate();
  }
}

void QMCSlater::setStartAndStopElectronPositions(int StartEl, int StopEl)
{
  Start = StartEl;
  Stop  = StopEl;

  allocate(StopEl-StartEl+1);
}

void QMCSlater::update_Ds()
{
  // The new idea here is to split all the work that the evaluate function used
  // to do.  This enables QMCFunction to control when the multiplication 
  // happens relative to the inverse

#ifdef QMC_GPU
  /* This will fill the resultsCollector with
  the results from the GPU. However, all of 
  resultsCollector's elements are really
  just pointers to the D object.*/
  gpuMatMult.getResults(resultsCollector);
  processInverse(0,gpu_D,gpu_D_inv,gpu_Laplacian_D,gpu_Grad_D);
#endif
  int gpp = Input->flags.getNumGPUWalkers();
  processInverse(gpp,D,D_inv,Laplacian_D,Grad_D);
}

template <class T>
void QMCSlater::processInverse(int start,
      Array2D< Array2D<T> > & psi, Array2D< Array2D<T> > & inv, 
      Array2D< Array2D<T> > & lap, Array3D< Array2D<T> > & grad)
{
  // The LU Decomposition inverse used here is O(1*N^3)
  // Updating one electron at a time is O(2*N^3-N^2)
  bool calcOK = true;
  for(int i=0; i<psi.dim1(); i++)
    {
      for (int j=0; j<WF->getNumberDeterminants(); j++)
        {
          psi(i,j).determinant_and_inverse(inv(i,j),(Psi(i+start))(j),&calcOK);
          (Singular(i+start))(j) = !calcOK;
        }

      if( !isSingular(i+start) )
        {
          calculate_DerivativeRatios(i,start,inv,lap,grad);

          if(printPsi)
            {
              printf("%4d: ",i+start);
              for (int j=0; j<WF->getNumberDeterminants(); j++)
                {
                  printf("psi # %s%18.15e ",(Psi(i))(j)<0?"":" ",(Psi(i+start))(j) );
                  printf("lap ratio # %s%18.15g\n",(Laplacian_PsiRatio(i+start))(j)<0?"":" ",(Laplacian_PsiRatio(i+start))(j) );
                }
            }
        }
    }
}

template void QMCSlater::processInverse<float>
(int start, Array2D< Array2D<float> > & psi, Array2D< Array2D<float> > & inv, 
            Array2D< Array2D<float> > & lap, Array3D< Array2D<float> > & grad);

template void QMCSlater::processInverse<double>
(int start, Array2D< Array2D<double> > & psi, Array2D< Array2D<double> > & inv, 
            Array2D< Array2D<double> > & lap, Array3D< Array2D<double> > & grad);

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
#ifdef QMC_GPU
void QMCSlater::gpuEvaluate(Array1D<Array2D<double>*> &X, int num)
{
  static double averageM = 0, timeM = 0;
  static double averageB = 0, timeB = 0;
  static double numT = -5;
  static int nBF = WF_coeffs(0).dim2();
  static int nOE = WF_coeffs(0).dim1();
  static int mat_multiplier = gpuMatMult.getNumIterations();
  static int bas_multiplier = gpuBF.getNumIterations();
  // 1 add + 1 mul per nBF * nOE * nOE
  double numOps = 5*num*(nBF * nOE * nOE * 2.0 - nOE * nOE) / 1000.0;

  GLuint texture;
  Stopwatch sw = Stopwatch();
  sw.reset();

  if(showTimings){sw.reset(); sw.start();}

  texture = gpuBF.runCalculation(X,num, Start, Stop);
  
  if(showTimings)
    {
      glFinish();
      sw.stop();
      timeB = sw.timeMS();
      if(numT >= 0) averageB += sw.timeMS();
    }

  if(showTimings){sw.reset(); sw.start();}

  gpuMatMult.runCalculation(num,texture);
  
  if(showTimings)
    {
      gpuMatMult.getResults(resultsCollector);
      glFinish();
      sw.stop();
      timeM = sw.timeMS();
      if(numT >= 0) averageM += sw.timeMS();
      if( num > 1) numT++;
      cout << "\ngpu bf: " << (int)(timeB/bas_multiplier+0.5) << " ( " << (int)(averageB/(numT*bas_multiplier)+0.5) << ") ";
      cout << "mm: " << (int)(timeM/mat_multiplier+0.5) << " (" << (int)(averageM/(numT*mat_multiplier)+0.5) << ") ";
      cout << "; mflops: " << (int)(mat_multiplier*numOps/timeM+0.5) <<
      " (" << (int)(mat_multiplier*numOps/(averageM/numT)+0.5) << ")\n";
    }
}
#endif

void QMCSlater::evaluate(Array1D<Array2D<double>*> &X, int num, int start)
{
  static double averageM = 0, timeM = 0;
  static double averageB = 0, timeB = 0;
  static double numT = -5;
  static int nBF = WF_coeffs(0).dim2();
  static int nOE = WF_coeffs(0).dim1();
  static int mat_multiplier = 1000;
  static int bas_multiplier = 1000;
  // 1 add + 1 mul per nBF * nOE * nOE
  double numOps = 5*num*(nBF * nOE * nOE * 2.0 - nOE * nOE) / 1000.0;

  Stopwatch sw = Stopwatch();
  sw.reset();

  timeB = 0; timeM = 0;
  for(int walker = 0; walker < num; walker++)
    {
      if(showTimings)
        {
          sw.reset(); sw.start();
        }
      BF->evaluateBasisFunctions(*X(walker+start),Start,Stop,
                                 Chi,
                                 Chi_gradient(0),
                                 Chi_gradient(1),
                                 Chi_gradient(2),
                                 Chi_laplacian);
      if(showTimings)
        {
          sw.stop();
          timeB += sw.timeUS();
          if(numT >= 0) averageB += sw.timeUS();
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
          timeM += sw.timeUS();
          //cout << "swtime " << sw.timeUS() << endl;
          if(numT >= 0) averageM += sw.timeUS();
        }
    }

  for(int walker = 0; walker < num; walker++)
    {
      if (Input->flags.calculate_bf_density == 1)
        {
          Chi_Density(walker+start) = 0.0;
          for (int i=0; i<WF->getNumberBasisFunctions(); i++)
            for (int j=0; j<D(0,0).dim1(); j++)
              {
                (Chi_Density(walker))(i) += Chi(j,i);
              }
        }
    }

  if(showTimings)
    {
      numT++;
      cout << "cpu bf: " << (int)(timeB/bas_multiplier+0.5) << " ( " << (int)(averageB/(numT*bas_multiplier)+0.5) << ") ";
      cout << "mm: " << (int)(timeM/mat_multiplier+0.5) << " (" << (int)(averageM/(numT*mat_multiplier)+0.5) << ") ";
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
template <class T>
void QMCSlater::calculate_DerivativeRatios(int k, int start, Array2D< Array2D<T> > & inv, 
      Array2D< Array2D<T> > & lap, Array3D< Array2D<T> > & grad)
{
  int numElectrons = Stop-Start+1;
  for(int i=0; i<WF->getNumberDeterminants(); i++)
    {
      double** grad_psiratioArray = Grad_PsiRatio(k+start).array()[i];

      (Laplacian_PsiRatio(k+start))(i) = (double)((lap(k,i)).dotAllElectrons(inv(k,i)));

      for(int j=0; j<numElectrons; j++)
        {
          grad_psiratioArray[j][0] = (double)((grad(k,i,0)).dotOneElectron(inv(k,i),j));
          grad_psiratioArray[j][1] = (double)((grad(k,i,1)).dotOneElectron(inv(k,i),j));
          grad_psiratioArray[j][2] = (double)((grad(k,i,2)).dotOneElectron(inv(k,i),j));
        }
    }
}

template void QMCSlater::calculate_DerivativeRatios<float>
(int k, int start, Array2D< Array2D<float> > & inv, 
 Array2D< Array2D<float> > & lap, Array3D< Array2D<float> > & grad);

template void QMCSlater::calculate_DerivativeRatios<double>
(int k, int start, Array2D< Array2D<double> > & inv, 
 Array2D< Array2D<double> > & lap, Array3D< Array2D<double> > & grad);

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
