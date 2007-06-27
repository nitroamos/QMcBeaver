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
    if (Input->flags.replace_electron_nucleus_cusps == 1)
        ElectronNucleusCusp = rhs.ElectronNucleusCusp;

    Singular = rhs.Singular;

    // This was done in this round about way to avoid some weird compiler errors
    // dealing with when const was used.  In the simplest format this would be
    // allocate( rhs.D.dim1() );
    D = rhs.D;
    allocate();
}

void QMCSlater::initialize(QMCInput *INPUT,
			   int startEl, int stopEl,
			   bool selectAlpha)
{
  Start = startEl;
  Stop = stopEl;
  isAlpha = selectAlpha;
  Input = INPUT;
  BF = &Input->BF;
  WF = &Input->WF;

  if(getNumberElectrons() == 0)
    return;
  
  allocate();
  
  for(int i=0; i<Singular.dim1(); i++)
    Singular(i) = false;
  
  if (Input->flags.replace_electron_nucleus_cusps == 1)
    for (int i=0; i<WF->getNumberDeterminants(); i++)
      {
	ElectronNucleusCusp(i).initialize(Input,startEl,stopEl,
					  * WF->getCoeff(i,isAlpha));
	ElectronNucleusCusp(i).fitReplacementOrbitals();
      }
  
#ifdef QMC_GPU
  // I should put some more effort into this question of copy
  // constructors/assignment
  GPUQMCBasisFunction temp(*BF, getNumberElectrons(), Input->flags.getNumGPUWalkers());
  gpuBF = temp;
  // I replaced WF_coeffs with getCoeff(ci,isAlpha), so GPU compiles will fail right here
  gpuMatMult = GPUQMCMatrix(Input, WF_coeffs,Input->flags.getNumGPUWalkers());
#endif
}

void QMCSlater::allocate()
{
  int N = getNumberElectrons();
    int ndet = WF->getNumberDeterminants();
    int nbasisfunc = WF->getNumberBasisFunctions();
    int wpp = Input->flags.walkers_per_pass;
    int gpp = Input->flags.getNumGPUWalkers();


    if(wpp - gpp > 0)
    {
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

    if (Input->flags.optimize_Orbitals == 1)
      {
	p_a.allocate(wpp);
	p2_xa.allocate(wpp);
	p3_xxa.allocate(wpp);
	
	for(int j=0; j<wpp; j++)
	  {
	    p_a(j).allocate(ndet);
	    p2_xa(j).allocate(ndet,N,3);
	    p3_xxa(j).allocate(ndet);

	    for(int pi=0; pi<ndet; pi++)
	      {
		(p_a(j))(pi).allocate(N,nbasisfunc);
		(p3_xxa(j))(pi).allocate(N,nbasisfunc);

		for(int pj=0; pj<N; pj++)
		  for(int pk=0; pk<3; pk++)
		    (p2_xa(j))(pi,pj,pk).allocate(N,nbasisfunc);		    
	      }
	  }
      }

    if (Input->flags.replace_electron_nucleus_cusps == 1)
        ElectronNucleusCusp.allocate(WF->getNumberDeterminants());

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

    if(Input->flags.optimize_Orbitals)
      {
	/*
	  If we're optimizing the orbitals, then we will
	  need to use this data to calculate the parameter
	  derivatives.
	*/
	Chi.allocate(wpp);
	Chi_laplacian.allocate(wpp);
	Chi_gradient.allocate(wpp);
      } else {
	/*
	  We don't need this data after we've multiplied these
	  by the WF_Coeffs
	*/
	Chi.allocate(1);
	Chi_laplacian.allocate(1);
	Chi_gradient.allocate(1);
      }

    for(int k=0; k<Chi.dim1(); k++)
      {
	Chi(k).allocate(N,nbasisfunc);
	Chi_laplacian(k).allocate(N,nbasisfunc);
	Chi_gradient(k).allocate(3);
	for(int j=0; j<3; j++)
	  (Chi_gradient(k))(j).allocate(N,nbasisfunc);
      }
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

	if (Input->flags.optimize_Orbitals == 1)
	  {
	    for(int j=0; j<p_a.dim1(); j++)
	      {
		for(int pi=0; pi<p_a(j).dim1(); pi++)
		  {
		    (p_a(j))(pi).deallocate();
		    (p3_xxa(j))(pi).deallocate();
		    
		    for(int pj=0; pj<p2_xa(j).dim2(); pj++)
		      for(int pk=0; pk<3; pk++)
			(p2_xa(j))(pi,pj,pk).deallocate();
		  }
		
		p_a(j).deallocate();
		p2_xa(j).deallocate();
		p3_xxa(j).deallocate();
	      }
	    p_a.deallocate();
	    p2_xa.deallocate();
	    p3_xxa.deallocate();
	  }
	
        Singular.deallocate();
        Psi.deallocate();
        Laplacian_PsiRatio.deallocate();
        Grad_PsiRatio.deallocate();

        if (Input->flags.calculate_bf_density == 1)
            Chi_Density.deallocate();

        if (Input->flags.replace_electron_nucleus_cusps == 1)
            ElectronNucleusCusp.deallocate();

#ifdef QMC_GPU
        //We do not need to deallocate the contents of
        //resultsCollector because the contents were
        //references to gpu_D, gpu_Laplacian, etc
        //and those were already deallocated.
        resultsCollector.deallocate();

        gpuMatMult.destroy();
#endif

	for(int k=0; k<Chi.dim1(); k++)
	  {
	    Chi(k).deallocate();
	    Chi_laplacian(k).deallocate();
	    for (int j=0; j<3; j++)
	      (Chi_gradient(k))(j).deallocate();
	    Chi_gradient(k).deallocate();
	  }
	Chi.deallocate();
	Chi_laplacian.deallocate();
	Chi_gradient.deallocate();
    }
}

void QMCSlater::update_Ds()
{
  if(getNumberElectrons() == 0)
    return;

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

    if(Input->flags.optimize_Orbitals)
      {
	/*
	  Let d |D| be shorthand for  \frac{ d |D| }{ d c_{jk} }
	  Let nabla |D| refer to derivatives with respect to position

	  Looking at the code probably doesn't make sense.

	  Jacobi's formula:
	  d|D| = |D| tr( D^{-1} . d D )

	  After figuring out how to take these derivatives, I took
	  shortcuts to calculate derivatives wrt all the c_{jk} at
	  once. Specifically,

	  Single derivatives:
	  d |D| = |D| D^{-1} . Chi

	  Mixed derivatives:
	  d \frac{ \nabla |D| }{ |D| }
	  = d tr( D^{-1} . \nabla D )
	  = tr( D^{-1} . d \nabla D ) + tr( d D^{-1} . \nabla D )
	  = tr( D^{-1} . d \nabla D ) - tr( D^{-1} . d D . D^{-1} . \nabla D )
	  = tr( D^{-1} . d \nabla D ) - tr( \nabla D . D^{-1} . d D . D^{-1} )
	*/
		
	// \nabla D . D^{-1} . d D
	Array1D< Array2D<double> > nabDinvDdD;
	nabDinvDdD.allocate(3);
	
	// \nabla^2 D . D^{-1} . d D
	Array2D<double> nab2DinvDdD;

	Array2D<double> temp;

	for(int w=0; w<D.dim1(); w++)
	  {
	    for(int ci=0; ci<WF->getNumberDeterminants(); ci++)
	      {
		Array2D<double> * pa = & (p_a(w+gpp))(ci);
		Array2D<double> * p3xxa = & (p3_xxa(w+gpp))(ci);
		* pa    = 0.0;
		* p3xxa = 0.0;
		
		// The data was stored transposed
		D_inv(w,ci).transpose();
		
		//calculate  d |D|
		D_inv(w,ci).gemm(Chi(w),*pa,false);
		(*pa) *= (Psi(w+gpp))(ci);
		
		for(int xyz=0; xyz<3; xyz++)
		  {
		    Grad_D(w,ci,xyz).gemm(D_inv(w,ci),temp,false);
		    temp.gemm(Chi(w),nabDinvDdD(xyz),false);
		    //temp = 0.0;
		  }
		Laplacian_D(w,ci).gemm(D_inv(w,ci),temp,false);
		temp.gemm(Chi(w),nab2DinvDdD,false);
		
		int nor = pa->dim1();
		int nbf = pa->dim2();
		for(int e=0; e<nor; e++)
		  {
		    for(int o=0; o<nor; o++)
		      {
			for(int bf=0; bf<nbf; bf++)
			  {
			    for(int xyz=0; xyz<3; xyz++)
			      {
				double d2xa = 0.0;
				d2xa += (D_inv(w,ci))(o,e) * ((Chi_gradient(w))(xyz))(e,bf);
				d2xa -= (D_inv(w,ci))(o,e) * (nabDinvDdD(xyz))(e,bf);
				((p2_xa(w+gpp))(ci,e,xyz))(o,bf) = d2xa;
			      }
			    double d3xxa = 0.0;
			    d3xxa += (D_inv(w,ci))(o,e) * (Chi_laplacian(w))(e,bf);
			    d3xxa -= (D_inv(w,ci))(o,e) * nab2DinvDdD(e,bf);
			    (*p3xxa)(o,bf) += d3xxa;			    
			  }
		      }
		  }
	      }
	  }
	
	for(int dim=0; dim<nabDinvDdD.dim1(); dim++)
	  nabDinvDdD(dim).deallocate();
	nabDinvDdD.deallocate();
	nab2DinvDdD.deallocate();
	temp.deallocate();
      }
}

template <class T>
void QMCSlater::processInverse(int start,
                               Array2D< Array2D<T> > & psi, Array2D< Array2D<T> > & inv,
                               Array2D< Array2D<T> > & lap, Array3D< Array2D<T> > & grad)
{
    // The LU Decomposition inverse used here is O(1*N^3)
    // Updating one electron at a time is O(2*N^3-N^2)

    /* SVD is a more careful way of calculating the inverse of a matrix.
       However, it seems like it usually doesn't make much of a difference.
       Lastly, it's guaranteed to be slower. */
    const static bool useSVD       = false;
    /* This will print out an error if the matrix inversion completely failed.*/
    const static bool checkInverse = false;

    bool calcOK = true;

    static Array2D<T> psiCopy;
    static Array2D<T> identity;
    static Array1D<T> W;
    static Array2D<T> r;

    if(useSVD){
        W.allocate(psi(0,0).dim1());
        r.allocate(psi(0,0).dim1(),psi(0,0).dim2());
    }

    for(int i=0; i<psi.dim1(); i++)
    {
        for (int j=0; j<WF->getNumberDeterminants(); j++)
        {

            if(checkInverse)
                psiCopy = psi(i,j);

            if(useSVD){
#if ! defined( USING_QSC )
                calcOK = SVDecompose(psi(i,j), W, r, 10) == 0 ? true : false;

                (Psi(i+start))(j) = 1.0;
                double smallest = W(0);
                double largest = W(0);
                for(int ii=0; ii<W.dim1(); ii++){
                    //NOTE: SVD will only provide the absolute value of the determinant
                    (Psi(i+start))(j) *= W(ii);
                    smallest = min(fabs((double)W(ii)),smallest);
                    largest = max(fabs((double)W(ii)),largest);

                    //the W values are measures of instability
                    //you *might* get away with setting the small ones
                    //to zero before back substitution
                    //printf("%20.10e ",W(ii));
                }
                //printf("\n");
                //cout << "Smallest = " << smallest << " largest = " << largest << " Condition number = " << largest/smallest << endl;

                SVDFwBackSubst(psi(i,j),W,r,inv(i,j));
#endif
            } else {

                //This choice will automatically use LAPACK if it has been linked.
                psi(i,j).determinant_and_inverse(inv(i,j),(Psi(i+start))(j),&calcOK);

                //This might improve the precision of the result
                //psi(i,j).mprove(psiCopy,inv(i,j));
            }

            if(checkInverse){
                psiCopy.gemm(inv(i,j),identity,!false);
                double error = 0;
                for(int ii=0; ii<identity.dim1(); ii++)
                    for(int jj=0; jj<identity.dim1(); jj++)
                        error += identity(ii,jj) * identity(ii,jj);
                error -= identity.dim1();
                if(abs(error) > 1e-10)
                    cerr << "ERROR: matrix inversion error of " << error << endl;
            }

            (Singular(i+start))(j) = !calcOK;
        }

        if( !isSingular(i+start) )
        {
            calculate_DerivativeRatios(i,start,inv,lap,grad);

            if(false)
            {
                double psiNorm = 0;
                double invNorm = 0;

                printf("%4d: ",i+start);
                for (int j=0; j<WF->getNumberDeterminants(); j++)
                {
                    invNorm = inv(i,j).pInfNorm();
                    printf("psi %20.10e ",  (Psi(i+start))(j) );
                    printf("lap %20.10e ", (Laplacian_PsiRatio(i+start))(j) );
                    printf("Pfn %20.10e ", psiNorm );
                    printf("Ifn %20.10e ", invNorm );
                    printf("CN  %20.10e ", psiNorm*invNorm);
                    printf("ICN %20.10e ", identity.pInfNorm());
                    printf("\n");
                    /*
                      for(int ii=0; ii<inv(i,j).dim1(); ii++){
                      for(int ij=0; ij<inv(i,j).dim2(); ij++){
                      printf("%18.10e",(inv(i,j))(ii,ij));
                      }
                      printf("\n");
                    }
                    */
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
  if(getNumberElectrons == 0)
    return;

    static double averageM = 0, timeM = 0;
    static double averageB = 0, timeB = 0;
    static double numT = -5;
    static int nBF = WF->getNumberBasisFunctions();
    static int nOE = getNumberElectrons();
    static int mat_multiplier = gpuMatMult.getNumIterations() * 1000;
    static int bas_multiplier = gpuBF.getNumIterations();
    // 1 add + 1 mul per nBF * nOE * nOE
    double numOps = 5*num*(nBF * nOE * nOE * 2.0 - nOE * nOE) / 1000.0;

    GLuint texture;
    Stopwatch sw = Stopwatch();
    sw.reset();

    if(showTimings)
    {
        sw.reset(); sw.start();
    }

    texture = gpuBF.runCalculation(X,num, Start, Stop);

    if(showTimings)
    {
        glFinish();
        sw.stop();
        timeB = sw.timeUS();
        if(numT >= 0) averageB += sw.timeUS();
    }

    if(showTimings)
    {
        sw.reset(); sw.start();
    }

    gpuMatMult.runCalculation(num,texture);

    if(showTimings)
    {
        //the more standard timing test does not include the download time
        //gpuMatMult.getResults(resultsCollector);
        glFinish();
        sw.stop();
        timeM = sw.timeUS();
        if(numT >= 0) averageM += sw.timeUS();

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
  if(getNumberElectrons() == 0)
    return;

    static double averageM = 0, timeM = 0;
    static double averageB = 0, timeB = 0;
    static double numT = -5;
    static int nBF = WF->getNumberBasisFunctions();
    static int nOE = getNumberElectrons();
    static int mat_multiplier = 1000;
    static int bas_multiplier = 1000;
    // 1 add + 1 mul per nBF * nOE * nOE
    double numOps = 5*num*(nBF * nOE * nOE * 2.0 - nOE * nOE) / 1000.0;
    int aveB=0, aveM=0, aveF=0;

    Stopwatch sw = Stopwatch();
    sw.reset();

    timeB = 0; timeM = 0;
    for(int walker = 0; walker < num; walker++)
    {
        if(showTimings)
        {
            sw.reset(); sw.start();
        }

	int Chi_i = 0;
	if(globalInput.flags.optimize_Orbitals)
	  {
	    /*
	      Some runs we might not be optimizing orbitals,
	      others we might. This is how we're going to make
	      sure that we have allocated the memory.
	    */
	    if(p_a.dim1() != Input->flags.walkers_per_pass)
	      allocate();
	    Chi_i = walker;	    
	  }

        BF->evaluateBasisFunctions(*X(walker+start),Start,Stop,
                                   Chi(Chi_i),
                                   (Chi_gradient(Chi_i))(0),
                                   (Chi_gradient(Chi_i))(1),
                                   (Chi_gradient(Chi_i))(2),
                                   Chi_laplacian(Chi_i));

        if(showTimings)
        {
	  sw.stop();
	  timeB += sw.timeUS();
	  if(numT >= 0) averageB += sw.timeUS();
	  sw.reset(); sw.start();
        }
        for(int i=0; i<WF->getNumberDeterminants(); i++)
	  {
	    Chi(Chi_i).gemm( * WF->getCoeff(i,isAlpha),
			     D(walker,i),true);
            Chi_laplacian(Chi_i).gemm( * WF->getCoeff(i,isAlpha),
				       Laplacian_D(walker,i),true);
            (Chi_gradient(Chi_i))(0).gemm( * WF->getCoeff(i,isAlpha),
					   Grad_D(walker,i,0),true);
            (Chi_gradient(Chi_i))(1).gemm( * WF->getCoeff(i,isAlpha),
					   Grad_D(walker,i,1),true);
            (Chi_gradient(Chi_i))(2).gemm( * WF->getCoeff(i,isAlpha),
					  Grad_D(walker,i,2),true);
	    
            if (Input->flags.replace_electron_nucleus_cusps == 1)
                ElectronNucleusCusp(i).replaceCusps(*X(walker+start),
						    D(walker,i),
						    Grad_D(walker,i,0),
						    Grad_D(walker,i,1),
						    Grad_D(walker,i,2),
						    Laplacian_D(walker,i));
        }
        if(showTimings)
        {
            sw.stop();
            timeM += sw.timeUS();
            if(numT >= 0) averageM += sw.timeUS();
        }
    }

    for(int walker = 0; walker < num; walker++)
      {
	int Chi_i = 0;
	if(Input->flags.optimize_Orbitals)
	  Chi_i = walker;
	
	if (Input->flags.calculate_bf_density == 1)
	  {
	    Chi_Density(walker+start) = 0.0;
	    for (int i=0; i<WF->getNumberBasisFunctions(); i++)
	      for (int j=0; j<D(0,0).dim1(); j++)
                {
		  (Chi_Density(walker))(i) += (Chi(Chi_i))(j,i);
                }
	  }
      }
    
    if(showTimings)
    {
        //if you don't check this, your averages might be off
        //when switching from equilibration calls to production calls
        if(num>1) numT++;
        if(numT > 0)
        {
            aveB = (int)(averageB/(numT*bas_multiplier)+0.5);
            aveM = (int)(averageM/(numT*mat_multiplier)+0.5);
            aveF = (int)(numOps/(averageM/(numT*mat_multiplier))+0.5);
        }
        else
        {
            aveB = 0;
            aveM = 0;
            aveF = 0;
        }

        cout << "cpu bf: " << (int)(timeB/bas_multiplier+0.5) << " (" << aveB << ") ";
        cout << "mm: "     << (int)(timeM/mat_multiplier+0.5) << " (" << aveM << ") ";
        cout << "; mflops: " << (int)(mat_multiplier*numOps/timeM+0.5) <<
        " (" << aveF << ")\n";
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
void QMCSlater::calculate_DerivativeRatios(int k, int start,
					   Array2D< Array2D<T> > & inv,
					   Array2D< Array2D<T> > & lap, 
					   Array3D< Array2D<T> > & grad)
{
  int numElectrons = getNumberElectrons();
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
  if(getNumberElectrons() != 0)
    return &Psi(i);

  if(Psi.dim1() == 0)
    {
      Psi.allocate(1);
      Psi(0).allocate(WF->getNumberDeterminants());
      Psi(0) = 1.0;
    }
  return &Psi(0);
}

Array1D<double>* QMCSlater::getLaplacianPsiRatio(int i)
{
  if(getNumberElectrons() != 0)
    return &Laplacian_PsiRatio(i);

  if(Laplacian_PsiRatio.dim1() == 0)
    {
      Laplacian_PsiRatio.allocate(1);
      Laplacian_PsiRatio(0).allocate(WF->getNumberDeterminants());
      Laplacian_PsiRatio(0) = 0.0;
    }
  return &Laplacian_PsiRatio(0);  
}

Array3D<double>* QMCSlater::getGradPsiRatio(int i)
{
  if(getNumberElectrons() != 0)
    return &Grad_PsiRatio(i);

  return 0;
}

Array1D<double>* QMCSlater::getChiDensity(int i)
{
  if(getNumberElectrons() != 0)
    return &Chi_Density(i);
  
  return 0;
}

bool QMCSlater::isSingular(int j)
{
  if(getNumberElectrons() == 0)
    return false;

  for (int i=0; i<WF->getNumberDeterminants(); i++)
    {
      if (Singular(j)(i)) return true;
    }
  return false;
}
