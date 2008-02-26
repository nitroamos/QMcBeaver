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
    {
      Array2D<qmcfloat> * coeffs = WF->getCoeff(isAlpha);
      ElectronNucleusCusp.initialize(Input, * coeffs);
      ElectronNucleusCusp.fitReplacementOrbitals();
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
      if(WF->getNumberDeterminants() > 1)
	ciDet = new Array2D<qmcfloat>();
      D.allocate(wpp-gpp);
      D_inv.allocate(wpp-gpp,ndet);
      Laplacian_D.allocate(wpp-gpp);
      Grad_D.allocate(wpp-gpp,3);
    }

#ifdef QMC_GPU
  gpu_D.allocate(gpp);
  gpu_D_inv.allocate(gpp);
  gpu_Laplacian_D.allocate(gpp);
  gpu_Grad_D.allocate(gpp,3);
  
  //this code is now out of sync when we move
  //one e at a time.
  assert(0);
  for(int j=0; j<gpu_D.dim1(); j++)
    {
      for (int i=0; i<gpu_D.dim2(); i++)
	gpu_D_inv(j,i).allocate(N,N);
      
      gpu_D(j).allocate(N,N);
      
      gpu_Laplacian_D(j).allocate(N,N);
      for(int k=0; k<3; k++)
	gpu_Grad_D(j,k).allocate(N,N);
    }
#endif

  Singular.allocate(wpp);
  Psi.allocate(wpp);
  Laplacian_PsiRatio.allocate(wpp);
  Grad_PsiRatio.allocate(wpp);

  pointer_det.allocate(wpp);
  pointer_lapPR.allocate(wpp);
  pointer_gradPR.allocate(wpp);

  if (Input->flags.calculate_bf_density == 1)
    Chi_Density.allocate(wpp);

  for(int w=0; w<wpp; w++)
    {
      Singular(w).allocate(ndet);
      Psi(w).allocate(ndet);
      Laplacian_PsiRatio(w).allocate(ndet);
      Grad_PsiRatio(w).allocate(ndet,N,3);
      if (Input->flags.calculate_bf_density == 1)
        Chi_Density(w).allocate(nbasisfunc);
    }

  if (Input->flags.optimize_Orbitals == 1)
    {
      p_a.allocate(wpp);
      p2_xa.allocate(wpp);
      p3_xxa.allocate(wpp);

      for(int w=0; w<wpp; w++)
        {
          p_a(w).allocate(ndet);
          p2_xa(w).allocate(ndet,N,3);
          p3_xxa(w).allocate(ndet);

          for(int pi=0; pi<ndet; pi++)
            {
              (p_a(w))(pi).allocate(N,nbasisfunc);
              (p3_xxa(w))(pi).allocate(N,nbasisfunc);

              for(int pj=0; pj<N; pj++)
                for(int pk=0; pk<3; pk++)
                  (p2_xa(w))(pi,pj,pk).allocate(N,nbasisfunc);
            }
        }
    }

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
}

void QMCSlater::allocateIteration(int whichE, int & start, int & stop)
{
  int nOR = WF->getNumberOrbitals();
  int nEL = getNumberElectrons();
  int nBF = WF->getNumberBasisFunctions();
  int wpp = Input->flags.walkers_per_pass;

  //nUp is the number of electrons updated per iteration.
  int nUP = whichE == -1 ? nEL : 1;

  if(whichE == -1)
    {
      start = Start;
      stop  = Stop;
    }
  else
    {
      start = whichE;
      stop  = whichE;
    }

  bool needToKeepChi = false;
  if(globalInput.flags.optimize_Orbitals)
    needToKeepChi = true;

  if(needToKeepChi)
    {
      /*
	If we're optimizing the orbitals, then we will
	need to keep Chi data a little bit longer than
	otherwise.
      */
      Chi.allocate(wpp);
      Chi_laplacian.allocate(wpp);
      Chi_gradient.allocate(wpp);
    }
  else
    {
      Chi.allocate(1);
      Chi_laplacian.allocate(1);
      Chi_gradient.allocate(1);
    }

  for(int w=0; w<Chi.dim1(); w++)
    {
      Chi(w).allocate(nUP,nBF);
      Chi_laplacian(w).allocate(nUP,nBF);
      Chi_gradient(w).allocate(3);
      for(int xyz=0; xyz<3; xyz++)
	(Chi_gradient(w))(xyz).allocate(nUP,nBF);
    }

  for(int w=0; w<D.dim1(); w++)
    {
      for (int ci=0; ci<WF->getNumberDeterminants(); ci++)
	D_inv(w,ci).allocate(nEL,nEL);
      
      D(w).allocate(nUP,nOR);
      
      Laplacian_D(w).allocate(nUP,nOR);
      for(int xyz=0; xyz<3; xyz++)
	Grad_D(w,xyz).allocate(nUP,nOR);
    }
}

QMCSlater::~QMCSlater()
{
  if (Start != Stop)
    {
      for(int j=0; j<D.dim1(); j++)
        {
          for (int ci=0; ci<WF->getNumberDeterminants(); ci++)
	    D_inv(j,ci).deallocate();

	  D(j).deallocate();
	  Laplacian_D(j).deallocate();
	  for(int k=0; k<3; k++)
	    Grad_D(j,k).deallocate();
        }
      if(ciDet != 0)
	delete ciDet;
      D.deallocate();
      D_inv.deallocate();
      Laplacian_D.deallocate();
      Grad_D.deallocate();

#ifdef QMC_GPU
      for(int j=0; j<gpu_D.dim1(); j++)
        {
          for (int i=0; i<WF->getNumberDeterminants(); i++)
	    gpu_D_inv(j,i).deallocate();
	  
          
	  gpu_D(j).deallocate();
	  
	  gpu_Laplacian_D(j).deallocate();
	  for(int k=0; k<3; k++)
	    gpu_Grad_D(j,k).deallocate();
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

void QMCSlater::update_Ds(Array1D< QMCWalkerData * > & walkerData)
{
  int whichE = walkerData(0)->whichE;
  int gpp = Input->flags.getNumGPUWalkers();
  int wpp = walkerData.dim1();

  if(getNumberElectrons() == 0)
    return;

  for(int w=0; w<wpp; w++)
    if(globalInput.flags.one_e_per_iter)
      {
	/*
	  Even if we didn't move any of our electrons, we will
	  still be asked for our data, so we need to assign the
	  pointers.
	*/
	if(isAlpha)
	  {
	    pointer_det(w) =    & walkerData(w+gpp)->PsiA;
	    pointer_gradPR(w) = & walkerData(w+gpp)->Grad_PsiRatioA;
	    pointer_lapPR(w) =  & walkerData(w+gpp)->Laplacian_PsiRatioA;
	  } else {
	    pointer_det(w) =    & walkerData(w+gpp)->PsiB;
	    pointer_gradPR(w) = & walkerData(w+gpp)->Grad_PsiRatioB;
	    pointer_lapPR(w) =  & walkerData(w+gpp)->Laplacian_PsiRatioB;
	  }
      } else {
	pointer_det(w) =    & Psi(w);
	pointer_gradPR(w) = & Grad_PsiRatio(w);
	pointer_lapPR(w) =  & Laplacian_PsiRatio(w);
      }
      
  if(whichE != -1)
    if(whichE < Start || whichE > Stop)
      return;

  for(int w=0; w<wpp; w++)
    {
      //walker offset index
      int off;	  
      Array2D<double> * inv;
      Array2D<qmcfloat> * gradx, * grady, * gradz, * lap;

#ifdef QMC_GPU
      //i changed a lot in qmcslater to move one e per iter
      //this bit of code will need to be updated when
      //i have the GPU stuff in front of me again.
      assert(0);
      
      /* This will fill the resultsCollector with
	 the results from the GPU. However, all of 
	 resultsCollector's elements are really
	 just pointers to the D object.*/
      off = w;
      gpuMatMult.getResults(resultsCollector);
      calculate_DerivativeRatios(0,ci,gpu_D,gpu_D_inv,
				 gpu_Laplacian_D,gpu_Grad_D);
      
#endif
      off = w + gpp;
      
      if(globalInput.flags.one_e_per_iter)
	{
	  //it might have been an all electron update
	  int numE = D(w).dim1();
	  //index of first updated electron
	  int iEl = whichE == -1 ? 0 : whichE - Start;
	  
	  //reload the data from the other electrons
	  if(isAlpha)
	    {
	      lap =   &  walkerData(off)->Laplacian_DA;
	      gradx = & (walkerData(off)->Grad_DA)(0);
	      grady = & (walkerData(off)->Grad_DA)(1);
	      gradz = & (walkerData(off)->Grad_DA)(2);
	    } else {
	      lap =   &  walkerData(off)->Laplacian_DB;
	      gradx = & (walkerData(off)->Grad_DB)(0);
	      grady = & (walkerData(off)->Grad_DB)(1);
	      gradz = & (walkerData(off)->Grad_DB)(2);
	    }

	  //and update the data for our electron
	  lap->setRows(iEl,0,numE,Laplacian_D(w));
	  gradx->setRows(iEl,0,numE,Grad_D(w,0));
	  grady->setRows(iEl,0,numE,Grad_D(w,1));
	  gradz->setRows(iEl,0,numE,Grad_D(w,2));
	  
	} else {
	  lap =   & Laplacian_D(w);
	  gradx = & Grad_D(w,0);
	  grady = & Grad_D(w,1);
	  gradz = & Grad_D(w,2);
	}
      
      for(int ci=0; ci<WF->getNumberDeterminants(); ci++)
	{
	  if(globalInput.flags.one_e_per_iter)
	    {
	      if(isAlpha) inv =   &  walkerData(off)->D_invA(ci);
	      else        inv =   &  walkerData(off)->D_invB(ci);

	      if(ci > 0 && whichE < 0)
		{
		  //We need to copy the previous inverse during all electron
		  //updates
		  if(isAlpha) (*inv) = walkerData(off)->D_invA(ci-1);
		  else        (*inv) = walkerData(off)->D_invB(ci-1);
		}
	    } else {
	      //this might only be used by orbital derivatives...
	      inv =   & D_inv(w,ci);
	      
	      //this part is used when updating the inverse,
	      //as opposed to calculating all of it
	      if(ci > 0)
		(*inv) = D_inv(w,ci-1);
	    }

	  (Singular(off))(ci) = calculate_DerivativeRatios(ci,whichE - Start,
							   D(w),*inv,*lap,
							   *gradx, *grady, *gradz,
							   *pointer_det(w),
							   *pointer_gradPR(w),
							   *pointer_lapPR(w));
	}
    }

  if(Input->flags.optimize_Orbitals)
    {
#if defined SINGLEPRECISION || defined QMC_GPU
      /*
	It's tedious to get qmcfloat variables to cooperate here,
	so i'm not going to program it. it's non essential.
	Just don't compile in single precision.
      */
      cout << "Error: Don't use single precision to optimize the orbitals.\n";
      //It's sure to crash eventually...
      exit(0);
#else
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
      Array2D<double> ciData(getNumberElectrons(),getNumberElectrons());

      for(int w=0; w<D.dim1(); w++)
        {
          for(int ci=0; ci<WF->getNumberDeterminants(); ci++)
            {
              Array2D<double> * pa = & (p_a(w+gpp))(ci);
              Array2D<double> * p3xxa = & (p3_xxa(w+gpp))(ci);
              * pa    = 0.0;
              * p3xxa = 0.0;

              //calculate  d |D|
              D_inv(w,ci).gemm(Chi(w),*pa,false);
              (*pa) *= (Psi(w+gpp))(ci);

              for(int xyz=0; xyz<3; xyz++)
                {
		  WF->getDataForCI(ci,isAlpha,Grad_D(w,xyz),ciData);
                  ciData.gemm(D_inv(w,ci),temp,false);
                  temp.gemm(Chi(w),nabDinvDdD(xyz),false);
                  //temp = 0.0;
                }
	      WF->getDataForCI(ci,isAlpha,Laplacian_D(w),ciData);
              ciData.gemm(D_inv(w,ci),temp,false);
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
#endif
    }
}

template <class T>
bool QMCSlater::calculate_DerivativeRatios(int ci, int row,
					   Array2D<T> & psi,
					   Array2D<double> & inv,
					   Array2D<T> & lap,
					   Array2D<T> & gradx,
					   Array2D<T> & grady,
					   Array2D<T> & gradz,
					   Array1D<double> & det,
					   Array3D<double> & gradPR,
					   Array1D<double> & lapPR)
{
  Array2D<int> * occ = WF->getOccupations(isAlpha);
  Array2D<int> * swap = WF->getDeterminantSwaps(isAlpha);

  bool calcOK = true;
  double ratio = 0;

  if(WF->getNumberDeterminants() == 1)
    {
      //ciDet is either used as a pointer to the ci matrix
      ciDet = & psi;
    } else {
      //or as storage for each ci matrix
      ciDet->allocate(psi.dim1(),inv.dim1());
      WF->getDataForCI(ci,isAlpha,psi,*ciDet);
    }

  if(ciDet->dim1() > 1 || row < 0)
    {
      if(ci == 0)
	{
	  /* The LU Decomposition inverse used here is O(1*N^3)
	     Notice: after this, psi is lost!
	  */
	  ciDet->determinant_and_inverse(inv,det(ci),&calcOK);
	} else {
	  /* Each determinant is different from the previous determinant
	     by probably only a couple of orbitals. So instead of performing
	     the entire inverse, we use the Sherman-Morrison formula to update
	     the previous inverse.
	     Note: In order for this to work, inv needs to already be set to the previous inv
	  */
	  det(ci) = det(ci-1);
	  for(int o=0; o<swap->dim2(); o++)	  
	    {
	      int so = swap->get(ci,o);
	      if(so != -1)
		{				
		  ratio = inv.inverseUpdateOneColumn(so,*ciDet,so);
		  
		  det[ci] *= ratio;
		  if(ratio == 0) break;
		}
	    }
	  
	  //The update might have failed, so we could try inverting from scratch...
	  if(det(ci) == 0)
	    {
	      calcOK = false;
	      ciDet->determinant_and_inverse(inv,det(ci),&calcOK);
	    }
	}
    }
  else
    {
      // Updating one electron at a time is O(2*N^3-N^2)
      ratio = inv.inverseUpdateOneRow(row,*ciDet,0);
      /* If ratio is 0, then we don't want to update the determinant, since
	 we'll just try another move
      */
      if(ratio != 0.0)
	det[ci] *= ratio;
      else
	calcOK = false;
    }

  if(WF->getNumberDeterminants() == 1)    
    ciDet = 0;

  if(!calcOK) return true;

  int numE = getNumberElectrons();
  double &   lap_PR = lapPR[ci];
  double ** grad_PR = gradPR.array()[ci];
  lap_PR = 0.0;

  for(int e=0; e<numE; e++)
    {
      grad_PR[e][0] = 0.0;
      grad_PR[e][1] = 0.0;
      grad_PR[e][2] = 0.0;
    }

  for(int o=0; o<numE; o++)
    {
      int orb = occ->get(ci,o);

      for(int e=0; e<numE; e++)
	{
	  lap_PR        +=   lap[e][orb] * inv[o][e];
	  grad_PR[e][0] += gradx[e][orb] * inv[o][e];
	  grad_PR[e][1] += grady[e][orb] * inv[o][e];
	  grad_PR[e][2] += gradz[e][orb] * inv[o][e];
	}
    }
  return false;
}

#if defined SINGLEPRECISION || defined QMC_GPU
template bool QMCSlater::calculate_DerivativeRatios<float>
(int ci, int row,
 Array2D<float> & psi,
 Array2D<double> & inv,
 Array2D<float> & lap,
 Array2D<float> & gradx,
 Array2D<float> & grady,
 Array2D<float> & gradz,
 Array1D<double> & det,
 Array3D<double> & gradPR,
 Array1D<double> & lapPR);
#endif

template bool QMCSlater::calculate_DerivativeRatios<double>
(int ci, int row,
 Array2D<double> & psi,
 Array2D<double> & inv,
 Array2D<double> & lap,
 Array2D<double> & gradx,
 Array2D<double> & grady,
 Array2D<double> & gradz,
 Array1D<double> & det,
 Array3D<double> & gradPR,
 Array1D<double> & lapPR);
  
/**
   This function contains the meat of the QMC calculation.
   1) Basis functions are calculated for each electron position
   for each determinant:
   2) The Slater matrix is calculated along with the associated gradient and 
   laplacian matricies:
   
   D(numElec x numOrb) = Chi(numElec x numBasisFunction) * 
                         WF_coeffs(numBasisFunction x numOrb)
 
  Note: the coefficient matricies are transposed.  Hand-coded matrix 
  multiplication is faster with a transposed matrix.
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

void QMCSlater::evaluate(Array1D<Array2D<double>*> &X,
                         int num, int start, int whichE)
{
  if(getNumberElectrons() == 0)
    return;

  //if this iteration doesn't involve updating
  //one of our electrons, we can return
  if(whichE != -1)
    if(whichE < Start || whichE > Stop)
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

      int Chi_i = Chi.dim1() == 1 ? 0 : walker;
      int startE, stopE;

      allocateIteration(whichE,startE,stopE);
      
      BF->evaluateBasisFunctions(*X(walker+start),
                                 startE,stopE,
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

      Array2D<qmcfloat> * coeffs = WF->getCoeff(isAlpha);

      Chi(Chi_i).gemm( *coeffs, D(walker),true);
      Chi_laplacian(Chi_i).gemm( *coeffs, Laplacian_D(walker),true);
      (Chi_gradient(Chi_i))(0).gemm( *coeffs, Grad_D(walker,0),true);
      (Chi_gradient(Chi_i))(1).gemm( *coeffs, Grad_D(walker,1),true);
      (Chi_gradient(Chi_i))(2).gemm( *coeffs, Grad_D(walker,2),true);

      if (Input->flags.replace_electron_nucleus_cusps == 1)
	ElectronNucleusCusp.replaceCusps(*X(walker+start),
					 startE,stopE,
					 D(walker),
					 Grad_D(walker,0),
					 Grad_D(walker,1),
					 Grad_D(walker,2),
					 Laplacian_D(walker));

      if(showTimings)
        {
          sw.stop();
          timeM += sw.timeUS();
          if(numT >= 0) averageM += sw.timeUS();
        }
    }

  for(int walker = 0; walker < num; walker++)
    {
      int Chi_i = Chi.dim1() == 1 ? 0 : walker;

      if (Input->flags.calculate_bf_density == 1)
        {
          Chi_Density(walker+start) = 0.0;
          for (int i=0; i<WF->getNumberBasisFunctions(); i++)
            for (int j=0; j<getNumberElectrons(); j++)
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

Array1D<double>* QMCSlater::getPsi(int i)
{
  if(getNumberElectrons() != 0)
    return pointer_det(i);

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
    return pointer_lapPR(i);
  return 0;
}

Array3D<double>* QMCSlater::getGradPsiRatio(int i)
{
  if(getNumberElectrons() != 0)
    return pointer_gradPR(i);

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
