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

#include "QMCSCFJastrow.h"

QMCSCFJastrow::QMCSCFJastrow()
{
  nalpha = 0;
  nbeta  = 0;
  wd     = 0;
  x      = 0;
}

QMCSCFJastrow::QMCSCFJastrow(QMCInput *INPUT, QMCHartreeFock *HF)
{
  initialize(INPUT,HF);
}

QMCSCFJastrow::QMCSCFJastrow(QMCInput *INPUT)
{
  QMCHartreeFock* HF = 0;
  initialize(INPUT,HF);
}

QMCSCFJastrow::QMCSCFJastrow(const QMCSCFJastrow & rhs )
{
  *this = rhs;
}

QMCSCFJastrow::~QMCSCFJastrow()
{
}

void QMCSCFJastrow::operator=(const QMCSCFJastrow & rhs )
{
  Input = rhs.Input;

  nalpha = rhs.nalpha;
  nbeta  = rhs.nbeta;

  Alpha = rhs.Alpha;
  Beta  = rhs.Beta;

  Jastrow = rhs.Jastrow;
  PE      = rhs.PE;

  Psi                    = rhs.Psi;
  Laplacian_PsiRatio     = rhs.Laplacian_PsiRatio;
  E_Local                = rhs.E_Local;
}

void QMCSCFJastrow::initialize(QMCInput *INPUT, QMCHartreeFock *HF)
{
  Input = INPUT;

  nalpha = Input->WF.getNumberElectrons(true);
  nbeta  = Input->WF.getNumberElectrons(false);

  Alpha.initialize(Input,
		   0,
		   nalpha-1,
		   true);
  
  Beta.initialize(Input,
		  nalpha,
		  Input->WF.getNumberElectrons()-1,
		  false);

  PE.initialize(Input,HF);
  Jastrow.initialize(Input);

  if(Input->flags.nuclear_derivatives != "none")
    nf.initialize(Input);

  wd     = 0;
  x      = 0;

  swTimers.allocate(1);
  swTimers(0).reset("Psi Quantities");
}

void QMCSCFJastrow::evaluate(Array1D<QMCWalkerData *> &walkerData,
			     Array1D<Array2D<double> * > &xData,
			     int num)
{
  /*
    For the GPU code, if the INT_FINISHES flags
    (set in each of the GPUQMC source files)
    are not set, then these commands will simply
    tell OpenGL what to start working
    on. That is, glFlush() as opposed to glFinish().
  */

  //When running without the GPU, gpp will recieve zero.
  int gpp = Input->flags.getNumGPUWalkers();

#ifdef QMC_GPU
  //These are pure GPU routines
  Alpha.gpuEvaluate(xData,gpp);
  Beta.gpuEvaluate(xData,gpp);
  Jastrow.setUpGPU(Alpha.gpuBF.getElectronicTexture(),
                   Beta.gpuBF.getElectronicTexture(), gpp);
#endif

  /*
    num-gpp is the number the GPU didn't calculate
    gpp is the first index the CPU needs to handle
  */
  //These are pure CPU routines
  Jastrow.evaluate(walkerData,xData,num-gpp,gpp);

  int whichE = walkerData(0)->whichE;
  if((whichE >= 0 && whichE < nalpha) || whichE == -1)
    Alpha.evaluate(xData,num-gpp,gpp,whichE);
  if(whichE >= nalpha || whichE == -1)
    Beta.evaluate(xData,num-gpp,gpp,whichE);

  if(Input->flags.nuclear_derivatives != "none")
    nf.evaluate(walkerData,xData,num);

  //These are CPU if no GPU, otherwise mixed CPU/GPU
  Alpha.update_Ds(walkerData);
  Beta.update_Ds(walkerData);

  PE.evaluate(xData,walkerData,num);
#ifdef QMC_GPU

  Jastrow.gpuEvaluate(xData,gpp);
#endif

  for(int i=0; i<num; i++)
  {
    /*
      Load the pointers; all the other functions
      will just look at these.
    */
    wd = walkerData(i);
    x  = xData(i);
    iWalker = i;

    swTimers(0).start();
    calculate_Psi_quantities();
    calculate_Modified_Grad_PsiRatio();
    calculate_E_Local(i);
    swTimers(0).lap();

    wd->localEnergy            = getLocalEnergy();
    wd->kineticEnergy          = getKineticEnergy();
    wd->potentialEnergy        = getPotentialEnergy(i);
    wd->neEnergy               = getEnergyNE(i);
    wd->eeEnergy               = getEnergyEE(i);
    wd->psi                    = getPsi();
    wd->singular              |= isSingular(i);

    wd->x2 = wd->y2 = wd->z2 = 0.0;
    for(int e=0; e<x->dim1(); e++){
      wd->x2 += x->get(e,0)*x->get(e,0);
      wd->y2 += x->get(e,1)*x->get(e,1);
      wd->z2 += x->get(e,2)*x->get(e,2);
    } 

    //checkParameterDerivatives(); 
  }

  /*
    Just to be sure nobody else tries to look at or
    modify the data behind these pointers.
  */
  wd      =  0;
  x       =  0;
  iWalker = -1;
}

void QMCSCFJastrow::evaluate(Array2D<double> &X, QMCWalkerData & data)
{
  Array1D< Array2D<double>* > temp(1);
  temp(0) = &X;

  Array1D<QMCWalkerData *> wdArray(1);
  wdArray(0) = &data;

  evaluate(wdArray,temp,1);
}

//  We end up doing a division with SCF_Psi. Therefore, it seems to be critical
//  to make sure we aren't dividing by too small a number
void QMCSCFJastrow::calculate_Psi_quantities()
{
  assert(wd);

  if(Alpha.isSingular(iWalker) || Beta.isSingular(iWalker))
    {
      /*
	Although maybe not all of the determinants are singular
	let's just assume that we don't want to actually
	salvage anything.
      */
      wd->singular = true;
      return;
    }

  term_AlphaGrad = 0;
  term_BetaGrad  = 0;
  
  alphaPsi       = Alpha.getPsi(iWalker);
  alphaGrad      = Alpha.getGradPsiRatio(iWalker);
  alphaLaplacian = Alpha.getLaplacianPsiRatio(iWalker);
  
  betaPsi        = Beta.getPsi(iWalker);
  betaGrad       = Beta.getGradPsiRatio(iWalker);
  betaLaplacian  = Beta.getLaplacianPsiRatio(iWalker);

  update_SCF();
  
  QMCDouble Jastrow_Psi = Jastrow.getJastrow(iWalker);

  Psi = wd->D * Jastrow_Psi;    

  // \frac{\nabla^2 D}{D} + \nabla^2 U
  Laplacian_PsiRatio = wd->D_xx + wd->U_xx;

  wd->gradPsiRatio           = 0.0;
  for (int i=0; i<Input->WF.getNumberElectrons(); i++)
    for (int j=0; j<3; j++)
      {
        wd->gradPsiRatio(i,j) = wd->D_x(i,j) + wd->U_x(i,j);
	
	// 2 \frac{\nabla D}{D} \cdot \nabla U
        Laplacian_PsiRatio += wd->U_x(i,j) * wd->D_x(i,j) * 2.0;
	// \nabla U \cdot \nabla U
	Laplacian_PsiRatio += wd->U_x(i,j) * wd->U_x(i,j);
      }

  if(Input->flags.calculate_Derivatives == 1)
    {
      int numAI = Input->getNumberAIParameters();
      
      if(wd->rp_a.dim1() != numAI)
	{
	  wd->rp_a.allocate(numAI);
	  wd->p3_xxa.allocate(numAI);
	}

      wd->p3_xxa = 0.0;
      
      int ai = 0;
      calculate_JastrowDerivatives(ai);
      calculate_CIDerivatives(ai);
      calculate_OrbitalDerivatives(ai);

      if(ai != numAI)
	{
	  clog << "Error: our parameter counting went bad...\n"
	       << "    ai = " << ai << endl
	       << " numAI = " << numAI << endl;
	  exit(0);
	}

      //Convert the derivative wrt to the laplacian into
      //the derivative wrt the kinetic energy
      wd->p3_xxa *= -0.5;
    }

  // Calculate the basis function density
  if (Input->flags.calculate_bf_density == 1)
  {
    wd->chiDensity = 0.0;
    if (nalpha > 0)
      wd->chiDensity = wd->chiDensity + *(Alpha.getChiDensity(iWalker));
    if (nbeta > 0)
      wd->chiDensity = wd->chiDensity + *(Beta.getChiDensity(iWalker));
  }
}

void QMCSCFJastrow::update_SCF()
{
  Array1D<QMCDouble> termPsiRatio(Input->WF.getNumberDeterminants());

  termPR.allocate(Input->WF.getNumberDeterminants());

  wd->D = 0.0;

  for (int ci=0; ci<Input->WF.getNumberDeterminants(); ci++)
    {
      termPsiRatio(ci) = QMCDouble(Input->WF.CI_coeffs(ci));
      termPsiRatio(ci).multiplyBy((*alphaPsi)(ci));
      termPsiRatio(ci).multiplyBy((*betaPsi)(ci));
      wd->D += termPsiRatio(ci);
    }
  
  if(!wd->D.isZero())
    termPsiRatio /= wd->D;
  
  for(int ci=0; ci<termPR.dim1(); ci++)
    termPR(ci) = (double)termPsiRatio(ci);
  termPsiRatio.deallocate();
  
  wd->D_xx = 0.0;  
  wd->D_x  = 0.0;
  for (int ci=0; ci<Input->WF.getNumberDeterminants(); ci++)
    {
      if(alphaGrad) term_AlphaGrad = alphaGrad->array()[ci];
      if(betaGrad)   term_BetaGrad =  betaGrad->array()[ci];
      
      if(nalpha > 0)
	wd->D_xx += termPR(ci) * (*alphaLaplacian)(ci);
      if(nbeta > 0)
	wd->D_xx += termPR(ci) * (*betaLaplacian)(ci);
      
      for (int e=0; e<nalpha; e++)
	for (int k=0; k<3; k++)
	  wd->D_x(e,k) += termPR(ci)*term_AlphaGrad[e][k];
      
      for (int e=0; e<nbeta; e++)
	for (int k=0; k<3; k++)
	  wd->D_x(e+nalpha,k) += termPR(ci)*term_BetaGrad[e][k];
    }
}

void QMCSCFJastrow::calculate_Modified_Grad_PsiRatio()
{
  // Call this after calculate_Grad_PsiRatio() is called

  double a = 0.0;
  double factor = 0.0;

  if( Input->flags.QF_modification_type == "none" )
    // do not modify the quantum force
    wd->modifiedGradPsiRatio = wd->gradPsiRatio;

  else if( Input->flags.QF_modification_type == "umrigar93_equalelectrons" )
  {
    // from Umrigar, Nightingale, and Runge JCP 99(4) 2865; 1993 eq 34.
    // The QF for each electron is changed by the same factor

    a = Input->flags.umrigar93_equalelectrons_parameter;

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
    for(int i=0; i<wd->gradPsiRatio.dim1(); i++)
      for(int j=0; j<wd->gradPsiRatio.dim2(); j++)
        magsqQF += wd->gradPsiRatio(i,j) * wd->gradPsiRatio(i,j);

    factor = ( -1.0 + sqrt( 1.0 + 2*a*magsqQF*Input->flags.dt )) /
             ( a*magsqQF*Input->flags.dt);

    wd->modifiedGradPsiRatio = wd->gradPsiRatio;
    wd->modifiedGradPsiRatio *= factor;
  }
  else if( Input->flags.QF_modification_type == "umrigar93_unequalelectrons" )
  {
    // from Umrigar, Nightingale, and Runge JCP 99(4) 2865; 1993 eq 35-36.
    // The QF for each electron is changed by a different factor

    wd->modifiedGradPsiRatio = wd->gradPsiRatio;
    Array1D<double> closest_nucleus_to_electron_norm_vector(3);
    Array1D<double> electron_velocity_norm_vector(3);
    for(int i=0; i<wd->modifiedGradPsiRatio.dim1(); i++)
    {
      // Find the closest nucleus to electron i
      int closest_nucleus_index  = Input->Molecule.findClosestNucleusIndex(wd->riA,i);
      // now we have the closest nucleus identified

      // Charge of this nucleus
      int closest_nucleus_Z =
        Input->Molecule.Z(closest_nucleus_index);

      double closest_nucleus_distance_squared = 0.0;
      double temp;

      // Unit vector from nearest nucleus to the electron
      for(int xyz=0; xyz<3; xyz++)
      {
        temp = (*x)(i,xyz) -
               Input->Molecule.Atom_Positions(closest_nucleus_index,xyz);

        closest_nucleus_to_electron_norm_vector(xyz) = temp;
        closest_nucleus_distance_squared += temp*temp;
      }
      closest_nucleus_to_electron_norm_vector *=
        1.0/sqrt(closest_nucleus_distance_squared);

      // Unit vector in the direction of electron velocity
      for(int xyz=0; xyz<3; xyz++)
        electron_velocity_norm_vector(xyz) = wd->gradPsiRatio(i,xyz);

      double electron_velocity_norm_squared =
        electron_velocity_norm_vector * electron_velocity_norm_vector;

      if (fabs(electron_velocity_norm_squared) < 1e-10)
        electron_velocity_norm_vector *= 0.0;
      else
        electron_velocity_norm_vector *= 1.0/sqrt(
                                           electron_velocity_norm_squared );

      // Now calculate the a(r) for the equation
      temp = closest_nucleus_Z * closest_nucleus_Z *
             closest_nucleus_distance_squared;

      double dotp = 0;
      for(int xyz=0; xyz<3; xyz++)
	dotp += closest_nucleus_to_electron_norm_vector(xyz)*
	  electron_velocity_norm_vector(xyz);
      a = 0.5*(1.0 + dotp) + temp/(10.0*(4.0+temp));

      //This should accomplish what the previous few lines do, but
      //valgrind thinks there's a problem. It might have to do with using
      //an array allocated with new to a fortran function...
      //a = 0.5*(1.0 + closest_nucleus_to_electron_norm_vector *
      //         electron_velocity_norm_vector) + temp/(10.0*(4.0+temp));

 
      // Now calculate the factor the QF of this electron is modified by
      if (fabs(electron_velocity_norm_squared) < 1e-10)
        factor = 1.0;
      else
        factor = (-1.0 + sqrt(1.0+2.0*a*electron_velocity_norm_squared*
                              Input->flags.dt))/
                 (a*electron_velocity_norm_squared*Input->flags.dt);

      // Now adjust the QF of electron I
      for(int xyz=0; xyz<3; xyz++)
        wd->modifiedGradPsiRatio(i,xyz) *= factor;
    }
  }
  else
  {
    cerr << "ERROR: Unknown value for QF_modification_type set!" << endl;
    exit(0);
  }

  Array2D<double> * tMGPR = & wd->modifiedGradPsiRatio;
  Array2D<double> * tGPR  = & wd->gradPsiRatio;
  
  double lengthGradTrialModified =
    sqrt((*tMGPR).dotAllElectrons(*tMGPR));
  double lengthGradTrialUnmodified =
    sqrt((*tGPR).dotAllElectrons(*tGPR));
  
  wd->modificationRatio = lengthGradTrialModified/lengthGradTrialUnmodified;
  
  if (IeeeMath::isNaN(lengthGradTrialModified) ||
      IeeeMath::isNaN(lengthGradTrialUnmodified) ||
      IeeeMath::isNaN(wd->modificationRatio) ||
      lengthGradTrialUnmodified == 0)
    {
      if(QMCPotential_Energy::printElec != -2){
	cerr << "WARNING: trial Grad Psi Ratio is NaN ";
	cerr << "   lengthGradTrialModified = " << lengthGradTrialModified << endl;
	cerr << "   lengthGradTrialUnmodified = " << lengthGradTrialUnmodified << endl;
      }
      wd->singular = true;	
    } 
}

bool QMCSCFJastrow::isSingular(int whichWalker)
{
  if(Alpha.isSingular(whichWalker) ||
     Beta.isSingular(whichWalker))
    return true;
  else
    return false;
}

void QMCSCFJastrow::calculate_JastrowDerivatives(int & ai)
{
  for(int ji=0; ji<Input->JP.getNumberJWParameters(); ji++)
    {
      wd->rp_a(ai)    = Jastrow.get_p_a_ln(iWalker,ji);
      wd->p3_xxa(ai)  = Jastrow.get_p3_xxa_ln(iWalker,ji);
      
      Array2D<double>* j_p_xa = Jastrow.get_p2_xa_ln(iWalker,ji);
      
      for (int i=0; i<Input->WF.getNumberElectrons(); i++)
	for (int j=0; j<3; j++)
	  {
	    wd->p3_xxa(ai) += 2.0 * (*j_p_xa)(i,j) * wd->D_x(i,j);
	    wd->p3_xxa(ai) += 2.0 * (*j_p_xa)(i,j) * wd->U_x(i,j);
	  }
      ai++;
    }
}

void QMCSCFJastrow::calculate_CIDerivatives(int & ai)
{
  if(Input->WF.getNumberCIParameters() <= 0)
    return;

  Array1D<double> rp_a(Input->WF.getNumberDeterminants());
  rp_a   = 0.0;
  Array1D<double> p3_xxa(Input->WF.getNumberDeterminants());
  p3_xxa = 0.0;

  Array1D<double> gradSum(Input->WF.getNumberDeterminants());
  gradSum = 0.0;

  for(int cj=0; cj<Input->WF.getNumberDeterminants(); cj++)
    {
      if(alphaGrad) term_AlphaGrad = alphaGrad->array()[cj];
      if(betaGrad)   term_BetaGrad =  betaGrad->array()[cj];
      for(int i=0; i<nalpha; i++)
	for(int j=0; j<3; j++)
	  gradSum(cj) += 2.0 * term_AlphaGrad[i][j] * wd->U_x(i,j);  
      for(int i=0; i<nbeta; i++)
	for(int j=0; j<3; j++)
	  gradSum(cj)  += 2.0 * term_BetaGrad[i][j] * wd->U_x(i+nalpha,j);

      if(nalpha > 0)
	gradSum(cj)  += (*alphaLaplacian)(cj);
      if(nbeta > 0)
	gradSum(cj)  += (*betaLaplacian)(cj);

    }

  //First, treat all the determinants as independent
  for(int ci=0; ci<Input->WF.getNumberDeterminants(); ci++)
    {
      rp_a(ci) = termPR(ci) / Input->WF.CI_coeffs(ci);

      p3_xxa(ci)  += (-termPR(ci) * termPR(ci) + termPR(ci)) * gradSum(ci); 

      for(int cj=0; cj<Input->WF.getNumberDeterminants(); cj++)
	{
	  if(ci == cj) continue;
	  p3_xxa(ci)  += -termPR(ci) * termPR(cj) * gradSum(cj);
	}
      p3_xxa(ci)  /= Input->WF.CI_coeffs(ci);
    }  
  //Some of the coefficients might not be independent
  for(int ci=0; ci<Input->WF.getNumberDeterminants(); ci++)
    {
      int c = Input->WF.CI_constraints(ci);
      if(c != -1){
	double const_coeff = Input->WF.CI_const_coeffs(ci);
	rp_a(c) += const_coeff * rp_a(ci);
	rp_a(ci) = 0.0;
	p3_xxa(c) += const_coeff * p3_xxa(ci);
	p3_xxa(ci) = 0.0;
      }
    }

  for(int ci=0; ci<Input->WF.getNumberDeterminants(); ci++)
    {
      if(fabs(rp_a(ci)) > 1e-30){
	//This is an independent parameter
	wd->rp_a(ai) = rp_a(ci);
	wd->p3_xxa(ai) = p3_xxa(ci);
	ai++;
      }
    }
}

void QMCSCFJastrow::calculate_OrbitalDerivatives(int & aic)
{
  int numDet = Input->WF.getNumberDeterminants();
  int unusedIndicator = Input->WF.getUnusedIndicator();
  Array1D<double> dtermPsi_dai(numDet);
  if(aic >= Input->getNumberAIParameters())
    return;

  Array1D<double> rp_a(Input->WF.OR_constraints.dim1());
  rp_a   = 0.0;
  Array1D<double> p3_xxa(Input->WF.OR_constraints.dim1());
  p3_xxa = 0.0;

  int ai = -1;
  for(int o=0; o<Input->WF.getNumberOrbitals(); o++)
    {
      bool used = false;
      for(int ci=0; ci<Input->WF.getNumberDeterminants(); ci++)
	{
	  if(Input->WF.AlphaOccupation(ci,o) != unusedIndicator ||
	     Input->WF.BetaOccupation(ci,o) != unusedIndicator)
	    {
	      used = true;
	      break;
	    }
	}
      
      //The orbital might not be used at all
      if(!used) continue;
      
      //If it is used, then we can add the BF to the loop
      for(int bf=0; bf<Input->WF.getNumberBasisFunctions(); bf++)
	{
	  ai++;
	  if(globalInput.WF.OR_constraints(ai) == -2)
	    continue;

	  dtermPsi_dai = 0.0;
	  double dSCF_Psi_dai = 0.0;
	  
	  for(int ci=0; ci<numDet; ci++)
	    {
	      int oa = Input->WF.AlphaOccupation(ci,o);
	      int ob = Input->WF.BetaOccupation(ci,o);
	      
	      double dtermPsi_a = 0.0;
	      if(oa != unusedIndicator)
		dtermPsi_a += Alpha.get_p_a(iWalker,ci)->get(oa,bf) * (*betaPsi)(ci);
	      if(ob != unusedIndicator)
		dtermPsi_a += Beta.get_p_a(iWalker,ci)->get(ob,bf) * (*alphaPsi)(ci);
	      
	      //This determinant might not depend upon this orbital
	      if(dtermPsi_a == 0) continue;
	      
	      dtermPsi_a *= Input->WF.CI_coeffs(ci);
	      dtermPsi_dai(ci) = dtermPsi_a;
	      dSCF_Psi_dai += dtermPsi_a;
	    }
	  
	  rp_a(ai) = dSCF_Psi_dai / wd->D;
	  
	  for(int ci=0; ci<numDet; ci++)
	    {
	      int oa = Input->WF.AlphaOccupation(ci,o);
	      int ob = Input->WF.BetaOccupation(ci,o);
	      
	      // d (N/D) = dN / D - N dD / D^2 ; N = termPsiRatio(ci), D = SCF_Psi
	      double dtermPsiRatio_dai;
	      dtermPsiRatio_dai = (dtermPsi_dai(ci) - termPR(ci) * dSCF_Psi_dai) / wd->D;
	      
	      // d ( N * (A + B) ) = dN * (A + B) + N * (dA + dB)
	      if(nalpha > 0)
		p3_xxa(ai) += dtermPsiRatio_dai * (*alphaLaplacian)(ci);
	      if(nbeta > 0)
		p3_xxa(ai) += dtermPsiRatio_dai * (*betaLaplacian)(ci);
	      
	      if(oa != unusedIndicator)
		p3_xxa(ai) += termPR(ci) * Alpha.get_p3_xxa(iWalker,ci)->get(oa,bf);
	      if(ob != unusedIndicator)
		p3_xxa(ai) += termPR(ci) * Beta.get_p3_xxa(iWalker,ci)->get(ob,bf);
	      
	      //And now add in the terms from the gradient
	      double Term_dSlater_a = 0.0;
	      if(oa != unusedIndicator)
		for(int i=0; i<nalpha; i++)
		  for(int j=0; j<3; j++)
		    Term_dSlater_a += Alpha.get_p2_xa(iWalker,ci,i,j)->get(oa,bf) * wd->U_x(i,j);
	      if(ob != unusedIndicator)
		for(int i=0; i<nbeta; i++)
		  for(int j=0; j<3; j++)
		    Term_dSlater_a += Beta.get_p2_xa(iWalker,ci,i,j)->get(ob,bf) * wd->U_x(i+nalpha,j);
	      Term_dSlater_a *= 2.0 * termPR(ci);
	      
	      if(alphaGrad) term_AlphaGrad = alphaGrad->array()[ci];
	      if(betaGrad)   term_BetaGrad =  betaGrad->array()[ci];
	      
	      double Slater_dTerm_a = 0.0;
	      for(int i=0; i<nalpha; i++)
		for(int j=0; j<3; j++)
		  Slater_dTerm_a += term_AlphaGrad[i][j] * wd->U_x(i,j);
	      for(int i=0; i<nbeta; i++)
		for(int j=0; j<3; j++)
		  Slater_dTerm_a += term_BetaGrad[i][j]  * wd->U_x(i+nalpha,j);
	      Slater_dTerm_a *= 2.0 * dtermPsiRatio_dai;
	      
	      p3_xxa(ai) += Slater_dTerm_a + Term_dSlater_a;
	    }
	}
    }

  /*
    Above, we stored the derivatives in lists which ignored the constraints.
  */
  for(int ori=0; ori<Input->WF.OR_constraints.dim1(); ori++)
    {
      int c = Input->WF.OR_constraints(ori);
      if(c > 0){
	rp_a(c) += rp_a(ori);
	rp_a(ori) = 0.0;
	p3_xxa(c) += p3_xxa(ori);
	p3_xxa(ori) = 0.0;
      }
    }

  for(int ori=0; ori<Input->WF.OR_constraints.dim1(); ori++)
    {
      if(fabs(rp_a(ori)) > 1e-30){
	//This is an independent parameter
	wd->rp_a(aic) = rp_a(ori);
	wd->p3_xxa(aic) = p3_xxa(ori);
	aic++;
      }
    }
}

void QMCSCFJastrow::checkParameterDerivatives()
{
  if(Input->flags.iseed == 0)
    {
      cout << "Error: to check derivatives, iseed must not be 0\n";
      exit(0);
    }
  if(Input->flags.calculate_Derivatives != 1)
    {
      cout << "Error: to check derivatives, you must calculate derivatives!\n";
      exit(0);
    }
  globalInput.printAISummary();
  
  int aStart = 0;
  int aStop  = Input->getNumberAIParameters();
  Array1D<double> params = Input->getAIParameters();
  globalInput.printAIParameters(cout,"Parameters:",25,params,false);
  printf("function %23c p   %25.15e p2_xx  %25.15e\n",
	 ' ',
	 (double)wd->psi,
	 wd->localEnergy);
  for(int ai=aStart; ai<aStop; ai++)
    printf("ai %3i %25.15e p_a %25.15e p3_xxa %25.15e\n",
	   ai,
	   params(ai),
	   (double)(Psi) * wd->rp_a(ai),
	   wd->p3_xxa(ai));
  exit(0);
}

void QMCSCFJastrow::calculate_CorrelatedSampling(Array1D<QMCWalkerData *> &walkerData,
						 Array1D<Array2D<double> * > &xData,
						 int num)
{
  for(int cs = globalInput.cs_Parameters.dim1()-1; cs >= 0; cs--)
    {
      globalInput.setAIParameters( globalInput.cs_Parameters(cs) );
      
      //int gpp = Input->flags.getNumGPUWalkers();
      int gpp = 0;

      if(globalInput.flags.optimize_Orbitals == 1)
	{
	  //We don't need to revaluate the basis functions, so there is room
	  //for efficiency improvement here
	  Alpha.evaluate(xData,num-gpp,gpp,-1);
	  Beta.evaluate(xData,num-gpp,gpp,-1);

	  Alpha.update_Ds(walkerData);
	  Beta.update_Ds(walkerData);
	}
      
      /*
	jastrow whichE needs to be -1 since i'm not going to program
	ways to save intermediate data
      */
      if(globalInput.JP.getNumberJWParameters() > 0)
	Jastrow.evaluate(walkerData,xData,num-gpp,gpp);
      
      for(int w=0; w<num; w++)
	{
	  iWalker = w;
	  wd      = walkerData(iWalker);
	  x       = xData(iWalker);

	  if(globalInput.flags.optimize_CI == 1 ||
	     globalInput.flags.optimize_Orbitals == 1)
	    {
	      alphaPsi       = Alpha.getPsi(iWalker);
	      alphaGrad      = Alpha.getGradPsiRatio(iWalker);
	      alphaLaplacian = Alpha.getLaplacianPsiRatio(iWalker);
	      
	      betaPsi        = Beta.getPsi(iWalker);
	      betaGrad       = Beta.getGradPsiRatio(iWalker);
	      betaLaplacian  = Beta.getLaplacianPsiRatio(iWalker);
	      
	      update_SCF();
	    }
	  
	  QMCDouble Jastrow_Psi = Jastrow.getJastrow(w);
	  QMCDouble cs_Psi = wd->D * Jastrow_Psi;

	  // \frac{\nabla^2 D}{D} + \nabla^2 U
	  double cs_LPR = wd->D_xx + wd->U_xx;
	  
	  for (int i=0; i<Input->WF.getNumberElectrons(); i++)
	    for (int j=0; j<3; j++)
	      {
		// 2 \frac{\nabla D}{D} \cdot \nabla U
		cs_LPR += wd->U_x(i,j) * wd->D_x(i,j) * 2.0;
		// \nabla U \cdot \nabla U
		cs_LPR += wd->U_x(i,j) * wd->U_x(i,j);
	      }

	  QMCDouble weight = cs_Psi / wd->psi;
	  weight *= weight;

	  if(wd->cs_LocalEnergy.dim1() != globalInput.cs_Parameters.dim1())
	    {
	      wd->cs_LocalEnergy.allocate(globalInput.cs_Parameters.dim1());
	      wd->cs_Weights.allocate(globalInput.cs_Parameters.dim1());
	    }

	  wd->cs_LocalEnergy(cs) = -0.5 * cs_LPR + wd->potentialEnergy;
	  wd->cs_Weights(cs)     = (double)(weight);
	}
    }

  /*
    Just to be sure nobody else tries to look at or
    modify the data behind these pointers.
  */
  wd      =  0;
  x       =  0;
  iWalker = -1;
}

int QMCSCFJastrow::getNumTimers()
{
  int numTimers = 1 +
    Alpha.getNumTimers() +
    Beta.getNumTimers() +
    PE.getNumTimers() +
    Jastrow.getNumTimers();

  return numTimers;
}

void QMCSCFJastrow::aggregateTimers(Array1D<Stopwatch> & timers,
				    int & idx)
{
  //int idxStart = idx;
  Alpha.aggregateTimers(timers,idx);
  //idx=idxStart;
  Beta.aggregateTimers(timers,idx);
  PE.aggregateTimers(timers,idx);
  Jastrow.aggregateTimers(timers,idx);

  for(int i=0; i<swTimers.dim1(); i++)
    timers(idx++).aggregateTimer(swTimers(i));
}
