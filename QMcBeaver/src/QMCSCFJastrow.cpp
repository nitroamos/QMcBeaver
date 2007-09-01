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

  nalpha = Input->WF.getNumberAlphaElectrons();
  nbeta  = Input->WF.getNumberBetaElectrons();

  Alpha.initialize(Input,
		   0,
		   Input->WF.getNumberAlphaElectrons()-1,
		   true);
  
  Beta.initialize(Input,
		  Input->WF.getNumberAlphaElectrons(),
		  Input->WF.getNumberElectrons()-1,
		  false);

  PE.initialize(Input,HF);
  Jastrow.initialize(Input);

  if(Input->flags.nuclear_derivatives != "none")
    nf.initialize(Input);

  wd     = 0;
  x      = 0;
}

void QMCSCFJastrow::evaluate(Array1D<QMCWalkerData *> &walkerData,
                            Array1D<Array2D<double> * > &xData,
                            int num)
{
  //These are some useful benchmarking variables
  //static const bool printKE = false;
  static const bool showTimings = !true;
  static double averageT = 0, timeT = 0;
  static double numT = -5;
  Stopwatch sw = Stopwatch();
  if(showTimings)
  {
    sw.reset();
    sw.start();
  }

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
  
  PE.evaluate(xData,walkerData,num);

  if(Input->flags.nuclear_derivatives != "none")
    nf.evaluate(walkerData,xData,num);

#ifdef QMC_GPU

  if(true)
  {
    Stopwatch finisher = Stopwatch();
    finisher.reset();
    finisher.start();
    glFinish();
    finisher.stop();
    printf("finish time: %15d\n", finisher.timeUS() );
  }
#endif

  //These are CPU if no GPU, otherwise mixed CPU/GPU
  Alpha.update_Ds(walkerData);
  Beta.update_Ds(walkerData);
#ifdef QMC_GPU

  Jastrow.gpuEvaluate(xData,gpp);
#endif

  if(showTimings)
  {
    sw.stop();
    timeT = sw.timeUS();
    if(numT >= 0)
      averageT += sw.timeUS();
    if( num > 1)
      numT++;
    cout << "total time: " << (int)(timeT+0.5) << " ( " << (int)(averageT/(numT)+0.5) << ")\n";
  }

  for(int i=0; i<num; i++)
  {
    /*
      Load the pointers; all the other functions
      will just look at these.
    */
    wd = walkerData(i);
    x  = xData(i);

    calculate_Psi_quantities(i);
    calculate_Modified_Grad_PsiRatio();
    calculate_E_Local(i);

    wd->localEnergy            = getLocalEnergy();
    wd->kineticEnergy          = getKineticEnergy();
    wd->potentialEnergy        = getPotentialEnergy(i);
    wd->neEnergy               = getEnergyNE(i);
    wd->eeEnergy               = getEnergyEE(i);
    wd->psi                    = getPsi();
    wd->lnJ                    = Jastrow.getLnJastrow(i);
    wd->isSingular            |= isSingular(i);
  }

  /*
    Just to be sure nobody else tries to look at or
    modify the data behind these pointers.
  */
  wd = 0;
  x  = 0;
}

void QMCSCFJastrow::evaluate(Array2D<double> &X, QMCWalkerData & data)
{
  Array1D< Array2D<double>* > temp(1);
  temp(0) = &X;

  Array1D<QMCWalkerData *> wdArray(1);
  wdArray(0) = &data;
  int whichE = data.whichE;

#ifdef QMC_GPU
  //These are pure GPU
  Alpha.gpuEvaluate(temp,1);
  Beta.gpuEvaluate(temp,1);
#else
  //These are pure CPU
  if((whichE >= 0 && whichE < nalpha) || whichE == -1)
    Alpha.evaluate(temp,1,0,whichE);
  if(whichE >= nalpha || whichE == -1)
    Beta.evaluate(temp,1,0,whichE);

#endif

  Alpha.update_Ds(wdArray);
  Beta.update_Ds(wdArray);
  //We just let the CPU initialize the Jastrow
  Jastrow.evaluate(wdArray,temp,1,0);
  PE.evaluate(temp,wdArray,1);

  wd  = & data;
  x   = & X;
  calculate_Psi_quantities(0);
  calculate_Modified_Grad_PsiRatio();
  calculate_E_Local(0);

  if(Input->flags.nuclear_derivatives != "none")
    nf.evaluate(data,X);

  data.isSingular = isSingular(0);
  data.localEnergy = getLocalEnergy();
  data.kineticEnergy = getKineticEnergy();
  data.potentialEnergy = getPotentialEnergy(0);
  data.eeEnergy = getEnergyEE(0);
  data.neEnergy = getEnergyNE(0);
  data.psi = getPsi();
}

//  We end up doing a division with SCF_Psi. Therefore, it seems to be critical
//  to make sure we aren't dividing by too small a number
void QMCSCFJastrow::calculate_Psi_quantities(int walker)
{
  assert(wd);

  if(Alpha.isSingular(walker) || Beta.isSingular(walker))
    {
      /*
	Although maybe not all of the determinants are singular
	let's just assume that we don't want to actually
	salvage anything.
      */
      wd->isSingular = true;
      return;
    }
  Array1D<QMCGreensRatioComponent> termPsiRatio(Input->WF.getNumberDeterminants());
  Array1D<double> termPR(Input->WF.getNumberDeterminants());

  double** term_AlphaGrad = 0;
  double** term_BetaGrad  = 0;
  
  Array1D<double>* alphaPsi       = Alpha.getPsi(walker);
  Array3D<double>* alphaGrad      = Alpha.getGradPsiRatio(walker);
  Array1D<double>* alphaLaplacian = Alpha.getLaplacianPsiRatio(walker);
  
  Array1D<double>* betaPsi        = Beta.getPsi(walker);
  Array3D<double>* betaGrad       = Beta.getGradPsiRatio(walker);
  Array1D<double>* betaLaplacian  = Beta.getLaplacianPsiRatio(walker);

  Array2D<double>* JastrowGrad    = Jastrow.getGradientLnJastrow(walker);
  
  QMCGreensRatioComponent Jastrow_Psi =
    QMCGreensRatioComponent(1.0,1.0,0.0,Jastrow.getLnJastrow(walker));

  QMCGreensRatioComponent SCF_Psi = 0.0;

  for (int ci=0; ci<Input->WF.getNumberDeterminants(); ci++)
    {
      termPsiRatio(ci) = QMCGreensRatioComponent(Input->WF.CI_coeffs(ci));
      termPsiRatio(ci).multiplyBy((*alphaPsi)(ci));
      termPsiRatio(ci).multiplyBy((*betaPsi)(ci));
      SCF_Psi += termPsiRatio(ci);
    }

  termPsiRatio /= SCF_Psi;
  
  for(int ci=0; ci<termPR.dim1(); ci++)
    termPR(ci) = (double)termPsiRatio(ci);
  termPsiRatio.deallocate();

  double ae = 0;
  double be = 0;
  wd->SCF_Laplacian_PsiRatio = 0.0;  
  wd->SCF_Grad_PsiRatio      = 0.0;
  for (int ci=0; ci<Input->WF.getNumberDeterminants(); ci++)
    {
      if(alphaGrad) term_AlphaGrad = alphaGrad->array()[ci];
      if(betaGrad)   term_BetaGrad =  betaGrad->array()[ci];
      
      if(nalpha > 0)
	wd->SCF_Laplacian_PsiRatio += termPR(ci) * (*alphaLaplacian)(ci);
      if(nbeta > 0)
	wd->SCF_Laplacian_PsiRatio += termPR(ci) * (*betaLaplacian)(ci);

      for (int e=0; e<nalpha; e++)
	for (int k=0; k<3; k++)
	  wd->SCF_Grad_PsiRatio(e,k) += termPR(ci)*term_AlphaGrad[e][k];

      for (int e=0; e<nbeta; e++)
	for (int k=0; k<3; k++)
	  wd->SCF_Grad_PsiRatio(e+nalpha,k) += termPR(ci)*term_BetaGrad[e][k];
    }
  
  Psi = SCF_Psi * Jastrow_Psi;    

  Laplacian_PsiRatio = wd->SCF_Laplacian_PsiRatio +
    Jastrow.getLaplacianLnJastrow(walker);

  //double PsuedoForce2_PsiRatio = 0;
  wd->gradPsiRatio           = 0.0;
  for (int i=0; i<Input->WF.getNumberElectrons(); i++)
    for (int j=0; j<3; j++)
      {
        wd->gradPsiRatio(i,j) = wd->SCF_Grad_PsiRatio(i,j) + (*JastrowGrad)(i,j);
	
        /*
          We are subtracting (*JastrowGrad)(i,j) because one too many
          was added in the 2*Grad_PsiRatio(i,j) term
          We are including the psuedo-force at this stage
          There are some terms that cancel when calculating
	  
          Laplacian_PsiRatio = del^2 log psi + (del log psi)^2 = 2F^2 - 4T
          KE = -0.5*Laplacian_PsiRatio = 2T - F^2
        */
        Laplacian_PsiRatio += (*JastrowGrad)(i,j) *
	  (2*wd->gradPsiRatio(i,j) - (*JastrowGrad)(i,j));
	
        //This is the 2 * F^2 term = ( del log Psi )^2
        //PsuedoForce2_PsiRatio += wd->gradPsiRatio(i,j) * wd->gradPsiRatio(i,j);
      }

  //This switches to KE = F
  //Laplacian_PsiRatio = -2.0*PsuedoForce2_PsiRatio;
  
  //This switches to KE = T
  //Laplacian_PsiRatio = 0.5*Laplacian_PsiRatio - PsuedoForce2_PsiRatio;

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
      for(int ji=0; ji<Input->JP.getNumberJWParameters(); ji++)
	{
	  wd->rp_a(ai)    = Jastrow.get_p_a_ln(walker,ji);
	  wd->p3_xxa(ai)  = Jastrow.get_p3_xxa_ln(walker,ji);
	  
	  Array2D<double>* j_p_xa = Jastrow.get_p2_xa_ln(walker,ji);
	  
	  for (int i=0; i<Input->WF.getNumberElectrons(); i++)
	    for (int j=0; j<3; j++)
	      {
		wd->p3_xxa(ai) += 2.0 * (*j_p_xa)(i,j) * wd->SCF_Grad_PsiRatio(i,j);
		wd->p3_xxa(ai) += 2.0 * (*j_p_xa)(i,j) * (*JastrowGrad)(i,j);
	      }
	  ai++;
	}
      
      for(int ci=0; ci<Input->WF.getNumberCIParameters(); ci++)
	{
	  wd->rp_a(ai) = termPR(ci) / Input->WF.CI_coeffs(ci);
	  
	  for(int cj=0; cj<Input->WF.getNumberDeterminants(); cj++)
	    {
	      double dtermPsiRatio_ci_dcj;
	      if(cj != ci)
		dtermPsiRatio_ci_dcj = -termPR(ci) * termPR(cj);
	      else
		dtermPsiRatio_ci_dcj = -termPR(ci) * termPR(ci) + termPR(ci);
	      dtermPsiRatio_ci_dcj /= Input->WF.CI_coeffs(ci);

	      if(alphaGrad) term_AlphaGrad = alphaGrad->array()[cj];
	      if(betaGrad)   term_BetaGrad =  betaGrad->array()[cj];
	      
	      if(nalpha > 0)
		wd->p3_xxa(ai)  += dtermPsiRatio_ci_dcj * (*alphaLaplacian)(cj);
	      if(nbeta > 0)
		wd->p3_xxa(ai)  += dtermPsiRatio_ci_dcj * (*betaLaplacian)(cj);

	      for(int i=0; i<nalpha; i++)
		  for(int j=0; j<3; j++)
		    wd->p3_xxa(ai) += 2.0 * dtermPsiRatio_ci_dcj * term_AlphaGrad[i][j] * (*JastrowGrad)(i,j);

	      for(int i=0; i<nbeta; i++)
		for(int j=0; j<3; j++)
		    wd->p3_xxa(ai) += 2.0 * dtermPsiRatio_ci_dcj * term_BetaGrad[i][j]  * (*JastrowGrad)(i+nalpha,j);
	    }
	  ai++;
	}

      int numDet = Input->WF.getNumberDeterminants();
      int unusedIndicator = Input->WF.getUnusedIndicator();
      Array1D<double> dtermPsi_dai(numDet);
      if(ai < numAI)
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
	      dtermPsi_dai = 0.0;
	      double dSCF_Psi_dai = 0.0;

	      for(int ci=0; ci<numDet; ci++)
		{
		  int oa = Input->WF.AlphaOccupation(ci,o);
		  int ob = Input->WF.BetaOccupation(ci,o);

		  double dtermPsi_a = 0.0;
		  if(oa != unusedIndicator)
		    dtermPsi_a += Alpha.get_p_a(walker,ci)->get(oa,bf) * (*betaPsi)(ci);
		  if(ob != unusedIndicator)
		    dtermPsi_a += Beta.get_p_a(walker,ci)->get(ob,bf) * (*alphaPsi)(ci);

		  //This determinant might not depend upon this orbital
		  if(dtermPsi_a == 0) continue;

		  dtermPsi_a *= Input->WF.CI_coeffs(ci);
		  dtermPsi_dai(ci) = dtermPsi_a;
		  dSCF_Psi_dai += dtermPsi_a;
		}

	      wd->rp_a(ai) = dSCF_Psi_dai / SCF_Psi;

	      for(int ci=0; ci<numDet; ci++)
		{
		  int oa = Input->WF.AlphaOccupation(ci,o);
		  int ob = Input->WF.BetaOccupation(ci,o);

		  // d (N/D) = dN / D - N dD / D^2 ; N = termPsiRatio(ci), D = SCF_Psi
		  double dtermPsiRatio_dai;
		  dtermPsiRatio_dai = (dtermPsi_dai(ci) - termPR(ci) * dSCF_Psi_dai) / SCF_Psi;

		  // d ( N * (A + B) ) = dN * (A + B) + N * (dA + dB)
		  if(nalpha > 0)
		    wd->p3_xxa(ai) += dtermPsiRatio_dai * (*alphaLaplacian)(ci);
		  if(nbeta > 0)
		    wd->p3_xxa(ai) += dtermPsiRatio_dai * (*betaLaplacian)(ci);

		  if(oa != unusedIndicator)
		    wd->p3_xxa(ai) += termPR(ci) * Alpha.get_p3_xxa(walker,ci)->get(oa,bf);
		  if(ob != unusedIndicator)
		    wd->p3_xxa(ai) += termPR(ci) * Beta.get_p3_xxa(walker,ci)->get(ob,bf);

		  //And now add in the terms from the gradient
		  double Term_dSlater_a = 0.0;
		  if(oa != unusedIndicator)
		    for(int i=0; i<nalpha; i++)
		      for(int j=0; j<3; j++)
			Term_dSlater_a += Alpha.get_p2_xa(walker,ci,i,j)->get(oa,bf) * (*JastrowGrad)(i,j);
		  if(ob != unusedIndicator)
		    for(int i=0; i<nbeta; i++)
		      for(int j=0; j<3; j++)
			Term_dSlater_a += Beta.get_p2_xa(walker,ci,i,j)->get(ob,bf) * (*JastrowGrad)(i+nalpha,j);
		  Term_dSlater_a *= 2.0 * termPR(ci);

		  if(alphaGrad) term_AlphaGrad = alphaGrad->array()[ci];
		  if(betaGrad)   term_BetaGrad =  betaGrad->array()[ci];

		  double Slater_dTerm_a = 0.0;
		  for(int i=0; i<nalpha; i++)
		    for(int j=0; j<3; j++)
		      Slater_dTerm_a += term_AlphaGrad[i][j] * (*JastrowGrad)(i,j);
		  for(int i=0; i<nbeta; i++)
		    for(int j=0; j<3; j++)
		      Slater_dTerm_a += term_BetaGrad[i][j]  * (*JastrowGrad)(i+nalpha,j);
		  Slater_dTerm_a *= 2.0 * dtermPsiRatio_dai;
		  
		  wd->p3_xxa(ai) += Slater_dTerm_a + Term_dSlater_a;
		}
	      ai++;
	    }
	}

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
      
      //Make sure that the derivative = (f(x + h) - f(x))/h
      //as h -> 0 by editing a ckmf file
      /*
      globalInput.printAISummary();

      //int aStart = Input->getNumberAIParameters() - Input->WF.getNumberORParameters();
      int aStart = 0;
      int aStop  = Input->getNumberAIParameters();
      Array1D<double> params = Input->getAIParameters();
      globalInput.printAIParameters(cout,"Parameters:",25,params,false);
      printf("function %23c p   %25.15e p2_xx  %25.15e\n",
	     ' ',
	     (double)Psi,
	     -0.5*Laplacian_PsiRatio);
      for(int ai=aStart; ai<aStop; ai++)
	printf("ai %3i %25.15e p_a %25.15e p3_xxa %25.15e\n",
	       ai,
	       params(ai),
	       (double)(Psi) * wd->rp_a(ai),
	       wd->p3_xxa(ai));
       exit(0);
      //*/
      termPR.deallocate();
    }

  // Calculate the basis function density
  if (Input->flags.calculate_bf_density == 1)
  {
    wd->chiDensity = 0.0;
    if (nalpha > 0)
      wd->chiDensity = wd->chiDensity + *(Alpha.getChiDensity(walker));
    if (nbeta > 0)
      wd->chiDensity = wd->chiDensity + *(Beta.getChiDensity(walker));
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
      int closest_nucleus_index  =
        Input->Molecule.findClosestNucleusIndex(*x,i);
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

      a = 0.5*(1.0 + closest_nucleus_to_electron_norm_vector *
               electron_velocity_norm_vector) + temp/(10.0*(4.0+temp));

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
      cerr << "WARNING: trial Grad Psi Ratio is NaN "; 
      cerr << "   lengthGradTrialModified = " << lengthGradTrialModified << endl;
      cerr << "   lengthGradTrialUnmodified = " << lengthGradTrialUnmodified << endl;
      wd->isSingular = true;	
    } 
}

void QMCSCFJastrow::calculate_E_Local(int i)
{
  // must calculate Laplacian_PsiRatio before calling this

  E_Local = -0.5 * Laplacian_PsiRatio + PE.getEnergy(i);
}

QMCGreensRatioComponent QMCSCFJastrow::getPsi()
{
  return Psi;
}

double QMCSCFJastrow::getLocalEnergy()
{
  return E_Local;
}

double QMCSCFJastrow::getKineticEnergy()
{
  return -0.5 * Laplacian_PsiRatio;
}

double QMCSCFJastrow::getPotentialEnergy(int i)
{
  return PE.getEnergy(i);
}

double QMCSCFJastrow::getEnergyEE(int i)
{
  return PE.getEnergyEE(i);
}

double QMCSCFJastrow::getEnergyNE(int i)
{
  return PE.getEnergyNE(i);
}

bool QMCSCFJastrow::isSingular(int walker)
{
  if(Alpha.isSingular(walker) ||
     Beta.isSingular(walker))
    return true;
  else
    return false;
}


