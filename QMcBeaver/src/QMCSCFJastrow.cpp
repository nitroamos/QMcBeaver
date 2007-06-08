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

  if (nalpha > 0)
    Alpha.initialize(Input,0,Input->WF.getNumberAlphaElectrons()-1,
                     &Input->WF.AlphaOccupation,&Input->WF.AlphaCoeffs);

  if (nbeta > 0)
    Beta.initialize(Input,Input->WF.getNumberAlphaElectrons(),
                    Input->WF.getNumberElectrons()-1,&Input->WF.BetaOccupation,
                    &Input->WF.BetaCoeffs);

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
  if (nalpha > 0)
    Alpha.gpuEvaluate(xData,gpp);
  if (nbeta > 0)
    Beta.gpuEvaluate(xData,gpp);
  Jastrow.setUpGPU(Alpha.gpuBF.getElectronicTexture(),
                   Beta.gpuBF.getElectronicTexture(), gpp);
#endif

  /*
    num-gpp is the number the GPU didn't calculate
    gpp is the first index the CPU needs to handle
  */
  //These are pure CPU routines
  Jastrow.evaluate(xData,num-gpp,gpp);
  if (nalpha > 0)
    Alpha.evaluate(xData,num-gpp,gpp);
  if (nbeta > 0)
    Beta.evaluate(xData,num-gpp,gpp);

  PE.evaluate(xData,num);

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
  if (nalpha > 0)
    Alpha.update_Ds();
  if (nbeta > 0)
    Beta.update_Ds();
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
  Array1D< Array2D<double>* > temp;
  temp.allocate(1);
  temp(0) = &X;

#ifdef QMC_GPU
  //These are pure GPU
  if (nalpha > 0)
    Alpha.gpuEvaluate(temp,1);
  if (nbeta > 0)
    Beta.gpuEvaluate(temp,1);
#else
  //These are pure CPU
  if (nalpha > 0)
    Alpha.evaluate(temp,1,0);
  if (nbeta > 0)
    Beta.evaluate(temp,1,0);
#endif

  if (nalpha > 0)
    Alpha.update_Ds();
  if (nbeta > 0)
    Beta.update_Ds();
  //We just let the CPU initialize the Jastrow
  Jastrow.evaluate(temp,1,0);
  PE.evaluate(temp,1);

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

//  We end up doing a division with SCF_sum. Therefore, it seems to be critical
//  to make sure we aren't dividing by too small a number
void QMCSCFJastrow::calculate_Psi_quantities(int walker)
{
  QMCGreensRatioComponent SCF_sum = 0.0;
  Laplacian_PsiRatio = 0.0;
  //double PsuedoForce2_PsiRatio = 0;

  wd->SCF_Laplacian_PsiRatio = 0.0;
  wd->SCF_Grad_PsiRatio      = 0.0;
  wd->gradPsiRatio           = 0.0;

  Array1D<QMCGreensRatioComponent> termPsi;
  termPsi.allocate(Input->WF.getNumberDeterminants());

  Array1D<double>* alphaPsi = 0;
  Array1D<double>* alphaLaplacian = 0;
  Array3D<double>* alphaGrad = 0;

  Array1D<double> tempPsi;
  Array1D<double> tempLaplacian;

  if (nalpha > 0)
  {
    alphaPsi       = Alpha.getPsi(walker);
    alphaLaplacian = Alpha.getLaplacianPsiRatio(walker);
    alphaGrad      = Alpha.getGradPsiRatio(walker);
  }
  else
  {
    tempPsi.allocate(Input->WF.getNumberDeterminants());
    tempPsi = 1.0;
    alphaPsi = &tempPsi;

    tempLaplacian.allocate(Input->WF.getNumberDeterminants());
    tempLaplacian = 0.0;
    alphaLaplacian = &tempLaplacian;
  }

  Array1D<double>* betaPsi = 0;
  Array1D<double>* betaLaplacian = 0;
  Array3D<double>* betaGrad = 0;

  if (nbeta > 0)
  {
    betaPsi       = Beta.getPsi(walker);
    betaLaplacian = Beta.getLaplacianPsiRatio(walker);
    betaGrad      = Beta.getGradPsiRatio(walker);
  }
  else
  {
    tempPsi.allocate(Input->WF.getNumberDeterminants());
    tempPsi = 1.0;
    betaPsi = &tempPsi;

    tempLaplacian.allocate(Input->WF.getNumberDeterminants());
    tempLaplacian = 0.0;
    betaLaplacian = &tempLaplacian;
  }

  double tempAlphaPsi = 0.0;
  double tempBetaPsi = 0.0;

  for (int i=0; i<Input->WF.getNumberDeterminants(); i++)
  {
    // If a determinant is NaN we exclude it from the sum.
    termPsi(i) = Input->WF.CI_coeffs(i);

    tempAlphaPsi = (*alphaPsi)(i);
    if (IeeeMath::isNaN(tempAlphaPsi))
      termPsi(i).multiplyBy(0.0);
    else
      termPsi(i).multiplyBy(tempAlphaPsi);

    tempBetaPsi = (*betaPsi)(i);
    if (IeeeMath::isNaN(tempBetaPsi))
      termPsi(i).multiplyBy(0.0);
    else
      termPsi(i).multiplyBy(tempBetaPsi);

    SCF_sum += termPsi(i);
  }

  double psiRatio = 0.0;

  for (int i=0; i<Input->WF.getNumberDeterminants(); i++)
  {
    tempAlphaPsi = (*alphaPsi)(i);
    tempBetaPsi = (*betaPsi)(i);

    if (IeeeMath::isNaN(tempBetaPsi) || IeeeMath::isNaN(tempAlphaPsi))
    {
      // If a determinant is NaN we exclude it from the gradient and
      // laplacian sums as well.
    }
    else
    {
      psiRatio = termPsi(i)/SCF_sum;

      if (nalpha > 0)
      {
        double** term_AlphaGrad = alphaGrad->array()[i];

        for (int j=0; j<nalpha; j++)
          for (int k=0; k<3; k++)
            wd->SCF_Grad_PsiRatio(j,k) += psiRatio*term_AlphaGrad[j][k];

	wd->SCF_Laplacian_PsiRatio += psiRatio * (*alphaLaplacian)(i);
      }

      if (nbeta > 0)
      {
        double** term_BetaGrad = betaGrad->array()[i];

        for (int j=0; j<nbeta; j++)
          for (int k=0; k<3; k++)
            wd->SCF_Grad_PsiRatio(j+nalpha,k)+=psiRatio*term_BetaGrad[j][k];

	wd->SCF_Laplacian_PsiRatio += psiRatio * (*betaLaplacian)(i);
      }
    }
  }
  
  QMCGreensRatioComponent jastrowValue =
    QMCGreensRatioComponent(1.0,1.0,0.0,Jastrow.getLnJastrow(walker));
  
  Psi = SCF_sum * jastrowValue;
  
  Laplacian_PsiRatio = wd->SCF_Laplacian_PsiRatio +
    Jastrow.getLaplacianLnJastrow(walker);

  Array2D<double>* JastrowGrad = Jastrow.getGradientLnJastrow(walker);

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
      wd->p3_xxa = 0.0;  
      for(int ai=0; ai<Jastrow.getNumAI(); ai++)
	{
	  //we never need dPsi/dai by itself, we only need the ratio. So this
	  //is dPsi/dai/Psi.
	  wd->rp_a(ai)    = Jastrow.get_p_a_ln(walker,ai);
	  wd->p3_xxa(ai)  = Jastrow.get_p3_xxa_ln(walker,ai);
	  
	  Array2D<double>* j_p_xa = Jastrow.get_p2_xa_ln(walker,ai);
	  
	  for (int i=0; i<Input->WF.getNumberElectrons(); i++)
	    {
	      for (int j=0; j<3; j++)
		{
		  wd->p3_xxa(ai) += 2 * (*j_p_xa)(i,j) * wd->SCF_Grad_PsiRatio(i,j);
		  wd->p3_xxa(ai) += 2 * (*j_p_xa)(i,j) * (*JastrowGrad)(i,j);
		}
	    }
	  wd->p3_xxa(ai) = -0.5*wd->p3_xxa(ai);
	}
      
      int shift = Jastrow.getNumAI();
      for(int ai=0; ai<Input->WF.getNumberDeterminants()-1; ai++)
	{
	  int ci = ai+1;
	  psiRatio = termPsi(ci)/SCF_sum;
	  double dPsiRatio_dci = (psiRatio - psiRatio * psiRatio) / Input->WF.CI_coeffs(ci);
	  
	  wd->rp_a(ai+shift)    = psiRatio / Input->WF.CI_coeffs(ci);

	  for(int cj=0; cj<Input->WF.getNumberDeterminants(); cj++)
	    {
	      double dPsiRatio_dcj;
	      if(cj != ci)
		dPsiRatio_dcj = -termPsi(cj) * psiRatio / SCF_sum / Input->WF.CI_coeffs(ci);
	      else
		dPsiRatio_dcj = dPsiRatio_dci;

	      double** term_AlphaGrad = 0;
	      double** term_BetaGrad = 0;
	      if(nalpha > 0) term_AlphaGrad = alphaGrad->array()[cj];
	      if(nbeta > 0)  term_BetaGrad = betaGrad->array()[cj];
	      
	      wd->p3_xxa(ai+shift)  += dPsiRatio_dcj * ((*alphaLaplacian)(cj) + (*betaLaplacian)(cj));
	      
	      for (int i=0; i<nalpha; i++)
		for (int j=0; j<3; j++)
		  wd->p3_xxa(ai+shift) += 2 * (*JastrowGrad)(i,j) * term_AlphaGrad[i][j] * dPsiRatio_dcj;
	      for (int i=0; i<nbeta; i++)
		for (int j=0; j<3; j++)
		  wd->p3_xxa(ai+shift) += 2 * (*JastrowGrad)(i+nalpha,j) * term_BetaGrad[i][j] * dPsiRatio_dcj;
	    }
	  wd->p3_xxa(ai+shift) = -0.5*wd->p3_xxa(ai+shift);
	}

      //This is how you hand check the derivatives
      /*
	printf("actual P=%25.15f  XX=%25.15f\n",(double)Psi,-0.5*Laplacian_PsiRatio);
	for(int ai=0; ai<wd->p3_xxa.dim1(); ai++)
	printf("ai=%2i RP=%25.15f XXA=%25.15f\n",ai,(double)(Psi * wd->rp_a(ai)),wd->p3_xxa(ai));
	exit(0);
      //*/
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

  termPsi.deallocate();
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
  if (nalpha > 0 && Alpha.isSingular(walker))
    return true;
  else if (nbeta > 0 && Beta.isSingular(walker))
    return true;
  else
    return false;
}

int QMCSCFJastrow::getNumAI()
{
  /*
    Since the wavefunction is invariant to normalization,
    there are only N-1 independent CI coefficients.
  */
  int numAI = Input->WF.getNumberDeterminants() - 1;

  if(Input->flags.calculate_Derivatives == 1)
    numAI += Jastrow.getNumAI();
  else
    numAI = 0;
  return numAI;
}



