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

#include "QMCFunctions.h"

QMCFunctions::QMCFunctions()
{
  nalpha = 0;
  nbeta = 0;
}

QMCFunctions::QMCFunctions(QMCInput *INPUT, QMCHartreeFock *HF)
{
  initialize(INPUT,HF);
}

QMCFunctions::QMCFunctions(QMCInput *INPUT)
{
  QMCHartreeFock* HF = 0;
  initialize(INPUT,HF);
}

QMCFunctions::QMCFunctions(const QMCFunctions & rhs )
{
  *this = rhs;
}

QMCFunctions::~QMCFunctions()
{
  SCF_Grad_PsiRatio.deallocate();
}

void QMCFunctions::operator=(const QMCFunctions & rhs )
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
  SCF_Laplacian_PsiRatio = rhs.SCF_Laplacian_PsiRatio;
  SCF_Grad_PsiRatio      = rhs.SCF_Grad_PsiRatio;
  E_Local                = rhs.E_Local;
}

void QMCFunctions::initialize(QMCInput *INPUT, QMCHartreeFock *HF)
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

  SCF_Grad_PsiRatio.allocate(Input->WF.getNumberElectrons(),3);    
}

void QMCFunctions::evaluate(Array1D<QMCWalkerData *> &walkerData,
                            Array1D<Array2D<double> * > &xData,
                            int num, bool writeConfigs)
{
  //These are some useful benchmarking variables
  static const bool printKE = false;
  static const bool showTimings = !true;
  static double averageT = 0, timeT = 0;
  static double numT = -5;  
  Stopwatch sw = Stopwatch();
  if(showTimings)
  {
      sw.reset(); sw.start();
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
  if (nalpha > 0) Alpha.gpuEvaluate(xData,gpp);
  if (nbeta > 0)  Beta.gpuEvaluate(xData,gpp);
  Jastrow.setUpGPU(Alpha.gpuBF.getElectronicTexture(),
                    Beta.gpuBF.getElectronicTexture(), gpp);
#endif

  /*
    num-gpp is the number the GPU didn't calculate
    gpp is the first index the CPU needs to handle
  */
  //These are pure CPU routines
  Jastrow.evaluate(xData,num-gpp,gpp);
  if (nalpha > 0) Alpha.evaluate(xData,num-gpp,gpp);
  if (nbeta > 0)  Beta.evaluate(xData,num-gpp,gpp);

  PE.evaluate(xData,num);

#ifdef QMC_GPU
  if(true){
    Stopwatch finisher = Stopwatch();
    finisher.reset(); finisher.start();
    glFinish();
    finisher.stop();
    printf("finish time: %15d\n", finisher.timeMS() );
  }
#endif

  //These are CPU if no GPU, otherwise mixed CPU/GPU
  if (nalpha > 0) Alpha.update_Ds();
  if (nbeta > 0)  Beta.update_Ds();
#ifdef QMC_GPU
  Jastrow.gpuEvaluate(xData,gpp);
#endif

  if(showTimings)
  {
    sw.stop();
    timeT = sw.timeMS();
    if(numT >= 0) averageT += sw.timeMS();
    if( num > 1) numT++;
    cout << "total time: " << (int)(timeT+0.5) << " ( " << (int)(averageT/(numT)+0.5) << ")\n";
  }

  for(int i=0; i<num; i++)
    {
      calculate_Psi_quantities(walkerData(i)->gradPsiRatio,
                               walkerData(i)->chiDensity, i);
      calculate_Modified_Grad_PsiRatio(*xData(i),
                                       walkerData(i)->modifiedGradPsiRatio,
                                       walkerData(i)->gradPsiRatio);
      calculate_E_Local(i);

      walkerData(i)->isSingular = isSingular(i);
      walkerData(i)->localEnergy = getLocalEnergy();
      walkerData(i)->kineticEnergy = getKineticEnergy();
      walkerData(i)->potentialEnergy = getPotentialEnergy(i);
      walkerData(i)->neEnergy = getEnergyNE(i);
      walkerData(i)->eeEnergy = getEnergyEE(i);
      walkerData(i)->psi = getPsi();

      if(writeConfigs)
        Input->outputer.writeCorrelatedSamplingConfiguration(*xData(i),
              SCF_Laplacian_PsiRatio,SCF_Grad_PsiRatio,Jastrow.getLnJastrow(i),
							      PE.getEnergy(i));
      if(printKE)
      {
        printf("%4d: ",i);
        printf("KE # %s%18.15e ",(*walkerData(i)).kineticEnergy<0?"":" ",(*walkerData(i)).kineticEnergy );
        printf("TE # %s%18.15g\n",(*walkerData(i)).localEnergy<0?"":" ",(*walkerData(i)).localEnergy );
      }
    }
}

void QMCFunctions::evaluate(Array2D<double> &X, QMCWalkerData & data)
{
  Array1D< Array2D<double>* > temp;
  temp.allocate(1);
  temp(0) = &X;

#ifdef QMC_GPU
  //These are pure GPU
  if (nalpha > 0) Alpha.gpuEvaluate(temp,1);
  if (nbeta > 0)  Beta.gpuEvaluate(temp,1);
#else
  //These are pure CPU
  if (nalpha > 0) Alpha.evaluate(temp,1,0);
  if (nbeta > 0)  Beta.evaluate(temp,1,0);
#endif

  if (nalpha > 0) Alpha.update_Ds();
  if (nbeta > 0)  Beta.update_Ds();
  //We just let the CPU initialize the Jastrow
  Jastrow.evaluate(temp,1,0);
  PE.evaluate(temp,1);

  calculate_Psi_quantities(data.gradPsiRatio,data.chiDensity,0);
  calculate_Modified_Grad_PsiRatio(X,data.modifiedGradPsiRatio,
				   data.gradPsiRatio);
  calculate_E_Local(0);

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
void QMCFunctions::calculate_Psi_quantities(Array2D<double> & Grad_PsiRatio,
    Array1D<double> & Chi_Density, int walker)
{
  QMCGreensRatioComponent SCF_sum = 0.0;
  SCF_Laplacian_PsiRatio = 0.0;
  Laplacian_PsiRatio = 0.0;

  SCF_Grad_PsiRatio = 0.0;
  Grad_PsiRatio = 0.0;

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

  for (int i=0; i<Input->WF.getNumberDeterminants(); i++)
    {
      termPsi(i) = Input->WF.CI_coeffs(i);
      termPsi(i).multiplyBy((*alphaPsi)(i)).multiplyBy((*betaPsi)(i));
      SCF_sum += termPsi(i);
    }

  for (int i=0; i<Input->WF.getNumberDeterminants(); i++)
    {
      double psiRatio = termPsi(i)/SCF_sum;

      if (nalpha > 0)
	{
	  double** term_AlphaGrad = alphaGrad->array()[i];

	  for (int j=0; j<nalpha; j++)
	    for (int k=0; k<3; k++)
	      SCF_Grad_PsiRatio(j,k) += psiRatio*term_AlphaGrad[j][k];
	}

      if (nbeta > 0)
	{
	  double** term_BetaGrad = betaGrad->array()[i];

	  for (int j=0; j<nbeta; j++)
	    for (int k=0; k<3; k++)
	      SCF_Grad_PsiRatio(j+nalpha,k) += psiRatio*term_BetaGrad[j][k];
	}

      SCF_Laplacian_PsiRatio += psiRatio *
                                ((*alphaLaplacian)(i) + (*betaLaplacian)(i));
    }

  QMCGreensRatioComponent jastrowValue =
    QMCGreensRatioComponent(1.0,1.0,0.0,Jastrow.getLnJastrow(walker));

  Psi = SCF_sum * jastrowValue;

  Array2D<double>* JastrowGrad = Jastrow.getGradientLnJastrow(walker);

  for (int i=0; i<Input->WF.getNumberElectrons(); i++)
    for (int j=0; j<3; j++)
      Grad_PsiRatio(i,j) = SCF_Grad_PsiRatio(i,j) + (*JastrowGrad)(i,j);

  Laplacian_PsiRatio = SCF_Laplacian_PsiRatio + 
    Jastrow.getLaplacianLnJastrow(walker);

  for(int i=0; i<Input->WF.getNumberElectrons(); i++)
    for(int j=0; j<3; j++)
      Laplacian_PsiRatio += (*JastrowGrad)(i,j) *
                            (2*Grad_PsiRatio(i,j) - (*JastrowGrad)(i,j));

  // Calculate the basis function density

  if (Input->flags.calculate_bf_density == 1)
    {
      Chi_Density = 0.0;
      if (nalpha > 0)
	Chi_Density = Chi_Density + *(Alpha.getChiDensity(walker));
      if (nbeta > 0)
	Chi_Density = Chi_Density + *(Beta.getChiDensity(walker));
    }
}

void QMCFunctions::calculate_Modified_Grad_PsiRatio(Array2D<double> &X,
    Array2D<double> &Modified_Grad_PsiRatio,
    Array2D<double> &Grad_PsiRatio)
{
  // Call this after calculate_Grad_PsiRatio() is called

  if( Input->flags.QF_modification_type == "none" )
    // do not modify the quantum force
    Modified_Grad_PsiRatio = Grad_PsiRatio;

  else if( Input->flags.QF_modification_type == "umrigar93_equalelectrons" )
    {
      // from Umrigar, Nightingale, and Runge JCP 99(4) 2865; 1993 eq 34.
      // The QF for each electron is changed by the same factor

      double a = Input->flags.umrigar93_equalelectrons_parameter;

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
      for(int i=0; i<Grad_PsiRatio.dim1(); i++)
	for(int j=0; j<Grad_PsiRatio.dim2(); j++)
	  magsqQF += Grad_PsiRatio(i,j) * Grad_PsiRatio(i,j);

      double factor = ( -1.0 + sqrt( 1.0 + 2*a*magsqQF*Input->flags.dt )) /
                      ( a*magsqQF*Input->flags.dt);

      Modified_Grad_PsiRatio = Grad_PsiRatio;
      Modified_Grad_PsiRatio *= factor;
    }
  else if( Input->flags.QF_modification_type == "umrigar93_unequalelectrons" )
    {
      // from Umrigar, Nightingale, and Runge JCP 99(4) 2865; 1993 eq 35-36.
      // The QF for each electron is changed by a different factor

      Modified_Grad_PsiRatio = Grad_PsiRatio;
      Array1D<double> closest_nucleus_to_electron_norm_vector(3);
      Array1D<double> electron_velocity_norm_vector(3);
      for(int i=0; i<Modified_Grad_PsiRatio.dim1(); i++)
        {
          // Find the closest nucleus to electron i
          int closest_nucleus_index  =
            Input->Molecule.findClosestNucleusIndex(X,i);

          double closest_nucleus_distance_squared = 0.0;

          for(int xyz=0; xyz<3; xyz++)
            {
              double temp = X(i,xyz) -
                     Input->Molecule.Atom_Positions(closest_nucleus_index,xyz);

              closest_nucleus_distance_squared += temp * temp;
            }

          // now we have the closest nucleus identified

          // Charge of this nucleus
          int closest_nucleus_Z =
            Input->Molecule.Z(closest_nucleus_index);

          // Unit vector from nearest nucleus to the electron
          for(int xyz=0; xyz<3; xyz++)
            {
              closest_nucleus_to_electron_norm_vector(xyz) =
                X(i,xyz) -
                Input->Molecule.Atom_Positions(closest_nucleus_index,xyz);
            }
          closest_nucleus_to_electron_norm_vector *= 1.0/sqrt(
                closest_nucleus_to_electron_norm_vector *
                closest_nucleus_to_electron_norm_vector);

          // Unit vector in the direction of electron velocity
          for(int xyz=0; xyz<3; xyz++)
            {
              electron_velocity_norm_vector(xyz) =
                Grad_PsiRatio(i,xyz);
            }
          double electron_velocity_norm_squared =
            electron_velocity_norm_vector * electron_velocity_norm_vector;
          electron_velocity_norm_vector *= 1.0/sqrt(
                                             electron_velocity_norm_squared );

          // Now calculate the a(r) for the equation
          double temp = closest_nucleus_Z * closest_nucleus_Z *
                        closest_nucleus_distance_squared;

          double a = 0.5*(1.0 + closest_nucleus_to_electron_norm_vector *
                          electron_velocity_norm_vector) +
                     temp/(10.0*(4.0+temp));

          // Now calculate the factor the QF of this electron is modified by
          double factor = (-1.0 + sqrt(1.0+
                                       2.0*a*electron_velocity_norm_squared*Input->flags.dt))/
                          (a*electron_velocity_norm_squared*Input->flags.dt);

          // Now adjust the QF of electron I
          for(int xyz=0; xyz<3; xyz++)
            {
              Modified_Grad_PsiRatio(i,xyz) *= factor;
            }
        }
    }
  else
    {
      cerr << "ERROR: Unknown value for QF_modification_type set!" << endl;
      exit(0);
    }
}

void QMCFunctions::calculate_E_Local(int i)
{
  // must calculate Laplacian_PsiRatio before calling this

  E_Local = -0.5 * Laplacian_PsiRatio + PE.getEnergy(i);
}

QMCGreensRatioComponent QMCFunctions::getPsi()
{
  return Psi;
}

double QMCFunctions::getLocalEnergy()
{
  return E_Local;
}

double QMCFunctions::getKineticEnergy()
{
  return -0.5 * Laplacian_PsiRatio;
}

double QMCFunctions::getPotentialEnergy(int i)
{
  return PE.getEnergy(i);
}

double QMCFunctions::getEnergyEE(int i)
{
  return PE.getEnergyEE(i);
}

double QMCFunctions::getEnergyNE(int i)
{
  return PE.getEnergyNE(i);
}

bool QMCFunctions::isSingular(int walker)
{
  if (nalpha > 0 && Alpha.isSingular(walker))
    return true;
  else if (nbeta > 0 && Beta.isSingular(walker))
    return true;
  else
    return false;
}





