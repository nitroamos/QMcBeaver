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
{}

QMCFunctions::QMCFunctions(QMCInput *INPUT)
{
  initialize(INPUT);
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

void QMCFunctions::initialize(QMCInput *INPUT)
{
  Input = INPUT;

  Alpha.initialize(Input,0,Input->WF.getNumberAlphaElectrons()-1,
		   &Input->WF.AlphaOccupation,&Input->WF.AlphaCoeffs);
  Beta.initialize(Input,Input->WF.getNumberAlphaElectrons(),
		  Input->WF.getNumberElectrons()-1,&Input->WF.BetaOccupation,
		  &Input->WF.BetaCoeffs);

  PE.initialize(Input);
  Jastrow.initialize(Input);

  SCF_Grad_PsiRatio.allocate(Input->WF.getNumberElectrons(),3);
}

void QMCFunctions::evaluate(Array1D<QMCWalkerData *> &walkerData,
                            Array1D<Array2D<double> * > &xData,
                            int num, bool writeConfigs)
{
  Alpha.update_Ds(xData,num);
  Beta.update_Ds(xData,num);
  Jastrow.evaluate(xData,num);
  PE.evaluate(xData,num);
  Alpha.evaluate(num);
  Beta.evaluate(num);

  for(int i=0; i<num; i++)
    {
      calculate_Psi_quantities((*walkerData(i)).gradPsiRatio,
                               (*walkerData(i)).chiDensity, i);
      calculate_Modified_Grad_PsiRatio(*xData(i),
                                       (*walkerData(i)).modifiedGradPsiRatio,
                                       (*walkerData(i)).gradPsiRatio);
      calculate_E_Local(i);

      (*walkerData(i)).isSingular = isSingular(i);
      (*walkerData(i)).localEnergy = getLocalEnergy();
      (*walkerData(i)).kineticEnergy = getKineticEnergy();
      (*walkerData(i)).potentialEnergy = getPotentialEnergy(i);
      (*walkerData(i)).psi = getPsi();

      if(writeConfigs)
        writeCorrelatedSamplingConfiguration(*(walkerData(i)->configOutput),i);
    }
}

void QMCFunctions::evaluate(Array2D<double> &X, QMCWalkerData & data)
{
  Array1D< Array2D<double>* > temp;
  temp.allocate(1);
  temp(0) = &X;

  Alpha.update_Ds(temp,1);
  Beta.update_Ds(temp,1);
  Alpha.evaluate(1);
  Beta.evaluate(1);

  Jastrow.evaluate(temp,1);
  PE.evaluate(temp,1);

  calculate_Psi_quantities(data.gradPsiRatio,data.chiDensity,0);
  calculate_Modified_Grad_PsiRatio(X, data.modifiedGradPsiRatio, data.gradPsiRatio);
  calculate_E_Local(0);

  data.isSingular = isSingular(0);
  data.localEnergy = getLocalEnergy();
  data.kineticEnergy = getKineticEnergy();
  data.potentialEnergy = getPotentialEnergy(0);
  data.psi = getPsi();
}

/*  we end up doing a division with SCF_sum. Therefore, it seems to be critical
    to make sure we aren't dividing by too small a number*/
void QMCFunctions::calculate_Psi_quantities(Array2D<double> & Grad_PsiRatio,
    Array1D<double> & Chi_Density, int walker)
{
  QMCGreensRatioComponent SCF_sum = 0.0;
  SCF_Laplacian_PsiRatio = 0.0;
  Laplacian_PsiRatio = 0.0;

  Array1D<QMCGreensRatioComponent> termPsi;
  termPsi.allocate(Input->WF.getNumberDeterminants());

  Array1D<double>* alphaPsi = Alpha.getPsi(walker);
  Array1D<double>* alphaLaplacian = Alpha.getLaplacianPsiRatio(walker);
  Array3D<double>* alphaGrad = Alpha.getGradPsiRatio(walker);

  Array1D<double>* betaPsi = Beta.getPsi(walker);
  Array1D<double>* betaLaplacian = Beta.getLaplacianPsiRatio(walker);
  Array3D<double>* betaGrad = Beta.getGradPsiRatio(walker);

  for (int i=0; i<Input->WF.getNumberElectrons(); i++)
    for (int j=0; j<3; j++)
      {
        SCF_Grad_PsiRatio(i,j) = 0.0;
        Grad_PsiRatio(i,j) = 0.0;
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
      double** term_AlphaGrad = (*alphaGrad).array()[i];
      double** term_BetaGrad  = (*betaGrad).array()[i];

      for (int j=0; j<Input->WF.getNumberAlphaElectrons(); j++)
        for (int k=0; k<3; k++)
          SCF_Grad_PsiRatio(j,k) += psiRatio*term_AlphaGrad[j][k];

      for (int j=0; j<Input->WF.getNumberBetaElectrons(); j++)
        for (int k=0; k<3; k++)
          SCF_Grad_PsiRatio(Input->WF.getNumberAlphaElectrons()+j,k) +=
            psiRatio*term_BetaGrad[j][k];

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

  Laplacian_PsiRatio = SCF_Laplacian_PsiRatio+Jastrow.getLaplacianLnJastrow(walker);
  for(int i=0; i<Input->WF.getNumberElectrons(); i++)
    for(int j=0; j<3; j++)
      Laplacian_PsiRatio += (*JastrowGrad)(i,j) *
                            (2*Grad_PsiRatio(i,j) - (*JastrowGrad)(i,j));

  // Calculate the basis function density

  if (Input->flags.calculate_bf_density == 1)
    Chi_Density = *(Alpha.getChiDensity(walker)) + *(Beta.getChiDensity(walker));
}

void QMCFunctions::calculate_Modified_Grad_PsiRatio(Array2D<double> &X,
    Array2D<double> &Modified_Grad_PsiRatio,
    Array2D<double> &Grad_PsiRatio)
{
  // Call this after calculate_Grad_PsiRatio() is called

  if( Input->flags.QF_modification_type == "none" )
    {
      // do not modify the quantum force
      Modified_Grad_PsiRatio = Grad_PsiRatio;
    }
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
        {
          for(int j=0; j<Grad_PsiRatio.dim2(); j++)
            {
              magsqQF += Grad_PsiRatio(i,j) * Grad_PsiRatio(i,j);
            }
        }

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

void QMCFunctions::writeCorrelatedSamplingConfiguration(stringstream& strm, int which)
{
  strm.str("");
  // Print the sum of the laplacians from the alpha and beta electrons
  strm << "D1" << endl;
  strm << "\t" << SCF_Laplacian_PsiRatio << endl;

  // print the GradPsiRatio excluding the Jastrow
  strm << "D2" << endl;

  for(int i=0; i<Input->WF.getNumberElectrons(); i++)
    {
      for(int j=0; j<3; j++)
        {
          strm << "\t" << SCF_Grad_PsiRatio(i,j);
        }
      strm << endl;
    }

  strm << "lnJ\t" << endl;
  strm << "\t" << Jastrow.getLnJastrow(which) << endl;
  strm << "PE\t" << endl;
  strm << "\t" << PE.getEnergy(which) << endl;
}

bool QMCFunctions::isSingular(int walker)
{
  if( Alpha.isSingular(walker) || Beta.isSingular(walker) )
    {
      return true;
    }
  else return false;
}





