// JTS 9-8-05, modification to QMcBeaver
// Hartree-Fock calculations using QMC, aimed at producing basis-independent 
// reference wavefunctions.
// Computes V_eff for each electron by averaging over multiple walkers and 
// configurations.
//
// Julius Su, jsu@caltech.edu

#include "QMCHartreeFock.h"

// define static variables
int QMCHartreeFock::numelecs; 
int QMCHartreeFock::numalphas; 
int QMCHartreeFock::numbetas;
int QMCHartreeFock::sample; 
int QMCHartreeFock::maxsamples; 
int QMCHartreeFock::numsamples;

// stores positions and weights of all samples of all electrons
Array2D<double> QMCHartreeFock::elec_x; 
Array2D<double> QMCHartreeFock::elec_y; 
Array2D<double> QMCHartreeFock::elec_z; 
Array2D<double> QMCHartreeFock::elec_weight;

QMCHartreeFock::QMCHartreeFock()
{
}

QMCHartreeFock::~QMCHartreeFock()
{
  elec_x.deallocate();
  elec_y.deallocate();
  elec_z.deallocate();
  elec_weight.deallocate();
}

void QMCHartreeFock::Initialize(QMCInput* IN)
{
  PsiPotential.Initialize(IN);

  numalphas = IN->WF.getNumberAlphaElectrons();
  numbetas = IN->WF.getNumberBetaElectrons();
  numelecs = numalphas + numbetas;
  maxsamples = IN->flags.hf_num_average;
  sample = 0; // current sample
  numsamples = 0; // number of samples in current queue

  // allocate electron arrays
  elec_x.allocate(numelecs, maxsamples);
  elec_y.allocate(numelecs, maxsamples);
  elec_z.allocate(numelecs, maxsamples);
  elec_weight.allocate(numelecs, maxsamples);
}

void QMCHartreeFock::AddElectron(int elec, double weight, double x, double y, 
				 double z)
{
  // add electron to array
  elec_x(elec, sample) = x;
  elec_y(elec, sample) = y;
  elec_z(elec, sample) = z;
  elec_weight(elec, sample) = weight;
}

void QMCHartreeFock::IncrementSample()
{
  // iterate for circular queue
  sample++;
  if (sample >= maxsamples) sample = 0;
  numsamples++;
  if (numsamples >= maxsamples) numsamples = maxsamples;
}

double QMCHartreeFock::GetVEff(int elec, double x, double y, double z)
{
  // Calculates V_eff by averaging over all other electrons in all samples.
  // This is the mixed estimator of 1/r.
  double dx, dy, dz, r;
  double sum1, sum2, sumV;

  sum1 = 0;
  sumV = 0;

  int noweights = 0;
  for (int i = 0; i < numelecs; i++)
  {
    if (i != elec)
    {
      sum2 = 0;
      double weights = 0;
      for (int j = 0; j < numsamples; j++)
      {
        dx = x - elec_x(i, j);
        dy = y - elec_y(i, j);
        dz = z - elec_z(i, j);
        r = sqrt(dx * dx + dy * dy + dz * dz);
        sum2 += elec_weight(i, j) / r;
        weights += elec_weight(i, j);
      }
      sum1 += sum2 / weights;
      if (weights == 0) noweights = 1;
      // electrons are ordered as follows -- a1 a2 .. an b1 b2 .. bm
      // select alpha or beta wavefunction as appropriate
      sumV += (i < numalphas) ? PsiPotential.GetPoten(i, 0, x, y, z) : PsiPotential.GetPoten(i - numalphas, 1, x, y, z);
    }
  }

  // Calculate the trial estimator of 1/r and extrapolate to the pure estimator
  // using 
  // pure = 2 * trial - mixed.
  double V_mixed = sum1 / 2.0;
  double V_trial = sumV / 2.0;
  double V_pure = 2 * V_mixed - V_trial;
  if (noweights) V_pure = V_trial;
  if (numsamples == 0) return V_trial; else return V_pure;
}
