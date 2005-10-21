#ifndef QMCHartreeFock_H
#define QMCHartreeFock_H

#include "Array2D.h"
#include "QMCPsiPotential.h"
#include "IeeeMath.h"
#include <math.h>
#include <stdio.h>

class QMCHartreeFock
{
 private:
  // number of electrons and max samples to accumulate in circular queue at a 
  // time.
  static int numalphas, numbetas, numelecs;
  static int sample, maxsamples, numsamples;

  // stores positions and weights of all samples of all electrons
  static Array2D<double> elec_x, elec_y, elec_z, elec_weight;

  QMCPsiPotential PsiPotential;

 public:
  QMCHartreeFock();
  virtual ~QMCHartreeFock();

  void Initialize(QMCInput* IN);
  void AddElectron(int elec,double weight,double x,double y,double z);
  void IncrementSample();
  double GetVEff(int elec, double x, double y, double z);
};

#endif
