#ifndef QMCPsiPotential_H
#define QMCPsiPotential_H

#include "QMCInput.h"
#include "MathFunctions.h"

class JTS_ORBITAL
{
public:
  int spin; double occ;
  int numgaussians;
  int *l, *m, *n;
  double *x, *y, *z;
  double *coeff, *coeff0, *exponent; 

  // coeff0 is the primative orbital coefficient
  // coeff is coeff0 times the wavefunction coefficients

  JTS_ORBITAL() {}
  virtual ~JTS_ORBITAL() {};
};

class QMCPsiPotential
{
 private:
  static QMCInput *qmc_input;
  static QMCBasisFunction *qmc_basis;
  static JTS_ORBITAL orbital;

  static double NormGaussian(int l, int m, int n, double exponent);
  static double OverlapVGaussians(int l1, int m1, int n1, double exponent1,
				  double x1, double y1, double z1, int l2,
				  int m2, int n2, double exponent2, double x2,
				  double y2, double z2, double xnuc,
				  double ynuc, double znuc);
  static double DFactorial(int n);
  static double intpow(double x, int n);
  static void binomial_f(int l, int m, double a, double b, double *f);
  static void g_function(int suml, double *f, double ig, double p, double *g);
  static int CountGaussians();
  static int AllocateOrbitals();
  static void LoadOrbitals(int numgaussians);
  static void UseWavefunction(int wfnum, int spin);

 public:
  QMCPsiPotential();
  virtual ~QMCPsiPotential();

  static void Initialize(QMCInput *s_qmc_input);
  static double GetPoten(int orbital, int spin, double x, double y, double z);
};

#endif

