// Compute potential energy from the trial wavefunction density.
// Used to extrapolate to the pure estimator in QMCHartreeFock of (1/r).
// Julius Su, jsu@caltech.edu

#include "QMCPsiPotential.h"

QMCInput *QMCPsiPotential::qmc_input;
QMCBasisFunction *QMCPsiPotential::qmc_basis;
JTS_ORBITAL QMCPsiPotential::orbital;

QMCPsiPotential::QMCPsiPotential()
{
}

QMCPsiPotential::~QMCPsiPotential()
{
}

void QMCPsiPotential::Initialize(QMCInput *s_qmc_input)
{
  qmc_input = s_qmc_input;
  qmc_basis = &(qmc_input->BF);
  LoadOrbitals(AllocateOrbitals());
}

int QMCPsiPotential::CountGaussians()
{
  int numgaussians = 0;
  for (int i = 0; i < qmc_basis->BFCoeffs.dim1(); i++)
    for (int j = 0; j < qmc_basis->BFCoeffs(i).getNumberBasisFunctions(); j++)
      numgaussians += qmc_basis->BFCoeffs(i).N_Gauss(j);
  return numgaussians;
}

int QMCPsiPotential::AllocateOrbitals()
{
  int numgaussians = CountGaussians();

  // Allocate gaussian primatives for the orbitals 
  orbital.l = (int *) malloc(sizeof(int) * numgaussians);
  orbital.m = (int *) malloc(sizeof(int) * numgaussians);
  orbital.n = (int *) malloc(sizeof(int) * numgaussians);
  orbital.x = (double *) malloc(sizeof(double) * numgaussians);
  orbital.y = (double *) malloc(sizeof(double) * numgaussians);
  orbital.z = (double *) malloc(sizeof(double) * numgaussians);
  orbital.coeff = (double *) malloc(sizeof(double) * numgaussians);
  orbital.coeff0 = (double *) malloc(sizeof(double) * numgaussians);
  orbital.exponent = (double *) malloc(sizeof(double) * numgaussians);

  return numgaussians;
}

void QMCPsiPotential::LoadOrbitals(int numgaussians)
{
  int idx = 0;
  for (int i = 0; i < qmc_basis->BFCoeffs.dim1(); i++)
  {
    double x, y, z;
    x = qmc_input->Molecule.Atom_Positions(i, 0);
    y = qmc_input->Molecule.Atom_Positions(i, 1);
    z = qmc_input->Molecule.Atom_Positions(i, 2);
    // loop over number of orbitals
    for (int j = 0; j < qmc_basis->BFCoeffs(i).getNumberBasisFunctions(); j++)
    {
      // loop over number of gaussians (1 for s, 3 for p, 5 for d, etc.)
      for (int k = 0; k < qmc_basis->BFCoeffs(i).N_Gauss(j); k++)
      {
       orbital.x[idx] = x;
       orbital.y[idx] = y;
       orbital.z[idx] = z;
       orbital.exponent[idx] = qmc_basis->BFCoeffs(i).Coeffs(j, k, 0);
       orbital.coeff0[idx] = qmc_basis->BFCoeffs(i).Coeffs(j, k, 1);
       orbital.l[idx] = qmc_basis->BFCoeffs(i).xyz_powers(j, 0);
       orbital.m[idx] = qmc_basis->BFCoeffs(i).xyz_powers(j, 1);
       orbital.n[idx] = qmc_basis->BFCoeffs(i).xyz_powers(j, 2);
       orbital.numgaussians = numgaussians;
       idx++;
      }
    }
  }
}

void QMCPsiPotential::UseWavefunction(int wfnum, int spin)
{
  // Apply wavefunction coefficients to current orbital set, multiplying 
  // coeff0 to get coeff.
  int idx = 0;

  // loop over number of atoms
  for (int i = 0; i < qmc_basis->BFCoeffs.dim1(); i++)
  {
    // loop over number of orbitals
    for (int j = 0; j < qmc_basis->BFCoeffs(i).getNumberBasisFunctions(); j++)
    {
      double factor = (spin == 0 ? qmc_input->WF.AlphaCoeffs(wfnum, j) : qmc_input->WF.BetaCoeffs(wfnum, j));
      // loop over number of gaussians (1 for s, 3 for p, 5 for d, etc.)
      for (int k = 0; k < qmc_basis->BFCoeffs(i).N_Gauss(j); k++)
      {
       orbital.coeff[idx] = orbital.coeff0[idx] * factor;
       idx++;
      }
    }
  }
}

double QMCPsiPotential::GetPoten(int wfnum, int spin, double x, double y, double z)
{
  int j, k;
  int numgaussians = orbital.numgaussians;
  
  // Load in desired wavefunction
  UseWavefunction(wfnum, spin);
 
  // Calculate potential
  double sum = 0;
  for (j = 0; j < numgaussians; j++)
  {
      int l1, m1, n1;
      double exponent1, coeff1;
      double x1, y1, z1;

      x1 = orbital.x[j]; y1 = orbital.y[j]; z1 = orbital.z[j];
      l1 = orbital.l[j]; m1 = orbital.m[j]; n1 = orbital.n[j];
      exponent1 = orbital.exponent[j];
      coeff1 = orbital.coeff[j];

      for (k = 0; k < numgaussians; k++)
      {
        int l2, m2, n2;
        double exponent2, coeff2;
        double x2, y2, z2;
      
        x2 = orbital.x[k]; y2 = orbital.y[k]; z2 = orbital.z[k];
        l2 = orbital.l[k]; m2 = orbital.m[k]; n2 = orbital.n[k];
        exponent2 = orbital.exponent[k];
        coeff2 = orbital.coeff[k];

        //printf("%f %f %i %i %i %f %f %f %f %i %i %i %f %f %f %f %f %f %f\n", coeff1, coeff2, l1, m1, n1, exponent1, x1, y1, z1, l2, m2, n2, exponent2, x2, y2, z2, x, y, z);
        sum += coeff1 * coeff2 * OverlapVGaussians(l1, m1, n1, exponent1, x1, y1, z1, l2, m2, n2, exponent2, x2, y2, z2, x, y, z);
      }
  }
  return sum;
}

double QMCPsiPotential::NormGaussian(int l, int m, int n, double exponent)
{
  return pow(2.0 * exponent / 3.1415926535, 0.75) * sqrt(pow(4.0 * exponent, l + m + n) / (DFactorial(2 * l - 1) * DFactorial(2 * m - 1) * DFactorial(2 * n - 1)));
}

double QMCPsiPotential::OverlapVGaussians(int l1, int m1, int n1, double exponent1, double x1, double y1, double z1, int l2, int m2, int n2, double exponent2, double x2, double y2, double z2, double xnuc, double ynuc, double znuc)
{
  /* Formula from equation 3.1, J. Phys. Soc. Japan 21(11) (1966) 2313-2324 
     Computes Integrate[gaussian1 gaussian2 1/rc] 
  */

  double gamma = exponent1 + exponent2;

  /* Find PA, PB vectors connecting A,B to the product gaussian center */
  double px, py, pz;
  double pax, pay, paz;
  double pbx, pby, pbz;

  px = (exponent1 * x1 + exponent2 * x2) / gamma;
  py = (exponent1 * y1 + exponent2 * y2) / gamma;
  pz = (exponent1 * z1 + exponent2 * z2) / gamma;
  pax = x1 - px; pay = y1 - py; paz = z1 - pz;
  pbx = x2 - px; pby = y2 - py; pbz = z2 - pz;

  int i;
  double f[9];

  double sum = 0;
  
  /* Vector connecting product gaussian center to nucleus */
  double pnx, pny, pnz, p2, p;
  pnx = xnuc - px;
  pny = ynuc - py;
  pnz = znuc - pz;
  p2 = pnx * pnx + pny * pny + pnz * pnz;
  p = sqrt(p2);

  /* Evaluate g functions over allowed ranges  */
  double gx[9], gy[9], gz[9];
  binomial_f(l1, l2, pax, pbx, f);
  g_function(l1 + l2, f, 1.0 / gamma, p, gx);
  binomial_f(m1, m2, pay, pby, f);
  g_function(m1 + m2, f, 1.0 / gamma, p, gy);
  binomial_f(n1, n2, paz, pbz, f);
  g_function(n1 + n2, f, 1.0 / gamma, p, gz);

  /* Evaluate auxiliary function over allowed range */
  double aux[27];
  for (i = 0; i <= l1 + l2 + m1 + m2 + n1 + n2; i++)
     aux[i] = MathFunctions::F_gamma((double) i, gamma * p2);

  /* Perform summation */
  int l, m, n;
   for (l = 0; l <= l1 + l2; l++)
     for (m = 0; m <= m1 + m2; m++)
       for (n = 0; n <= n1 + n2; n++)
         sum += gx[l] * gy[m] * gz[n] * aux[l + m + n];

  /* Calculate constant prefactor and return value */
  double dx, dy, dz;
  dx = x2 - x1; dy = y2 - y1; dz = z2 - z1;
  
  return (2.0 * 3.1415926535 / gamma) * exp(-exponent1 * exponent2 * (dx * dx + dy * dy + dz * dz) / gamma) * sum;
}

double QMCPsiPotential::DFactorial(int n)
{
  /* Double factorial function */
  /* n * (n - 2) * (n - 4) etc. */
  double x = 1;
  while (n > 0)
  {
    x *= n;
    n -= 2;
  }
  return (double) x;
}

double QMCPsiPotential::intpow(double x, int n)
{
  double y = 1;
  while (n > 0)
  {
    y *= x; n--;
  }
  return y;
}

void QMCPsiPotential::binomial_f(int l, int m, double a, double b, double *f)
{
  /* Modifies array f such that 
       f[j] = the coefficient of x^j in the expansion of (x+a)^l (x+b)^m
     Works for l, m = 0 .. 4 (s,p,d,f,g functions).
  */

  int n;
  double c;

  /* Make sure l < m, else swap. */
  if (m < l) 
  { 
    n = m; m = l; l = n; 
    c = a; a = b; b = c;
  }

  /* Make sure we're within the valid range of the function */
  if (l < 0 || l > 4) {printf("Invalid l value: %i\n", l); exit(1);}
  if (m < 0 || m > 4) {printf("Invalid m value: %i\n", m); exit(1);}
 
  if (l == 0)
  {
    if (m == 0)
    {
      f[0] = 1;
    }
    else if (m == 1)
    {
      f[0] = b;
      f[1] = 1;
    }
    else if (m == 2)
    {
      f[0] = b * b;
      f[1] = 2 * b;
      f[2] = 1;
    }
    else if (m == 3)
    {
      f[0] = b * b * b;
      f[1] = 3 * b * b;
      f[2] = 3 * b;
      f[3] = 1;
    }
    else if (m == 4)
    {
      f[0] = b * b * b * b;
      f[1] = 4 * b * b * b;
      f[2] = 6 * b * b;
      f[3] = 4 * b;
      f[4] = 1;
    }
  }
  else if (l == 1)
  {
    if (m == 1)
    {
      f[0] = a * b;
      f[1] = a + b;
      f[2] = 1;
    }
    else if (m == 2)
    {
      f[0] = a * b * b;
      f[1] = (2 * a + b) * b;
      f[2] = a + 2 * b;
      f[3] = 1;
    }
    else if (m == 3)
    {
      f[0] = a * b * b * b;
      f[1] = b * b * (3 * a + b); 
      f[2] = 3 * b * (a + b);
      f[3] = a + 3 * b;
      f[4] = 1;
    }
    else if (m == 4)
    {
      f[0] = a * b * b * b * b;
      f[1] = b * b * b * (4 * a + b);
      f[2] = 2 * b * b * (3 * a + 2 * b);
      f[3] = 2 * b * (2 * a + 3 * b);
      f[4] = a + 4 * b;
      f[5] = 1;
    }
  }
  else if (l == 2)
  {
    if (m == 2)
    {
      f[0] = a * a * b * b;
      f[1] = 2 * a * b * (a + b);
      f[2] = a * a + 4 * a * b + b * b;
      f[3] = 2 * (a + b);
      f[4] = 1;
    }
    else if (m == 3)
    {
      f[0] = a * a * b * b * b;
      f[1] = a * b * b * (3 * a + 2 * b);
      f[2] = b * (3 * a * a + 6 * a * b + b * b);
      f[3] = a * a + 6 * a * b + 3 * b * b;
      f[4] = 2 * a + 3 * b;
      f[5] = 1;
    }
    else if (m == 4)
    {
      f[0] = a * a * b * b * b * b;
      f[1] = 2 * a * b * b * b * (2 * a + b);
      f[2] = b * b * (6 * a * a + 8 * a * b + b * b);
      f[3] = 4 * b * (a * a + 3 * a * b + b * b);
      f[4] = a * a + 8 * a * b + 6 * b * b;
      f[5] = 2 * (a + 2 * b);
      f[6] = 1;
    }
  }
  else if (l == 3)
  {
    if (m == 3)
    {
      f[0] = a * a * a * b * b * b;
      f[1] = 3 * a * a * b * b * (a + b);
      f[2] = 3 * a * b * (a * a + 3 * a * b + b * b);
      f[3] = a * a * a + 9 * a * a * b + 9 * a * b * b + b * b * b;
      f[4] = 3 * (a * a + 3 * a * b + b * b);
      f[5] = 3 * (a + b);
      f[6] = 1;
    }
    else if (m == 4)
    {
      f[0] = a * a * a * b * b * b * b;
      f[1] = a * a * b * b * b * (4 * a + 3 * b);
      f[2] = 3 * a * b * b * (2 * a * a + 4 * a * b + b * b);
      f[3] = b * (4 * a * a * a + 18 * a * a * b + 12 * a * b * b + b * b * b);
      f[4] = a * a * a + 12 * a * a * b + 18 * a * b * b + 4 * b * b * b;
      f[5] = 3 * (a * a + 4 * a * b + 2 * b * b);
      f[6] = 3 * a + 4 * b;
      f[7] = 1;
    }
  }
  else if (l == 4)
  {
    if (m == 4)
    {
      f[0] = a * a * a * a * b * b * b * b;
      f[1] = 4 * a * a * a * b * b * b * (a + b);
      f[2] = 2 * a * a * b * b * (3 * a * a + 8 * a * b + 3 * b * b);
      f[3] = 4 * a * b * (a * a * a + 6 * a * a * b + 6 * a * b * b + b * b * b);
      f[4] = a * a * a * a + 16 * a * a * a * b + 36 * a * a * b * b + 16 * a * b * b * b + b * b * b * b;
      f[5] = 4 * (a * a * a + 6 * a * a * b + 6 * a * b * b + b * b * b);
      f[6] = 6 * a * a + 16 * a * b + 6 * b * b;
      f[7] = 4 * (a + b);
      f[8] = 1;
    }
  }
}

void QMCPsiPotential::g_function(int suml, double *f, double ig, double p, double *g)
{
  /* G function as defined in "Gaussian-Expansion Methods for Molecular Integrals" 
     J. Phys. Soc. Japan 21(11) 1966, p. 2319, eqn. 3.2.  Helps evaluate nuclear
     attraction integral. 

     suml = l1 + l2 (can range from 0 to 8, so that l1, l2 = s, p, d, f, g.
     f = binomial_f array as defined above.
     ig = 1 / (exponent1 + exponent2), inverse gamma.
     p = distance between nuclear center and product gaussian center.
     g = output of g_function, values from 0 to suml.
  */

  /* Make sure we're within the valid range of the function */
  if (suml > 8) {printf("Invalid suml value: %i\n", suml); exit(1);}

  if (suml == 0)
  {
    g[0] = 1;
  }
  else if (suml == 1)
  {
    g[0] = f[0];
    g[1] = -p;
  }
  else if (suml == 2)
  {
    g[0] = f[0] + 0.5 * ig;
    g[1] = -p * f[1] - 0.5 * ig;
    g[2] = p * p;
  }
  else if (suml == 3)
  {
    g[0] = f[0] + 0.5 * ig * f[2];
    g[1] = -p * f[1] - 0.5 * ig * f[2] - 1.5 * ig * p;
    g[2] = p * p * f[2] + 1.5 * ig * p;
    g[3] = -p * p * p;
  }
  else if (suml == 4)
  {
    g[0] = f[0] + 0.5 * ig * f[2] + 0.75 * ig * ig;
    g[1] = -p * f[1] - 0.5 * ig * f[2] - 1.5 * ig * p * f[3] - 1.5 * ig * ig;
    g[2] = p * p * f[2] + 1.5 * ig * p * f[3] + 0.75 * ig * ig + 3 * ig * p * p;
    g[3] = -p * p * p * f[3] - 3 * ig * p * p;
    g[4] = p * p * p * p;
  }
  else if (suml == 5)
  {
    g[0] = f[0] + 0.5 * ig * f[2] + 0.75 * ig * ig * f[4];
    g[1] = -p * f[1] - 0.5 * ig * f[2] - 1.5 * ig * p * f[3] - 1.5 * ig * ig * f[4] - 3.75 * ig * ig * p;
    g[2] = p * p * f[2] + 1.5 * ig * p * f[3] + 0.75 * ig * ig * f[4] + 3 * ig * p * p * f[4] + 7.5 * ig * ig * p;
    g[3] = -p * p * p * f[3] - 3 * ig * p * p * f[4] - 3.75 * ig * ig * p * f[5] - 5 * ig * p * p * p;
    g[4] = p * p * p * p * f[4] + 5 * ig * p * p * p;
    g[5] = -p * p * p * p * p;
  }
  else if (suml == 6)
  {
    g[0] = f[0] + 0.5 * ig * f[2] + 0.75 * ig * ig * f[4] + 1.875 * ig * ig * ig;
    g[1] = -p * f[1] - 0.5 * ig * f[2] - 1.5 * ig * p * f[3] - 1.5 * ig * ig * f[4] - 3.75 * ig * ig * p * f[5] - 5.625 * ig * ig * ig;
    g[2] = p * p * f[2] + 1.5 * ig * p * f[3] + 0.75 * ig * ig * f[4] + 3 * ig * p * p * f[4] + 7.5 * ig * ig * p * f[5] + 5.625 * ig * ig * ig + 11.25 * ig * ig * p * p;
    g[3] = -p * p * p * f[3] - 3 * ig * p * p * f[4] - 3.75 * ig * ig * p * f[5] - 5 * ig * p * p * p * f[5] - 1.875 * ig * ig * ig * f[6] - 22.5 * ig * ig * p * p;
    g[4] = p * p * p * p * f[4] + 5 * ig * p * p * p * f[5] + 11.25 * ig * ig * p * p + 7.5 * ig * p * p * p * p;
    g[5] = -p * p * p * p * p * f[5] - 7.5 * ig * p * p * p * p;
    g[6] = p * p * p * p * p * p;
  }
  else if (suml == 7)
  {
    g[0] = f[0] + 0.5 * ig * f[2] + 0.75 * ig * ig * f[4] + 1.875 * ig * ig * ig * f[6];
    g[1] = -p * f[1] - 0.5 * ig * f[2] - 1.5 * ig * p * f[3] - 1.5 * ig * ig * f[4] - 3.75 * ig * ig * p * f[5] - 5.625 * ig * ig * ig * f[6] - 13.125 * ig * ig * ig * p;
    g[2] = p * p * f[2] + 1.5 * ig * p * f[3] + 0.75 * ig * ig * f[4] + 3 * ig * p * p * f[4] + 7.5 * ig * ig * p * f[5] + 5.625 * ig * ig * ig * f[6] + 11.25 * ig * ig * p * p * f[6] + 39.375 * ig * ig * ig * p;
    g[3] = -p * p * p * f[3] - 3 * ig * p * p * f[4] - 3.75 * ig * ig * p * f[5] - 5 * ig * p * p * p * f[5] - 1.875 * ig * ig * ig * f[6] - 22.5 * ig * ig * p * p * f[6] - 39.375 * ig * ig * ig * p - 26.25 * ig * ig * p * p * p;
    g[4] = p * p * p * p * f[4] + 5 * ig * p * p * p * f[5] + 11.25 * ig * ig * p * p * f[6] + 7.5 * ig * p * p * p * p * f[6] + 13.125 * ig * ig * ig * p + 52.5 * ig * ig * p * p * p;
    g[5] = -p * p * p * p * p * f[5] - 7.5 * ig * p * p * p * p * f[6] - 26.25 * ig * ig * p * p * p - 10.5 * ig * p * p * p * p * p;
    g[6] = p * p * p * p * p * p * f[6] + 10.5 * ig * p * p * p * p * p;
    g[7] = -p * p * p * p * p * p * p; 
  }
  else if (suml == 8)
  {
    g[0] = f[0] + 0.5 * ig * f[2] + 0.75 * ig * ig * f[4] + 1.875 * ig * ig * ig * f[6] + 6.5625 * ig * ig * ig * ig * f[8];
    g[1] = -p * f[1] - 0.5 * ig * f[2] - 1.5 * ig * p * f[3] - 1.5 * ig * ig * f[4] - 3.75 * ig * ig * p * f[5] - 5.625 * ig * ig * ig * f[6] - 13.125 * ig * ig * ig * p * f[7] - 26.25 * ig * ig * ig * ig;
    g[2] = p * p * f[2] + 1.5 * ig * p * f[3] + 0.75 * ig * ig * f[4] + 3 * ig * p * p * f[4] + 7.5 * ig * ig * p * f[5] + 5.625 * ig * ig * ig * f[6] + 11.25 * ig * ig * p * p * f[6] + 39.375 * ig * ig * ig * p * f[7] + 39.375 * ig * ig * ig * ig + 52.5 * ig * ig * ig * p * p;
    g[3] = -p * p * p * f[3] - 3 * ig * p * p * f[4] - 3.75 * ig * ig * p * f[5] - 5 * ig * p * p * p * f[5] - 1.875 * ig * ig * ig * f[6] - 22.5 * ig * ig * p * p * f[6] - 39.375 * ig * ig * ig * p * f[7] - 26.25 * ig * ig * p * p * p * f[7] - 26.25 * ig * ig * ig * ig - 157.5 * ig * ig * ig * p * p;
    g[4] = p * p * p * p * f[4] + 5 * ig * p * p * p * f[5] + 11.25 * ig * ig * p * p * f[6] + 7.5 * ig * p * p * p * p * f[6] + 13.125 * ig * ig * ig * p * f[7] + 52.5 * ig * ig * p * p * p * f[7] + 6.5625 * ig * ig * ig * ig * f[8] + 157.5 * ig * ig * ig * p * p + 52.5 * ig * ig * p * p * p * p;
    g[5] = -p * p * p * p * p * f[5] - 7.5 * ig * p * p * p * p * f[6] - 26.25 * ig * ig * p * p * p * f[7] - 10.5 * ig * p * p * p * p * p * f[7] - 52.5 * ig * ig * ig * p * p - 105 * ig * ig * p * p * p * p;
    g[6] = p * p * p * p * p * p * f[6] + 10.5 * ig * p * p * p * p * p * f[7] + 52.5 * ig * ig * p * p * p * p + 14 * ig * p * p * p * p * p * p;
    g[7] = -p * p * p * p * p * p * p * f[7] - 14 * ig * p * p * p * p * p * p;
    g[8] = p * p * p * p * p * p * p * p;
  } 
}
