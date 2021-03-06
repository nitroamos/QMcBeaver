//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2000-2.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

#include "MathFunctions.h"
#include "Random.h"

double MathFunctions::erf(double x)
{
  // Series on [0,1]
  int ERFC_COEF_LENGTH = 14;
  double ERFC_COEF[] = 
  {
    -.490461212346918080399845440334e-1,
    -.142261205103713642378247418996e0,
    .100355821875997955757546767129e-1,
    -.576876469976748476508270255092e-3,
    .274199312521960610344221607915e-4,
    -.110431755073445076041353812959e-5,
    .384887554203450369499613114982e-7,
    -.118085825338754669696317518016e-8,
    .323342158260509096464029309534e-10,
    -.799101594700454875816073747086e-12,
    .179907251139614556119672454866e-13,
    -.371863548781869263823168282095e-15,
    .710359900371425297116899083947e-17,
    -.126124551191552258324954248533e-18
  };

  double  ans;
  double  y = fabs(x);
  
  if (y <= 1.49012e-08) 
    {
      // 1.49012e-08 = Math.sqrt(2*EPSILON_SMALL)
      ans = 2 * x / 1.77245385090551602729816748334;
    } 
  else if (y <= 1) 
    {
      ans = x * (1 + csevl(2 * x * x - 1, ERFC_COEF, ERFC_COEF_LENGTH));
    } 
  else if (y < 6.013687357) 
    {
      // 6.013687357 = 
      // Math.sqrt(-Math.log(1.77245385090551602729816748334 * EPSILON_SMALL))
      ans = sign(1 - erfc(y), x);
    } 
  else 
    {
      ans = sign(1, x);
    }

  return ans;
}

double MathFunctions::erfc(double x)
{
  // Series on [0,1]
  int ERFC_COEF_LENGTH = 14;
  static double ERFC_COEF[] = 
  {
    -.490461212346918080399845440334e-1,
    -.142261205103713642378247418996e0,
    .100355821875997955757546767129e-1,
    -.576876469976748476508270255092e-3,
    .274199312521960610344221607915e-4,
    -.110431755073445076041353812959e-5,
    .384887554203450369499613114982e-7,
    -.118085825338754669696317518016e-8,
    .323342158260509096464029309534e-10,
    -.799101594700454875816073747086e-12,
    .179907251139614556119672454866e-13,
    -.371863548781869263823168282095e-15,
    .710359900371425297116899083947e-17,
    -.126124551191552258324954248533e-18
  };

  // Series on [0.25,1.00]
  int ERFC2_COEF_LENGTH = 27;
  static double ERFC2_COEF[] = 
  {
    -.69601346602309501127391508262e-1,
    -.411013393626208934898221208467e-1,
    .391449586668962688156114370524e-2,
    -.490639565054897916128093545077e-3,
    .715747900137703638076089414183e-4,
    -.115307163413123283380823284791e-4,
    .199467059020199763505231486771e-5,
    -.364266647159922287393611843071e-6,
    .694437261000501258993127721463e-7,
    -.137122090210436601953460514121e-7,
    .278838966100713713196386034809e-8,
    -.581416472433116155186479105032e-9,
    .123892049175275318118016881795e-9,
    -.269063914530674343239042493789e-10,
    .594261435084791098244470968384e-11,
    -.133238673575811957928775442057e-11,
    .30280468061771320171736972433e-12,
    -.696664881494103258879586758895e-13,
    .162085454105392296981289322763e-13,
    -.380993446525049199987691305773e-14,
    .904048781597883114936897101298e-15,
    -.2164006195089607347809812047e-15,
    .522210223399585498460798024417e-16,
    -.126972960236455533637241552778e-16,
    .310914550427619758383622741295e-17,
    -.766376292032038552400956671481e-18,
    .190081925136274520253692973329e-18
  };

  // Series on [0,0.25]
  int ERFCC_COEF_LENGTH = 29;
  static double ERFCC_COEF[] = 
  {
    .715179310202924774503697709496e-1,
    -.265324343376067157558893386681e-1,
    .171115397792085588332699194606e-2,
    -.163751663458517884163746404749e-3,
    .198712935005520364995974806758e-4,
    -.284371241276655508750175183152e-5,
    .460616130896313036969379968464e-6,
    -.822775302587920842057766536366e-7,
    .159214187277090112989358340826e-7,
    -.329507136225284321486631665072e-8,
    .72234397604005554658126115389e-9,
    -.166485581339872959344695966886e-9,
    .401039258823766482077671768814e-10,
    -.100481621442573113272170176283e-10,
    .260827591330033380859341009439e-11,
    -.699111056040402486557697812476e-12,
    .192949233326170708624205749803e-12,
    -.547013118875433106490125085271e-13,
    .158966330976269744839084032762e-13,
    -.47268939801975548392036958429e-14,
    .14358733767849847867287399784e-14,
    -.444951056181735839417250062829e-15,
    .140481088476823343737305537466e-15,
    -.451381838776421089625963281623e-16,
    .147452154104513307787018713262e-16,
    -.489262140694577615436841552532e-17,
    .164761214141064673895301522827e-17,
    -.562681717632940809299928521323e-18,
    .194744338223207851429197867821e-18
  };

  double ans;
  double y = fabs(x);
  
  if (x <= -6.013687357) 
    {
      // -6.013687357 = 
      // -Math.sqrt(-Math.log(1.77245385090551602729816748334 * EPSILON_SMALL))
      ans = 2;
    } 
  else if (y < 1.49012e-08) 
    {
      // 1.49012e-08 = Math.sqrt(2*EPSILON_SMALL)
      ans = 1 - 2*x/1.77245385090551602729816748334;
    } 
  else 
    {
      double ysq = y*y;
    
      if (y < 1) 
	{
	  ans = 1 - x*(1+csevl(2*ysq-1,ERFC_COEF,ERFC_COEF_LENGTH));
	} 
      else if (y <= 4.0) 
	{
	  ans = exp(-ysq)/y*(0.5+csevl((8.0/ysq-5.0)/3.0,ERFC2_COEF,
				       ERFC2_COEF_LENGTH));
	  if (x < 0)  ans = 2.0 - ans;
	  if (x < 0)  ans = 2.0 - ans;
	  if (x < 0)  ans = 2.0 - ans;
	} 
      else 
	{
	  ans = exp(-ysq)/y*(0.5+csevl(8.0/ysq-1,ERFCC_COEF,
				       ERFCC_COEF_LENGTH));
	  if (x < 0)  ans = 2.0 - ans;
	}
    }

  return ans;
}

void MathFunctions::fitU(double Z, double Frk, double r, double delta,
			 double & zeta, double & a, double & I)
{
  //Now we solve for a and z (looking at eqns 41 & 42)
  double t1 = Frk + Z;
  double discriminant = t1*(r*t1 - 4.0)/r;
  zeta = 0.0;
  if(discriminant > 0.0)
    {
      double root1 = Z - Frk;
      double root2 = Z - Frk;
      double sq = sqrt(discriminant);	  
      root1 -= sq;
      root2 += sq;
      root1 /= 2.0;
      root2 /= 2.0;
      
      if(root1*root2 <= 0.0)
	{
	  //return the postive root
	  zeta = root1 > root2 ? root1 : root2;
	} else if(root1 <= 0.0)
	  {
	    //both were negative
	    zeta = 0.0;	     
	  } else {
	    //both positive, return the smallest
	    zeta = root1 < root2 ? root1 : root2;
	  }
    }
  
  if(zeta <= 1e-50)
    {
      if(Frk < 0.0)
	zeta = -Frk;
      else
	zeta = 1.0;
    }
  
  a = -(Frk + zeta)/(-1.0 + r*(Frk + zeta));  

  /*
    This is from eqns on pg 158 (In the appendix).
    Note: there is a typo in the formulas. Umrigar says the nodes
    are at -a, but as is easily verified, the nodes are actually at -1/a.
    
    Because we're integrating the absolute value of the function, we need
    to separate this into two integrals if the integrand goes through a node.
  */
  double node = -1.0/a;
  if(fabs(a) > 1e-50 && r/delta < node && node < r*delta)
    {
      I = fabs(accel_G(1.5,zeta*r*delta,-zeta/a,a,zeta)) +
          fabs(accel_G(1.5,-zeta/a,zeta*r/delta,a,zeta));
    } else {
      I = fabs(accel_G(1.5,zeta*r*delta,zeta*r/delta,a,zeta));
    }

  /*
  printf("a=%20.10f zeta=%20.10f r=%20.10f root=%20.10f L=%20.10f U=%20.10f I=%20.10f\n",
	 a,zeta,r,node,r/delta,r*delta,I);
  */
}

double MathFunctions::accel_G(double n, double x, double y,
			      double a, double zeta)
{
  double den = pow(zeta,-n-1);
  double t1 = (gamma_inc(n,x)   - gamma_inc(n,y)  )*den*zeta;
  double t2 = (gamma_inc(n+1,x) - gamma_inc(n+1,y))*den;
  return t1 + a*t2;
}

double MathFunctions::csevl(double x, double * coef, int lengthCoef)
{
  double b0, b1, b2, twox;
  int	 i;

  b1 = 0.0;
  b0 = 0.0;
  b2 = 0.0;
  twox = 2.0*x;

  for (i = lengthCoef-1;  i >= 0;  i--) 
    {
      b2 = b1;
      b1 = b0;
      b0 = twox*b1 - b2 + coef[i];
    }

  return 0.5*(b0-b2);
}

double MathFunctions::sign(double x, double y)
{
  double abs_x = ((x < 0) ? -x : x);
  return (y < 0.0) ? -abs_x : abs_x;
}

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

double MathFunctions::F_gamma(double v, double t)
{
  /* Computes the function 
      Integrate[u^(2 v) Exp(-t u^2), {u, 0, 1}]
      = (1/2) t^(-1/2 - v) * Gamma(1/2 + v, t)
  */
  if (t < 0.000001) 
    return 1.0 / (2.0 * v + 1.0); /* limiting case */
  else
    return 0.5 * pow(t, -0.5 - v) * gamma_inc(0.5 + v, t);
}

double MathFunctions::gamma_inc(double a, double x)
{
  /* From Numerical Recipes in Fortran, gammp.f 
     Returns the incomplete gamma function 
        Integrate[Exp[-t] t^(a - 1), {t, 0, x}] / Gamma[a].
     Picks the series or continued fraction algorithm depending on input. 
  */

  if (x < 0 || a < 0) {printf("Bad argument in gamma function.\n"); exit(1);}
  if (x < a + 1)
    return gamma_series(a, x);
  else
    return gamma_cf(a, x);
}

double MathFunctions::gamma_log(double x2)
{
  /* From Numerical Recipes in C, gammln.c
     Returns the value Log(gamma(x)), x > 0
  */
  
  double x, y, tmp, ser;
  static double cof[] = {76.18009172947146, -86.50532032941677, 
			 24.01409824083091, -1.231739572450155,
			 0.1208650973866179E-2, -0.5395239384953E-5, 
			 2.5066282746310005};
  int j;

  y = x = x2;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.000000000190015;

  for (j = 0; j < 6; j++)
    {
      y += 1.0;
      ser +=  cof[j] / y;
    }
  return -tmp + log(2.5066282746310005 * ser/x);
}

double MathFunctions::gamma_series(double a, double x)
{
  /* From Numerical Recipes in C, gser.c
     Returns the incomplete gamma function
        Integrate[Exp[-t] t^(a - 1), {t, 0, x}].
     Uses a series expansion.
  */

  if (x <= 0.0) {printf("x less than 0 in routine gamma_series.\n"); exit(1);}
 
  int n; 
  double sum, del, ap;

  ap = a;
  del = sum = 1.0 / a;
  for (n = 1; n <= ITMAX; n++)
  {
    ++ap;
    del *= x / ap;
    sum += del;
    if (fabs(del) < fabs(sum) * EPS)
      return sum * exp(-x + a * log(x));
  }
  printf("a too large, ITMAX too small in routine gamma_series.\n"); exit(1);
  return 0.0;
}

double MathFunctions::gamma_cf(double a, double x)
{
  /* From Numerical Recipes in C, gcf.c
     Returns the incomplete gamma function
        Integrate[Exp[-t] t^(a - 1), {t, 0, x}] / Gamma[a].
     Uses a continued fraction expansion.
  */

  double an, b, c, d, del, h;

  b = x + 1.0 - a;
  c = 1.0 / FPMIN;
  d = 1.0 / b;
  h = d;
  int i;
  for (i = 1; i <= ITMAX; i++)
  {
    an = -i * (i - a);
    b += 2.0;
    d = an * d + b;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = b + an / c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    del = d * c;
    h *= del;
    if (fabs(del - 1.0) < EPS) break;
  }
  if (i > ITMAX) 
    {
      printf("a too large, ITMAX too small in gamma_cf\n"); 
      exit(1);
    }
  return exp(gamma_log(a)) - exp(-x + a * log(x)) * h;
}

double MathFunctions::legendre(int l, double x)
{
  if ( l == 0 )
    return 1.0;
  if ( l == 1 )
    return x;
  
  return ( (2.0*l-1.0)*x*legendre(l-1, x) - (l-1.0)*legendre(l-2, x) ) / l;
}

void MathFunctions::randomlyRotate(Array2D<double> & points)
{
  double phi = ran.unidev()*2*pi;
  double theta = ran.sindev();
  double angle = ran.unidev()*2*pi;
  
  Array1D<double> axis(3);
  axis(0) = sin(theta)*cos(phi);
  axis(1) = sin(theta)*sin(phi);
  axis(2) = cos(theta);
  
  for (int i=0; i<points.dim1(); i++)
    {
      Array1D<double> temp_coords(3);
      for (int j=0; j<3; j++) temp_coords(j) = points(i,j);
      temp_coords.rotate(axis,angle);
      for (int k=0; k<3; k++) points(i,k) = temp_coords(k);
    }
}

double MathFunctions::rij(Array2D<double> &position, int i, int j)
{
  double r;
  r = sqrt((position(i,0)-position(j,0))*
           (position(i,0)-position(j,0))+
           (position(i,1)-position(j,1))*
           (position(i,1)-position(j,1))+
           (position(i,2)-position(j,2))*
           (position(i,2)-position(j,2)));
  return r;
}

double MathFunctions::rij(Array2D<double> &positioni,
			  Array2D<double> &positionj,
			  int i, int j)
{
  double r;
  r = sqrt((positioni(i,0)-positionj(j,0))*
           (positioni(i,0)-positionj(j,0))+
           (positioni(i,1)-positionj(j,1))*
           (positioni(i,1)-positionj(j,1))+
           (positioni(i,2)-positionj(j,2))*
           (positioni(i,2)-positionj(j,2)));
  return r;
}

#define EXP_A (1048576/M_LN2) 
#define EXP_C 60801 	
double MathFunctions::fast_exp(double y)
{
  if(y < -600.0) return 0.0;
  union 
  { 
    double d; 
#ifdef LITTLE_ENDIAN 
    struct { int j, i; } n; 
#else 
    struct { int i, j; } n; 
#endif 
  } 
  _eco; 
  _eco.n.i = (int)(EXP_A*(y)) + (1072693248 - EXP_C); 
  _eco.n.j = 0; 
  return  _eco.d; 
} 
