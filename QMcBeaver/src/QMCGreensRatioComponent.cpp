#include "QMCGreensRatioComponent.h"

QMCGreensRatioComponent::QMCGreensRatioComponent()
{
  initialize();
}
QMCGreensRatioComponent::QMCGreensRatioComponent(double value)
{
  initialize();
  k = value;
}
QMCGreensRatioComponent::QMCGreensRatioComponent(double w, double x, double y,\
						 double z)
{
  initialize();
  k = w;
  a = x;
  b = y;
  c = z;
}
QMCGreensRatioComponent::QMCGreensRatioComponent( const \
						  QMCGreensRatioComponent & rhs )
{
  *this = rhs;
}
QMCGreensRatioComponent::~QMCGreensRatioComponent()
{
  // does nothing
}
void QMCGreensRatioComponent::operator=( const QMCGreensRatioComponent & rhs )
{
  k = rhs.k;
  a = rhs.a;
  b = rhs.b;
  c = rhs.c;
}
void QMCGreensRatioComponent::initialize()
{
  k = 1.0;
  a = 1.0;
  b = 1.0;
  c = 0.0;
}
double QMCGreensRatioComponent::getValue()
{
  return k*pow(a,b)*exp(c);
}
double QMCGreensRatioComponent::DivideBy(QMCGreensRatioComponent & denom)
{
  // Get resulting exp(stuff) -> get "stuff left over
  double EXP_EXP = c-denom.c;
  // Handle pow(a,b)/pow(d.a,d.b)
  double POW_A, POW_B;
  SimplifyRatioPowers(a,b,denom.a,denom.b,POW_A,POW_B);
  // Get resulting BASE_LEFT(EXP_LEFT)
  double BASE_LEFT, EXP_LEFT;
  SimplifyRatioPowers(POW_A,POW_B,exp(1.0),-1.0*EXP_EXP,BASE_LEFT,EXP_LEFT);
      
  /* If the code has successfully handled largeness/smallness in other parts
     of this class, then because forward/reverse green's functions have similar
     magnitudes, there should be no need to protect any function here.*/
  return k/denom.k * pow(BASE_LEFT,EXP_LEFT);
}

void QMCGreensRatioComponent::MultiplyBy(QMCGreensRatioComponent &X)
{
  double temp_k, temp_a, temp_b, temp_c;

  /*If k and X.k are both very large or very small, their product could overflow or underflow
    the temp_k variable. This is actually likely because several MultiplyBy iterations are called
    for molecules with several molecules. If the configuration isn't equilibrated, then very
    large/small values will start accumulating in the k degree of freedom. So, if this happens,
    we need to shift the large/small exponents over to the c degree of freedom.*/
  temp_k = k * X.k;
  if(fabs(temp_k) > 1e300 || fabs(temp_k) < 1e-300){
    c += log(k);
    X.c += log(X.k);
    k = 1.0;
    X.k = 1.0;
    temp_k = 1.0;
  }

  temp_c = c + X.c;

  if ( fabs(a-1.0) < 1e-15 )
    {
      if ( fabs(X.a-1.0) < 1e-15 )
	{
	  temp_a = 1.0;
	  temp_b = 0.0;
	}
      else
	{
	  temp_a = X.a;
	  temp_b = X.b;
	}
    }
  else if ( fabs(X.a-1.0) < 1e-15 )
    {
      if ( fabs(a-1.0) < 1e-15 )
	{
	  temp_a = 1.0;
	  temp_b = 0.0;
	}
      else
	{
	  temp_a = a;
	  temp_b = b;
	}
    }
  else if ( fabs(a-X.a) < 1e-15 )
    {
      temp_a = a;
      temp_b = b + X.b;
    }
  else if (a < X.a)
    {
      temp_a = X.a;
      temp_b = X.b + b*log(a)/log(X.a);
    }
  else
    {
      temp_a = a;
      temp_b = b + X.b*log(X.a)/log(a);
    }
  k = temp_k;
  a = temp_a;
  b = temp_b;
  c = temp_c;
}
QMCGreensRatioComponent QMCGreensRatioComponent::operator + \
(const QMCGreensRatioComponent &X)
{
  if(k == 0){
    return X;
  } else if(X.k == 0){
    return *this;
  }
  
  QMCGreensRatioComponent result;
  double POW_A, POW_B, expArg = X.c - c;
  double temp, value, dumb;
  SimplifyRatioPowers(X.a,X.b,a,b,POW_A,POW_B);
 
  //exp(-690) ~= 1e-300
  /* This if statement attempts to handle extreme cases by avoiding the explicit
     evaluation of exp on very large/small numbers. It definitely represents the
     "long way around", but I couldn't get the QSC at LANL to run without this
     apparent (apparently not) extremism. I wanted to avoid correcting the
     overflow/underflow possibilities by setting exp(large) = 690 because the
     exponent is possibly important to the rest of the number. Anyway, as we try
     larger and larger molecules, this code needs to be very robust.*/
  if(fabs(expArg) > 690){

    /* This next little section is very (i think) dumb. This is to correct a code
       crash (on the QSC) when somehow, a value smaller than 1e-308 was assigned to
       X.k (but somehow not caught by the X.k == 0 check above) which would crash
       the evaluation of log(X.k). I don't understand how this is possible... but
       this fixes it.*/
    if(fabs(X.k) < 1e-300){
      dumb = 1e100 * X.k;
      expArg += log(dumb) + log(1e-100);
    } else if(fabs(X.k) > 1e300){
      dumb = 1e-100 * X.k;
      expArg += log(dumb) + log(1e100);
    } else {
      expArg += log(X.k);
    }

    expArg += log(pow(POW_A,POW_B));
    temp = log(k);
    /* I want to allow for any possiblity in the comparison between temp
       and expArg*/
    if(fabs(temp - expArg) < 15){
      result.k = 1.0 + exp(temp - expArg);
      result.c = c + expArg;
    } else if(expArg > temp){
      result.k = 1.0;
      result.c = c + expArg;
    } else {
      result.k = 1.0;
      result.c = c + temp;
    }
  } else {
    result.k = k + X.k * pow(POW_A,POW_B) * exp(expArg);
    result.c = c;
  }

  result.a = a;
  result.b = b;
  return result;
}
  
void QMCGreensRatioComponent::toXML(ostream & strm)
{
  strm << "<QMCGreensRatioComponent>" << endl;
  strm << "\t<k>\t" << k << "\t</k>" << endl;
  strm << "\t<a>\t" << a << "\t</a>" << endl;
  strm << "\t<b>\t" << b << "\t</b>" << endl;
  strm << "\t<c>\t" << c << "\t</c>" << endl;
  strm << "</QMCGreensRatioComponent>\n" << endl;
}
// ra^rb = na^nb/da^db
void QMCGreensRatioComponent::SimplifyRatioPowers(double na, double nb, \
						  double da, double db, double &ra, double &rb)
{
  double power1, power2;
  int method = -1;
#ifdef CRAY_X1_DEBUG
  cout << "----------------------------------------------------------" << endl;
  cout << "na:\t" << na << "\tnb:\t" << nb << "\tda:\t" << da << "\tdb:\t";
  cout << db << endl;
  cout << "&&&\t" << fabs(na-da) << endl;
#endif
  if ( fabs(na-da) < 1e-15 )
    {
      // They already have the same base
      method = 0;
    }
  
  else if ( na>0.0 && da>0.0 )
    {
      if ( fabs(log(na)) > fabs(log(da)) )
	{
	  method = 1;
	}
      else
	{
	  method = 2;
	}
    }
#ifdef CRAY_X1_DEBUG
  cout << "method:\t" << method << endl;
#endif
  if (method == 0)
    {
      ra = na;
      rb = nb-db;
    }
  else if (method == 1)
    {
      ra = na;
      power1 = nb;
      power2 = log(da)/log(na)*db;
      rb = power1 - power2;
    }
  
  else if (method == 2)
    {
      ra = da;
      power1 = log(na)/log(da)*nb;
      power2 = db;
      rb = power1 - power2;
    }
  else
    {
      cerr << "Error in QMCGreensRatioComponent::SimplifyRatioPowers()";
      cerr << endl;
      cerr << "Your query is not yet supported." << endl;
      cerr << "na:\t" << na << endl;
      cerr << "nb:\t" << nb << endl;
      cerr << "da:\t" << da << endl;
      cerr << "db:\t" << db << endl;
      exit(1);
    }
#ifdef CRAY_X1_DEBUG
  cout << "RESULTS:\tra:\t" << ra << "\trb:\t" << rb << endl;
  cout << "----------------------------------------------------------" << endl;
#endif
}