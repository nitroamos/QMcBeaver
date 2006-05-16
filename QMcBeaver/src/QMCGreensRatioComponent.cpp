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

void QMCGreensRatioComponent::initialize()
{
  k = 1.0;
  a = 1.0;
  b = 0.0;
  c = 0.0;
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
void QMCGreensRatioComponent::SimplifyRatioPowers(double na, double nb,
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
    // They already have the same base
    method = 0;
  
  else if ( na>0.0 && da>0.0 )
    {
      if ( fabs(log(na)) > fabs(log(da)) )
	method = 1;
      else
	method = 2;
    }
  if ( (fabs(na-1.0) < 1e-15) || (fabs(nb) < 1e-15) )
    method = 3;
  if ( (fabs(da-1.0) < 1e-15) || (fabs(db) < 1e-15) )
    method = 4;

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
      power2 = (log(da)/log(na))*db;
      rb = power1 - power2;
    }
  
  else if (method == 2)
    {
      ra = da;
      power1 = (log(na)/log(da))*nb;
      power2 = db;
      rb = power1 - power2;
    }
  else if (method == 3)
    {
      ra = da;
      rb = -db;
    }
  else if (method == 4)
    {
      ra = na;
      rb = nb;
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

double QMCGreensRatioComponent::getValue() const
{
  if (a<0.0)
    {
      cerr << "Error in QMCGreensRatioComponent::getValue():" << endl;
      cerr << "attempting to take a power of a negative number.";
      cerr << endl;
      cerr << "k:\t" << k << endl;
      cerr << "a:\t" << a << endl;
      cerr << "b:\t" << b << endl;
      cerr << "c:\t" << c << endl;
      exit(1);
    }      
  if (IeeeMath::isNaN(k) || IeeeMath::isNaN(a) || IeeeMath::isNaN(b) || IeeeMath::isNaN(c) )
    {
      cerr << "Error in QMCGreensRatioComponent::getValue():" << endl;
      cerr << "k:\t" << k << endl;
      cerr << "a:\t" << a << endl;
      cerr << "b:\t" << b << endl;
      cerr << "c:\t" << c << endl;
      exit(1);
    }      
  return k*pow(a,b)*exp(c);
}

QMCGreensRatioComponent & QMCGreensRatioComponent::divideBy(const QMCGreensRatioComponent & denom)
{
  // Handle pow(a,b)/pow(d.a,d.b)
  double POW_A, POW_B;
  SimplifyRatioPowers(a,b,denom.a,denom.b,POW_A,POW_B);
  
  k /= denom.k;
  a = POW_A;
  b = POW_B;
  c -= denom.c;

  // Here we even out the values if necessary.

  if ( (c < -10.0) && (log(fabs(k)) > 10.0) )
    {
      c += log(fabs(k));
      k = k < 0 ? -1.0 : 1.0;
    }
  else if ( (c > 10.0) && (log(fabs(k)) < -10.0) )
    {
      c += log(fabs(k));
      k = k < 0 ? -1.0 : 1.0;
    }

  if ( (c < -10.0) && (b*log(fabs(a)) > 10.0) )
    {
      c += b*log(fabs(a));
      a = 1.0;
      b = 0.0;
    }
  else if ( (c > 10.0) && (b*log(fabs(a)) < -10.0) )
    {
      c += b*log(fabs(a));
      a = 1.0;
      b = 0.0;
    }

  return *this;
}

QMCGreensRatioComponent & QMCGreensRatioComponent::multiplyBy(const QMCGreensRatioComponent &X)
{
  if (IeeeMath::isNaN(k) || IeeeMath::isNaN(a) || IeeeMath::isNaN(b) || IeeeMath::isNaN(c))
    {
      cerr << "Error in QMCGreensRatioComponent::multiplyBy():" << endl;
      cerr << "k:\t" << k << "\ta:\t" << a << "\tb:\t" << b << "\tc:\t" << c;
      cerr << endl;
      exit(1);
    }

  if (IeeeMath::isNaN(X.k) || IeeeMath::isNaN(X.a) || IeeeMath::isNaN(X.b) || IeeeMath::isNaN(X.c))
    {
      cerr << "Error in QMCGreensRatioComponent::multiplyBy():" << endl;
      cerr << "X.k:\t" << X.k << "\tX.a:\t" << X.a << "\tX.b:\t" << X.b;
      cerr << "\tX.c:\t" << X.c << endl;
      exit(1);
    }

  if (a < 0.0 || X.a < 0.0)
    {
      cerr << "Error in QMCGreensRatioComponent::multiplyBy():" << endl;
      cerr << "Attempting to take a power of a negative number." << endl;
      cerr << "a:\t" << a << "\tb:\t" << b << endl;
      cerr << "X.a\t" << X.a << "\tX.b:\t" << X.b << endl;
      exit(1);
    }

  double temp_k, temp_a, temp_b;
  
  // If k and X.k are both very large or very small, their product could 
  //  overflow or underflow the temp_k variable. This is actually likely 
  // because several MultiplyBy iterations are called for molecules with 
  // several atoms.  If the configuration isn't equilibrated, then very 
  // large/small values will start accumulating in the k degree of freedom. So,
  // if this happens, we need to shift the large/small exponents over to the c
  // degree of freedom.

  temp_k = k * X.k;
  if(fabs(temp_k) > 1e200 || fabs(temp_k) < 1e-200)
    {
      if (fabs(k) < 1e-100)
	c += log(fabs(k)*1e100) + log(1e-100);
      else if (fabs(k) > 1e100)
	c += log(fabs(k)*1e-100) + log(1e100);
      else
	c += log(fabs(k));

      if (fabs(X.k) < 1e-100)
	c += log(fabs(X.k)*1e100) + log(1e-100);
      else if (fabs(X.k) > 1e100)
	c += log(fabs(X.k)*1e-100) + log(1e100);
      else
	c += log(fabs(X.k));

      k = k < 0 ? -1.0 : 1.0;
      k = X.k < 0 ? -k : k;
    }
  else
    k = temp_k;

  c += X.c;
  
  if ( (fabs(a-1.0) < 1e-15) || (fabs(b) < 1e-15) )
    {
      if ( (fabs(X.a-1.0) < 1e-15) || (fabs(X.b) < 1e-15) )
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
  else if ( (fabs(X.a-1.0) < 1e-15) || (fabs(X.b) < 1e-15) )
    { 
      temp_a = a;
      temp_b = b;
    }	
  else if ( fabs(a-X.a) < 1e-15 )
    {
      temp_a = a;
      temp_b = b + X.b;
    } 
  else if (a < X.a)
    {
      temp_a = X.a;

      temp_b = b;
      if (fabs(a) < 1e-100)
	temp_b *= log(fabs(a)*1e100) + log(1e-100);
      else if (fabs(a) > 1e100)
	temp_b *= log(fabs(a)*1e-100) + log(1e100);
      else
	temp_b *= log(fabs(a));

      if (fabs(X.a) < 1e-100)
	temp_b /= log(fabs(X.a)*1e100) + log(1e-100);
      else if (fabs(X.a) > 1e100)
	temp_b /= log(fabs(X.a)*1e-100) + log(1e100);
      else
	temp_b /= log(fabs(X.a));

      temp_b += X.b;
    }	
  else 
    {
      temp_a = a;
      
      temp_b = X.b;
      if (fabs(X.a) < 1e-100)
	temp_b *= log(fabs(X.a)*1e100) + log(1e-100);
      else if (fabs(X.a) > 1e100)
	temp_b *= log(fabs(X.a)*1e-100) + log(1e100);
      else
	temp_b *= log(fabs(X.a));

      if (fabs(a) < 1e-100)
	temp_b /= log(fabs(a)*1e100) + log(1e-100);
      else if (fabs(a) > 1e100)
	temp_b /= log(fabs(a)*1e-100) + log(1e100);
      else
	temp_b /= log(fabs(a));

      temp_b += b;
    }
  a = temp_a;
  b = temp_b;

  // Here we even out the values if necessary.

  if ( (c < -10.0) && (log(fabs(k)) > 10.0) )
    {
      c += log(fabs(k));
      k = k < 0 ? -1.0 : 1.0;
    }
  else if ( (c > 10.0) && (log(fabs(k)) < -10.0) )
    {
      c += log(fabs(k));
      k = k < 0 ? -1.0 : 1.0;
    }

  if ( (c < -10.0) && (b*log(fabs(a)) > 10.0) )
    {
      c += b*log(fabs(a));
      a = 1.0;
      b = 0.0;
    }
  else if ( (c > 10.0) && (b*log(fabs(a)) < -10.0) )
    {
      c += b*log(fabs(a));
      a = 1.0;
      b = 0.0;
    }

  return *this;
}

QMCGreensRatioComponent & QMCGreensRatioComponent::add(const QMCGreensRatioComponent &X)
{
  // This function has been causing a lot of problems on QSC.
  // First we make sure that none of the elements is NaN.

  if (IeeeMath::isNaN(k) || IeeeMath::isNaN(a) || IeeeMath::isNaN(b) || IeeeMath::isNaN(c))
    {
      cerr << "Error in QMCGreensRatioComponent::add():" << endl;
      cerr << "k:\t" << k << "\ta:\t" << a << "\tb:\t" << b << "\tc:\t" << c;
      cerr << endl;
      exit(1);
    }

  if (IeeeMath::isNaN(X.k) || IeeeMath::isNaN(X.a) || IeeeMath::isNaN(X.b) || IeeeMath::isNaN(X.c))
    {
      cerr << "Error in QMCGreensRatioComponent::add():" << endl;
      cerr << "X.k\t" << X.k << "\tX.a:\t" << X.a << "\tX.b:\t" << X.b;
      cerr << "\tX.c:\t" << X.c << endl;      
      exit(1);
    }

  if (a < 0.0 || X.a < 0.0)
    {
      cerr << "Error in QMCGreensRatioComponent::add():" << endl;
      cerr << "Attempting to take a power of a negative number." << endl;
      cerr << "a:\t" << a << "\tb:\t" << b << endl;
      cerr << "X.a\t" << X.a << "\tX.b:\t" << X.b << endl;
    }

  if(k == 0.0)
    {
      *this = X;
      return *this;
    } 
  else if(X.k == 0.0)
    return *this;

  // If the k value is extreme, we shift it into the c degree of freedom.
  if (fabs(k) < 1e-100)
    {
      c += log(fabs(k)*1e100) + log(1e-100);
      k = k > 0.0 ? 1.0 : -1.0;
    }
  else if (fabs(k) > 1e100)
    {
      c += log(fabs(k)*1e-100) + log(1e100);
      k = k > 0.0 ? 1.0 : -1.0;
    }

  double log_middle_term;
  double expArg = X.c - c;
  double temp, Xtemp, signK = k, signXk = X.k;

  if (fabs(a) < 1e-100)
    temp = b*(log(a*1e100) + log(1e-100));

  else if (fabs(a) > 1e100)
    temp = b*(log(a*1e-100) + log(1e100));

  else
    temp = b*log(a);

  if (fabs(X.a) < 1e-100)
    Xtemp = X.b*(log(X.a*1e100) + log(1e-100));

  else if (fabs(X.a) > 1e100)
    Xtemp = X.b*(log(X.a*1e-100) + log(1e100));

  else
    Xtemp = X.b*log(X.a);

  log_middle_term = Xtemp - temp;

  // exp(-690) ~ 1e-300
  // This if statement attempts to handle extreme cases by avoiding the 
  // explicit evaluation of exp on very large/small numbers. It definitely 
  // represents the "long way around", but I couldn't get the QSC at LANL to 
  // run without this apparent (apparently not) extremism. I wanted to avoid 
  // correcting the overflow/underflow possibilities by setting 
  // exp(large) = 690 because the exponent is possibly important to the rest of
  // the number. Anyway, as we try larger and larger molecules, this code needs
  // to be very robust.
  if(fabs(expArg) > 690)
    {
      // If the value of X.k is extreme, we avoid evaluating its log directly.
      if (fabs(X.k) < 1e-100)
	expArg += log(fabs(X.k)*1e100) + log(1e-100);
      else if (fabs(X.k) > 1e100)
	expArg += log(fabs(X.k)*1e-100) + log(1e100);
      else
	expArg += log(fabs(X.k));
      
      expArg += log_middle_term;

      // We already know that the value of k is all right.
      temp = log(fabs(k));

      // allow for any possiblity in the comparison of temp and expArg.
      // the goal is to evaluate (e^a + e^b) or equivalently, (1 + e^(a-b))*e^b
      // case 1: if e^a and e^b are within 15 orders of magnitude 
      // (computer precision) then we should go ahead and evaluate the term as
      // stated.
      // case 2 and 3: if the 'if' fell through, then we know that one of the 
      // terms is much larger than the other, and that we only have to include 
      // the larger term to the final answer. the reason for going through this
      // trouble is that if we evalate the final answer through cases 2 or 3, 
      // then we don't actually have to calculate the exp of a potentially very
      // large or very small number -- instead we immediately shift the 
      // magnitude into the c degree of freedom.
      // ln(1e15) = 34

      if (IeeeMath::isNaN(temp))
	{
	  cerr << "Error in QMCGreensRatioComponent::add()" << endl;
	  cerr << "Attempting to calculate exp(" << temp << ")" << endl;
	  exit(1);
	}
      else if(fabs(temp - expArg) < 34)
	{
	  if (temp > expArg)
	    {
	      temp -= expArg;
	      k = exp(temp);
	      k = signK < 0 ? -k : k;
	      k += signXk < 0 ? -1.0 : 1.0;
	      c += expArg;
	    }
	  else if (expArg > temp)
	    {
	      expArg -= temp;
	      k = exp(expArg);
	      k = signXk < 0 ? -k : k;
	      k += signK < 0 ? -1.0 : 1.0;
	      c += temp;
	    }
	}
      else if (expArg > temp)
	{
	  k = 1.0;
	  k = signXk  < 0 ? -k : k;
	  c += expArg;
	} 
      else 
	{
	  k = 1.0;
	  k = signK < 0 ? -k : k;
	  c += temp;
	}
    } 
  else 
    {
      // An exponent is calculated here.
      expArg += log_middle_term;

      if (IeeeMath::isNaN(expArg))
	{
	  cerr << "Error in QMCGreensRatioComponent::add()" << endl;
	  cerr << "Attempting to calculate exp(" << expArg << ")" << endl;
	  exit(1);
	}
      k += X.k * exp(expArg);
    }

  // Here we even out the values if necessary.

  if ( (c < -10.0) && (log(fabs(k)) > 10.0) )
    {
      c += log(fabs(k));
      k = k < 0 ? -1.0 : 1.0;
    }
  else if ( (c > 10.0) && (log(fabs(k)) < -10.0) )
    {
      c += log(fabs(k));
      k = k < 0 ? -1.0 : 1.0;
    }

  if ( (c < -10.0) && (b*log(fabs(a)) > 10.0) )
    {
      c += b*log(fabs(a));
      a = 1.0;
      b = 0.0;
    }
  else if ( (c > 10.0) && (b*log(fabs(a)) < -10.0) )
    {
      c += b*log(fabs(a));
      a = 1.0;
      b = 0.0;
    }

  return *this;
}

QMCGreensRatioComponent QMCGreensRatioComponent::operator + ( const QMCGreensRatioComponent & rhs ) const {
  QMCGreensRatioComponent result = *this;
  return result.add(rhs);
}

QMCGreensRatioComponent QMCGreensRatioComponent::operator - ( const QMCGreensRatioComponent & rhs ) const {
  QMCGreensRatioComponent result = *this;
  return result.add(rhs*-1.0);
}

QMCGreensRatioComponent QMCGreensRatioComponent::operator * ( const QMCGreensRatioComponent & rhs ) const {
  QMCGreensRatioComponent result = *this;
  return result.multiplyBy(rhs);
}

QMCGreensRatioComponent QMCGreensRatioComponent::operator / ( const QMCGreensRatioComponent & rhs ) const {
  QMCGreensRatioComponent result = *this;
  return result.divideBy(rhs);
}

void QMCGreensRatioComponent::operator += ( const QMCGreensRatioComponent & rhs ){
  add(rhs);
}

void QMCGreensRatioComponent::operator -= ( const QMCGreensRatioComponent & rhs ){
  add(rhs*-1.0);
}  

void QMCGreensRatioComponent::operator *= ( const QMCGreensRatioComponent & rhs ){
  multiplyBy(rhs);
}

void QMCGreensRatioComponent::operator /= ( const QMCGreensRatioComponent & rhs ){
  divideBy(rhs);
}

void QMCGreensRatioComponent::operator = ( const QMCGreensRatioComponent & rhs ){
  k = rhs.k;
  a = rhs.a;
  b = rhs.b;
  c = rhs.c;
}

QMCGreensRatioComponent QMCGreensRatioComponent::operator + ( const double & rhs ) const {
  QMCGreensRatioComponent result = *this;
  return result.add(QMCGreensRatioComponent(rhs));
}

QMCGreensRatioComponent QMCGreensRatioComponent::operator - ( const double & rhs ) const {
  QMCGreensRatioComponent result = *this;
  return result.add(QMCGreensRatioComponent(-1.0*rhs));
}

QMCGreensRatioComponent QMCGreensRatioComponent::operator * ( const double & rhs ) const {
  QMCGreensRatioComponent result = *this;
  return result.multiplyBy(QMCGreensRatioComponent(rhs));
}

QMCGreensRatioComponent QMCGreensRatioComponent::operator / ( const double & rhs ) const {
  QMCGreensRatioComponent result = *this;
  return result.divideBy(QMCGreensRatioComponent(rhs));
}

void QMCGreensRatioComponent::operator += ( const double & rhs ){
  add(QMCGreensRatioComponent(rhs));
}

void QMCGreensRatioComponent::operator -= ( const double & rhs ){
  add(QMCGreensRatioComponent(-1.0*rhs));
}

void QMCGreensRatioComponent::operator *= ( const double & rhs ){
  multiplyBy(QMCGreensRatioComponent(rhs));
}

void QMCGreensRatioComponent::operator /= ( const double & rhs ){
  divideBy(QMCGreensRatioComponent(rhs));
}

void QMCGreensRatioComponent::operator = ( double rhs ){
  operator = (QMCGreensRatioComponent(rhs));
}

QMCGreensRatioComponent::operator double() const {
  return getValue();
}





