#include "QMCDouble.h"
#include <iomanip>

QMCDouble::QMCDouble()
{
  initialize();
}

QMCDouble::QMCDouble(double value)
{
  initialize();
  k = value;
}

QMCDouble::QMCDouble(double w, double x, double y,\
		     double z) 
{
  initialize();
  k = w;
  a = x;
  b = y;
  c = z;
}

QMCDouble::QMCDouble( const QMCDouble & rhs )
{
  *this = rhs;
}

QMCDouble::~QMCDouble()
{
  // does nothing
}

void QMCDouble::initialize()
{
  k = 1.0;
  a = 1.0;
  b = 0.0;
  c = 0.0;
}

ostream& operator << (ostream& strm, const QMCDouble &rhs)
{
  int w = 15;
  strm.setf(ios::scientific);
  //*
  //k (a^b) Exp[c]
  strm << setw(w) << rhs.k;
  strm << " ( " << setw(w) << rhs.a;
  strm << " ^ " << setw(w) << rhs.b;
  strm << " ) Exp[ " << setw(w) << rhs.c << " ]";
  /*/
  strm << "k: " << setw(w) << rhs.k;
  strm << " a: " << setw(w) << rhs.a;
  strm << " b: " << setw(w) << rhs.b;
  strm << " c: " << setw(w) << rhs.c;
  //*/
  return strm;
}

void QMCDouble::toXML(ostream & strm)
{
  strm << "<QMCDouble>" << endl;
  strm << "\t<k>\t" << k << "\t</k>" << endl;
  strm << "\t<a>\t" << a << "\t</a>" << endl;
  strm << "\t<b>\t" << b << "\t</b>" << endl;
  strm << "\t<c>\t" << c << "\t</c>" << endl;
  strm << "</QMCDouble>\n" << endl;
}

// ra^rb = na^nb/da^db
void QMCDouble::SimplifyRatioPowers(double na, double nb,
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
      cerr << "Error in QMCDouble::SimplifyRatioPowers()";
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

double QMCDouble::getValue() const
{
  if (a<0.0)
    {
      cerr << "Error in QMCDouble::getValue():" << endl;
      cerr << "attempting to take a power of a negative number." << endl;
      cerr << "value = " << *this << endl;
      exit(1);
    }      
  double t1 = 1.0;
  double t2 = 1.0;

  if(b != 0.0) t1 = pow(a,b);
  if(c > 500)  t2 = 1e250;
  else if(c != 0.0) t2 = exp(c);

  return k*t1*t2;
}

bool QMCDouble::isNotValid() const
{
  //k a^b exp(c)  
  if(IeeeMath::isNaN(k) || IeeeMath::isNaN(a) || IeeeMath::isNaN(b) || IeeeMath::isNaN(c))
    {
      return true;
    }
#ifdef QMC_DEBUG
  if(a < 0.0)
    {
      /*
	Unless b is an integer, a^b is imaginary.
	It's almost certainly a programming error, so let's
	not waste time explicitly checking it.
      */
      return true;
    }
#endif
  return false;
}

bool QMCDouble::isZero() const
{
  if(fabs(k) <  1e-250 || c < -690.0 )
    {
      return true;
    }
  return false;
}

QMCDouble & QMCDouble::divideBy(const QMCDouble & denom)
{
  if(isZero())
    return *this;

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

  if (isNotValid())
    {
      cerr << "Error in QMCDouble::divideBy():" << endl;
      cerr << "dividend = " << *this << endl;
      cerr << " divisor = " << denom << endl;
      exit(1);
    }      

  return *this;
}

QMCDouble & QMCDouble::multiplyBy(const QMCDouble &X)
{
  if (isNotValid() || X.isNotValid())
    {
      cerr << "Error in QMCDouble::multiplyBy():" << endl;
      cerr << "multiplicand = " << *this << endl;
      cerr << "  multiplier = " << X << endl;
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
  if( (fabs(temp_k) > 1e200 || fabs(temp_k) < 1e-200) &&
      (fabs(k) > 1e-300 && fabs(X.k) > 1e-300) )
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

  if (isNotValid() || X.isNotValid())
    {
      cerr << "Error in QMCDouble::multiplyBy():" << endl;
      cerr << "multiplicand = " << *this << endl;
      cerr << "  multiplier = " << X << endl;
      *this = 0.0;
    }
  
  return *this;
}

QMCDouble & QMCDouble::add(const QMCDouble &X)
{
  // This function has been causing a lot of problems on QSC.
  // First we make sure that none of the elements is NaN.

  if (isNotValid() || X.isNotValid())
    {
      cerr << "Error in QMCDouble::add():" << endl;
      cerr << "augend = " << *this << endl;
      cerr << "addend = " << X << endl;
      exit(1);
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
	  cerr << "Error in QMCDouble::add()" << endl;
	  cerr << "Attempting to calculate exp(" << temp << ")" << endl;
	  cerr << "augend = " << *this << endl;
	  cerr << "addend = " << X << endl;
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
	  cerr << "Error in QMCDouble::add()" << endl;
	  cerr << "Attempting to calculate exp(" << expArg << ")" << endl;
	  cerr << "augend = " << *this << endl;
	  cerr << "addend = " << X << endl;
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

QMCDouble QMCDouble::operator + ( const QMCDouble & rhs ) const {
  QMCDouble result = *this;
  return result.add(rhs);
}

QMCDouble QMCDouble::operator - ( const QMCDouble & rhs ) const {
  QMCDouble result = *this;
  return result.add(rhs*-1.0);
}

QMCDouble QMCDouble::operator * ( const QMCDouble & rhs ) const {
  QMCDouble result = *this;
  return result.multiplyBy(rhs);
}

QMCDouble QMCDouble::operator / ( const QMCDouble & rhs ) const {
  QMCDouble result = *this;
  return result.divideBy(rhs);
}

void QMCDouble::operator += ( const QMCDouble & rhs ){
  add(rhs);
}

void QMCDouble::operator -= ( const QMCDouble & rhs ){
  add(rhs*-1.0);
}  

void QMCDouble::operator *= ( const QMCDouble & rhs ){
  multiplyBy(rhs);
}

void QMCDouble::operator /= ( const QMCDouble & rhs ){
  divideBy(rhs);
}

void QMCDouble::operator = ( const QMCDouble & rhs ){
  k = rhs.k;
  a = rhs.a;
  b = rhs.b;
  c = rhs.c;
}

QMCDouble QMCDouble::operator + ( const double & rhs ) const {
  QMCDouble result = *this;
  return result.add(QMCDouble(rhs));
}

QMCDouble QMCDouble::operator - ( const double & rhs ) const {
  QMCDouble result = *this;
  return result.add(QMCDouble(-1.0*rhs));
}

QMCDouble QMCDouble::operator * ( const double & rhs ) const {
  QMCDouble result = *this;
  return result.multiplyBy(QMCDouble(rhs));
}

QMCDouble QMCDouble::operator / ( const double & rhs ) const {
  QMCDouble result = *this;
  return result.divideBy(QMCDouble(rhs));
}

void QMCDouble::operator += ( const double & rhs ){
  add(QMCDouble(rhs));
}

void QMCDouble::operator -= ( const double & rhs ){
  add(QMCDouble(-1.0*rhs));
}

void QMCDouble::operator *= ( const double & rhs ){
  multiplyBy(QMCDouble(rhs));
}

void QMCDouble::operator /= ( const double & rhs ){
  divideBy(QMCDouble(rhs));
}

void QMCDouble::operator = ( double rhs ){
  operator = (QMCDouble(rhs));
}

QMCDouble::operator double() const {
  return getValue();
}





