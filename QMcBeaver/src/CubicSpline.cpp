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

#include "CubicSpline.h"

CubicSpline::CubicSpline(){;}


void CubicSpline::operator=( const CubicSpline & rhs )
{
  yp0=rhs.yp0;
  ypend=rhs.ypend;
  n=rhs.n;
  f=rhs.f;
  dfdx=rhs.dfdx;
  ddfddx=rhs.ddfddx;

  x_list=rhs.x_list;
  y_list=rhs.y_list;
  yp_list=rhs.yp_list;
  a0_list=rhs.a0_list;
  a1_list=rhs.a1_list;
  a2_list=rhs.a2_list;
  a3_list=rhs.a3_list;
}


void CubicSpline::initializeWithFunctionValues(Array1D<double> &x_input, 
                Array1D<double> &y_input, double y_prime_0, double y_prime_end)
{
  //initialize the spline tables of values

  n=x_input.dim1();

  x_list.allocate(n);
  y_list.allocate(n);
  yp_list.allocate(n);
  a0_list.allocate(n-1);
  a1_list.allocate(n-1);
  a2_list.allocate(n-1);
  a3_list.allocate(n-1);

  for(int i=0;i<n;i++)
    {
      x_list(i)=x_input(i);
      y_list(i)=y_input(i);
      yp_list(i)=0.0;
    }

  for(int i=0;i<n-1;i++)
    {
      a0_list(i)=0.0;
      a1_list(i)=0.0;
      a2_list(i)=0.0;
      a3_list(i)=0.0;
    }

  //set BC's 

  yp0   = y_prime_0;
  ypend = y_prime_end;

  //if we look at the intervals xj-1,xj,xj+1
  //xj-1 <---h---> xj <---d---> xj+1 are the spacings

  //The equation we will use to get the yp_list from the y_list values is:
  //3((fj+1)/d^2+(fj)(1/h^2-1/d^2)-(fj-1)/h^2)
  //   =(f'j+1)/d+(f'j)(2/h+2/d)+(f'j-1)/h

  //Make tridiagonal matrix

  double h,d,h2,d2;
  row_right.allocate(n);
  row_middle.allocate(n);
  row_left.allocate(n);
  col_rhs.allocate(n);

  for(int i=1;i<n-1;i++)
    {
      h=x_list(i)-x_list(i-1);
      d=x_list(i+1)-x_list(i);

      row_right(i)    = 1/d;
      row_middle(i) = 2/h+2/d;
      row_left(i) = 1/h;
    }

  //do the edges of the matrix

  row_right(0)    = 0;
  row_middle(0)   = 1;
  row_left(0)     = 0;
  row_right(n-1)  = 0;
  row_middle(n-1) = 1;
  row_left(n-1)   = 0;

  //fill in the right hand side

  for(int i=1;i<n-1;i++)
    {
      h=x_list(i)-x_list(i-1);
      d=x_list(i+1)-x_list(i);
      h2=h*h;
      d2=d*d;

      col_rhs(i) = 3.0*(  y_list(i+1)/d2 
			       + y_list(i)  *(1/h2-1/d2) 
			       - y_list(i-1)/h2);
    }

  //do the edges of the matrix

  col_rhs(0)   = yp0;
  col_rhs(n-1) = ypend;

  //Solve tridiagonal matrix for the yp_list values

  solve_tridiagonal_system();

  //put the resulting y primes into the yp_list

  for(int i=0;i<n;i++)
    {
      yp_list(i)=col_rhs(i);
    }

  //The spline function, s(x) has this form in each interval:
  //sj(x)=a0+a1(x-xj)+a2(x-xj)^2+a3(x-xj)^3
  //with a different set of a0,...,a3 coeff for each interval

  //To get the coeff lists values we need to evaluate the following
  //a0=fj
  //a1=f'j
  //a2=3/d^2(fj+1-fj)-1/d((f'j+1)+2(f'j))
  //a3=2/d^3(fj-fj+1)+1/d^2((f'j)+(f'j+1))

  double fj,fj1,dfj,dfj1;

  for(int i=0;i<n-1;i++)
    {
      d    = x_list(i+1)-x_list(i);
      fj   = y_list(i);
      fj1  = y_list(i+1);
      dfj  = yp_list(i);
      dfj1 = yp_list(i+1);
      a0_list(i)=fj;
      a1_list(i)=dfj;
      a2_list(i)=3/d/d*(fj1-fj)-1/d*(dfj1+2*dfj);
      a3_list(i)=2/d/d/d*(fj-fj1)+1/d/d*(dfj+dfj1);
    }
}


//given the dervs at each point and an initial value at y_0 we find the spline
//with the last two function values equal
void CubicSpline::initializeWithDerivativeValues(Array1D<double> &x_input, 
				        Array1D<double> &yp_input, double y_0)
{
  //initialize the spline tables of values

  n=x_input.dim1();

  x_list.allocate(n);
  y_list.allocate(n);
  yp_list.allocate(n);
  a0_list.allocate(n-1);
  a1_list.allocate(n-1);
  a2_list.allocate(n-1);
  a3_list.allocate(n-1);

  for(int i=0;i<n;i++)
    {
      x_list(i)=x_input(i);
      yp_list(i)=yp_input(i);
      y_list(i)=0.0;
    }

  for(int i=0;i<n-1;i++)
    {
      a0_list(i)=0.0;
      a1_list(i)=0.0;
      a2_list(i)=0.0;
      a3_list(i)=0.0;
    }

  //set BC's 

  y0   = y_0;

  //if we look at the intervals xj-1,xj,xj+1
  //xj-1 <---h---> xj <---d---> xj+1 are the spacings

  //The equation we will use to get the y_list from the yp_list values is:
  //3((fj+1)/d^2+(fj)(1/h^2-1/d^2)-(fj-1)/h^2)
  //   =(f'j+1)/d+(f'j)(2/h+2/d)+(f'j-1)/h

  //Make tridiagonal matrix

  double h,d,h2,d2;
  row_right.allocate(n);
  row_middle.allocate(n);
  row_left.allocate(n);
  col_rhs.allocate(n);

  for(int i=1;i<n-1;i++)
    {
      h=x_list(i)-x_list(i-1);
      d=x_list(i+1)-x_list(i);
      h2=h*h;
      d2=d*d;

      row_right(i)    = 3/d2;
      row_middle(i)   = 3/h2-3/d2;
      row_left(i)     = -3/h2;
    }

  //do the edges of the matrix

  row_right(0)    = 0;
  row_middle(0)   = 1;
  row_left(0)     = 0;
  row_right(n-1)  = 0;
  row_middle(n-1) = 1;
  row_left(n-1)   = -1;

  //fill in the right hand side

  for(int i=1;i<n-1;i++)
    {
      h=x_list(i)-x_list(i-1);
      d=x_list(i+1)-x_list(i);

      col_rhs(i) = 3.0*(  yp_list(i+1)/d 
			       + yp_list(i)  *(2/h+2/d) 
			       + yp_list(i-1)/h);
    }

  //do the edges of the matrix

  col_rhs(0)   = y0;
  col_rhs(n-1) = 0.0;

  //Solve tridiagonal matrix for the y_list values

  solve_tridiagonal_system2();

  //put the resulting y's into the y_list

  for(int i=0;i<n;i++)
    {
      y_list(i)=col_rhs(i);
    }

  //The spline function, s(x) has this form in each interval:
  //sj(x)=a0+a1(x-xj)+a2(x-xj)^2+a3(x-xj)^3
  //with a different set of a0,...,a3 coeff for each interval

  //To get the coeff lists values we need to evaluate the following
  //a0=fj
  //a1=f'j
  //a2=3/d^2(fj+1-fj)-1/d((f'j+1)+2(f'j))
  //a3=2/d^3(fj-fj+1)+1/d^2((f'j)+(f'j+1))

  double fj,fj1,dfj,dfj1;

  for(int i=0;i<n-1;i++)
    {
      d    = x_list(i+1)-x_list(i);
      fj   = y_list(i);
      fj1  = y_list(i+1);
      dfj  = yp_list(i);
      dfj1 = yp_list(i+1);
      a0_list(i)=fj;
      a1_list(i)=dfj;
      a2_list(i)=3/d/d*(fj1-fj)-1/d*(dfj1+2*dfj);
      a3_list(i)=2/d/d/d*(fj-fj1)+1/d/d*(dfj+dfj1);
    }
}

//not stable if one of the middle row elements starts out as zero
void CubicSpline::solve_tridiagonal_system()
{  
  //forward part
  //eliminate row_left

  double factor;

  for(int i=0;i<=(n-2);i++)
    {
      factor = row_left(i+1)/row_middle(i);
      row_middle(i+1) = row_middle(i+1) - 
	row_right(i)*factor;
      col_rhs(i+1) = col_rhs(i+1) - col_rhs(i)*factor;
      row_left(i+1)   = 0.0;
    }

  //backwards part
  //eliminate the row_right

  for(int i=n-1;i>=1;i--)
    {
      factor = row_right(i-1)/row_middle(i);
      col_rhs(i-1) = col_rhs(i-1) - col_rhs(i)*factor;
      row_right(i-1)   = 0.0;
    }

  //make sure everything is normalized

  for(int i=0;i<n;i++)
    {
      col_rhs(i)    = col_rhs(i)/row_middle(i);
      row_middle(i) = 1.0;
    }
}


//stable if one of the middle row elements starts out as zero and
//the last row has 0,0,0,0,...0,x,y,0
void CubicSpline::solve_tridiagonal_system2()
{  
  //backward part
  //eliminate row_right

  double factor,b,c,d,e,f,g;

  for(int i=n-1;i>=1;i--)
    {
      b=row_middle(i-1);
      c=row_right(i-1);
      d=col_rhs(i-1);
      e=row_left(i);
      f=row_middle(i);
      g=col_rhs(i);
      
      factor = c/f;
      row_middle(i-1) = b - e*factor;
      col_rhs(i-1)    = d - g*factor;
      row_right(i-1)  = 0.0;
    }

  //forwards part
  //eliminate the row_left

  for(int i=0;i<=n-2;i++)
    {
      b=row_middle(i);
      d=col_rhs(i);
      e=row_left(i+1);
      f=row_middle(i+1);
      g=col_rhs(i+1);
      factor = e/b;

      col_rhs(i+1)  = g-d*factor;
      row_left(i+1) = 0.0;
    }

  //make sure everything is normalized

  for(int i=0;i<n;i++)
    {
      col_rhs(i)    = col_rhs(i)/row_middle(i);
      row_middle(i) = 1.0;
    }
}

void CubicSpline::evaluate(double x)
{
  //find which index in the table by means of bisection
  int klo,khi,k;

  //push any value outside the spline range into the last value of the spline

  if(x>=x_list(n-1))
    {
      x=x_list(n-1);
    }
  else if(x<=x_list(0))
    {
      x=x_list(0);
    }

  klo=0;
  khi=n-1;

  while (khi-klo > 1) 
    {
      k=(khi+klo) >> 1;
      if (x_list(k) > x) khi=k;
      else klo=k;
    }

  evaluate(x,klo);
}

void CubicSpline::evaluate(double x, int index)
{
  int klo = index;
  int khi = index + 1;

  if(khi>n || khi<1 || klo>(n-1) || klo<0)
    {
      cerr << "Bad value for index in CubicSpline::evaluate(double x, "
	   << "int index)" << endl;
      cerr << "index= " << index << endl;
      cerr << "x= " << x << endl;
      exit(1);
    }

  double r,r2,a0,a1,a2,a3;

  r=x-x_list(klo);
  r2=r*r;
  a0=a0_list(klo);
  a1=a1_list(klo);
  a2=a2_list(klo);
  a3=a3_list(klo);

  f = a0+a1*r+a2*r2+a3*r2*r;

  dfdx = a1+2.0*a2*r+3.0*a3*r2;

  ddfddx = 2.0*a2+6.0*a3*r;
}

double CubicSpline::getFunctionValue()
{
  return f;
}

double CubicSpline::getFirstDerivativeValue()
{
  return dfdx;
}
 
double CubicSpline::getSecondDerivativeValue()
{
  return ddfddx;
}

void CubicSpline::toXML(ostream& strm)
{
  strm << "<CubicSpline>" << endl;
  strm << "<x_list>" << endl;
  for(int i=0;i<x_list.dim1();i++)
    {
      strm << x_list(i) << "\t";  
    }
  strm << "</x_list>" << endl;

  strm << "<y_list>" << endl;
  for(int i=0;i<y_list.dim1();i++)
    {
      strm << y_list(i) << "\t";  
    }
  strm << "</y_list>" << endl;

  strm << "<yp_list>" << endl;
  for(int i=0;i<yp_list.dim1();i++)
    {
      strm << yp_list(i) << "\t";  
    }
  strm << "</yp_list>" << endl;

  strm << "<a0_list>" << endl;
  for(int i=0;i<a0_list.dim1();i++)
    {
      strm << a0_list(i) << "\t";  
    }
  strm << "</a0_list>" << endl;

  strm << "<a1_list>" << endl;
  for(int i=0;i<a1_list.dim1();i++)
    {
      strm << a1_list(i) << "\t";  
    }
  strm << "</a1_list>" << endl;

  strm << "<a2_list>" << endl;
  for(int i=0;i<a2_list.dim1();i++)
    {
      strm << a2_list(i) << "\t";  
    }
  strm << "</a2_list>" << endl;

  strm << "<a3_list>" << endl;
  for(int i=0;i<a3_list.dim1();i++)
    {
      strm << a3_list(i) << "\t";  
    }
  strm << "</a3_list>" << endl;

  strm << "<f>\t" << f << "\t</f>" << endl;
  strm << "<dfdx>\t" << dfdx << "\t</dfdx>" << endl;
  strm << "<ddfddx>\t" << ddfddx << "\t</ddfddx>" << endl;
  strm << "</CubicSpline>" << endl;
}

ostream& operator<<(ostream & strm, CubicSpline & rhs)
{
  rhs.toXML(strm);
  return strm;
}


