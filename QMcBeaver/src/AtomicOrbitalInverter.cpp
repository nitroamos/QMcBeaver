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

#include "AtomicOrbitalInverter.h"

AtomicOrbitalInverter::AtomicOrbitalInverter()
{
  cutoff1  = 1.0;  //inner high density of points
  cutoff2  = 4.0;  //medium distance
  cutoff3  = 15.0; //long range
  npoints_r     = 1000;
  npoints_theta = 1000;
  npoints_phi   = 1000;
}
 
void AtomicOrbitalInverter::operator=( const AtomicOrbitalInverter & rhs )
{
  xn      = rhs.xn;
  yn      = rhs.yn;
  zn      = rhs.zn;
  cutoff1 = rhs.cutoff1;
  cutoff2 = rhs.cutoff2;
  cutoff3 = rhs.cutoff3;
  npoints_r     = rhs.npoints_r;
  npoints_theta = rhs.npoints_theta;
  npoints_phi   = rhs.npoints_phi;
  r_form     = rhs.r_form;
  theta_form = rhs.theta_form;
  phi_form   = rhs.phi_form;
}


void AtomicOrbitalInverter::set_npoints(int points_r,int points_theta,int points_phi)
{
  npoints_r     = points_r;
  npoints_theta = points_theta;
  npoints_phi   = points_phi;
}

void AtomicOrbitalInverter::initialize(int Xn, int Yn, int Zn,Array1D<double> B,Array1D<double> C)
{
  xn = Xn;
  yn = Yn;
  zn = Zn;
  
  b.allocate(B.dim1());
  c.allocate(C.dim1());

  for(int i=0;i<B.dim1();i++)
    {
      b(i)=B(i);
      c(i)=C(i);
    }

  Array1D<double> r_inp;
  Array1D<double> R_inp;
  Array1D<double> theta_inp;
  Array1D<double> THETA_inp;
  Array1D<double> phi_inp;
  Array1D<double> PHI_inp;
  r_inp.allocate(npoints_r);
  R_inp.allocate(npoints_r);
  theta_inp.allocate(npoints_theta);
  THETA_inp.allocate(npoints_theta);
  phi_inp.allocate(npoints_phi);
  PHI_inp.allocate(npoints_phi);
  
  double rend1,rend2,rend3,dr1,dr2,dr3,dtheta,dphi,r,f;
  r = 0.0;
  
  rend1 =  cutoff1;
  rend2 =  cutoff2;
  rend3 =  cutoff3;
  
  int np1,np2,np3;
  np1=npoints_r/3;
  np2=npoints_r/3;
  np3=npoints_r-np1-np2;

  dr1=(rend1)/np1;
  dr2=(rend2-rend1)/np2;
  dr3=(rend3-rend2)/np3;
  dtheta=PI/npoints_theta;
  dphi=2.0*PI/npoints_phi;
  
  //  if((xn+yn+zn)>0)  
  if(true)   //just making the spline no matter what right now
    {
      r_inp(0)=0.0;
      for(int i=1;i<npoints_r;i++)
	{
	  if(r_inp(i-1)<rend1)
	    {
	      //left part of region 3
	      r=r_inp(i-1)+dr1;
	    }
	  else if(r_inp(i-1)<rend2)
	    {
	      //left part of region 2
	      r=r_inp(i-1)+dr2;
	    }
	  else if(r_inp(i-1)<rend3)
	    {
	      //region 1
	      r=r_inp(i-1)+dr3;
	    }
	  
	  r_inp(i)=r;
	  
	  //initialize r part
	  f=pow(r,xn+yn+zn)*eval_gaussians(r);
	  f=f*f*r*r;    //extra r^2 is for spherical volume element
	  R_inp(i)=f;
	}
      
      for(int i=0;i<npoints_theta;i++)
	{
	  theta_inp(i)=i*dtheta;
	  THETA_inp(i)=pow(sin(theta_inp(i)),2*xn+2*yn)*pow(cos(theta_inp(i)),2*zn)*sin(theta_inp(i)); //extra sin(theta) is for spherical volume element
	}
      
      //if((xn+yn)>0)
      if(true)   //just making the spline no matter what right now
	{
	  for(int i=0;i<npoints_phi;i++)
	    {
	      phi_inp(i)=i*dphi;
	      PHI_inp(i)=pow(cos(phi_inp(i)),2*xn)*pow(sin(phi_inp(i)),2*yn);
	    }
	}
      
      //initialize r part
      r_form.initialize(r_inp,R_inp);
      //initialize theta part
      theta_form.initialize(theta_inp,THETA_inp);
      
      //if((xn+yn)>0)
      if(true)   //just making the spline no matter what right now
	{
	  //initialize phi part
	  phi_form.initialize(phi_inp,PHI_inp);
	}
    }
}


//NOT DONE
//void AtomicOrbitalInverter::initialize(QMCBasisFunctionCoefficients &BSC)
//{
  //Fill in the coefficients into the initialization form below
  //initialize(Xn, Yn, Zn, B, C);
//}
//END NOT DONE

void AtomicOrbitalInverter::get_xyz(long &iseed, double &x, double &y, double &z)
{
  double r,theta,phi;
  //if((xn+yn+zn)>0)
  if(true)   //just making the spline no matter what right now
    {
      r     = r_form.random(iseed);
      theta = theta_form.random(iseed);
      //if((xn+yn)>0)
      if(true)   //just making the spline no matter what right now
	{
	  phi = phi_form.random(iseed);
	}
      else
	{
	  //phi is uniform
	  phi=ran1(&iseed)*2.0*PI;
	}
    }
  else
    {
      //theta, phi, and R are uniform so just invert the gaussians
      r     = invert_gaussians(iseed);
      theta = ran1(&iseed)*PI;
      phi   = ran1(&iseed)*2.0*PI;
    }

  x = r * sin(theta)*cos(phi);
  y = r * sin(theta)*sin(phi);
  z = r * cos(theta);
}

double AtomicOrbitalInverter::eval_gaussians(double r)
{
  double sum=0.0;
  for(int i=0;i<b.dim1();i++)
    {
      sum = sum + c(i)*exp(-b(i)*r*r);
    }
  return sum;
}


//not used in the code on 10/30/01
double AtomicOrbitalInverter::invert_gaussians(long & iseed)
{
  //analytically invert the sum of the gaussians squared (another sum of g's)
  
  //sum is the integral of the atomic orbital squared 
  double sum_a=0.0;
  for(int i=0;i<b.dim1();i++)
    {
      for(int j=0;j<b.dim1();j++)
	{
	  sum_a=sum_a+c(i)*c(j)*sqrt(PI/(b(i)+b(j)));
	}
    }
  double R=ran1(&iseed);
  double sum_b=0.0;
  int I,J;
  I = 0;
  J = 0;
  for(int i=0;i<b.dim1();i++)
    {
      for(int j=0;j<b.dim1();j++)
	{
	  sum_b=sum_b+c(i)*c(j)*sqrt(PI/(b(i)+b(j)))/sum_a;
	  if(sum_b>R)
	    {
	      I=i;
	      J=j;
	      i=b.dim1();
	      j=b.dim1();
	    }
	}
    }
  
  double b_new=b(I)+b(J);
  double sigma=sqrt(1/2.0/b_new);
  return fabs(gasdev(&iseed)*sigma);
}
