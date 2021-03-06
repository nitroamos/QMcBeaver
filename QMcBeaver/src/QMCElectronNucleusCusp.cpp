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

#include "QMCElectronNucleusCusp.h"

QMCElectronNucleusCusp::QMCElectronNucleusCusp()
{
}

void QMCElectronNucleusCusp::operator=(const QMCElectronNucleusCusp& rhs)
{
  Input = rhs.Input;
  Molecule = rhs.Molecule;
  BF = rhs.BF;

  ORWF_coeffs = rhs.ORWF_coeffs;

  natoms = rhs.natoms;
  norbitals = rhs.norbitals;

  ORParams = rhs.ORParams;
}

void QMCElectronNucleusCusp::initialize(QMCInput* input, const Array2D<qmcfloat>& WFCoeffs)
{
  Input = input;
  Molecule = &Input->Molecule;
  BF = &Input->BF;
  ORWF_coeffs = WFCoeffs;

  natoms = Molecule->getNumberAtoms();
  norbitals = ORWF_coeffs.dim1();

  ORParams.allocate(natoms,norbitals);
}

void QMCElectronNucleusCusp::fitReplacementOrbitals()
{
  int BFindex = 0;
  int nSGaussians = 0;
  int counter = 0;
  double xnuc,ynuc,znuc;
  double xrel,yrel,zrel,rsqrel,xyz;
  int a,b,c;
  double temp,p0,p1;
  bool repOrbitalNeeded;

  for (int i=0; i<natoms; i++)
    {
      BF_coeffs = BF->getBFCoeffs(i);

      // We find how many s Gaussians are centered on this nucleus.
      nSGaussians = 0;
      for (int k=0; k<BF_coeffs->getNumberBasisFunctions(); k++)
	if (BF_coeffs->Type(k) == "s")
	  nSGaussians += BF_coeffs->N_Gauss(k);

      sTypeCoeffs.allocate(nSGaussians,2);

      // If there are no s Gaussians centered on this nucleus, we can skip it.
      if (nSGaussians == 0)
	{
	  QMCElectronNucleusCuspParameters TRASH;
	  TRASH.set_rc(0.0);
	  for (int m=0; m<norbitals; m++)
	    ORParams(i,m) = TRASH;
	  continue;
	}

      // Then we find where the orbital coeffs for this atom start.
      BFindex = 0;
      for (int j=0; j<i; j++)
	BFindex += BF->getNumberBasisFunctions(j);

      Z = Molecule->Z(i);

      xnuc = Molecule->Atom_Positions(i,0);
      ynuc = Molecule->Atom_Positions(i,1);
      znuc = Molecule->Atom_Positions(i,2);

      for (int m=0; m<norbitals; m++)
	{
	  counter = 0;
	  BF_coeffs = BF->getBFCoeffs(i);
	  for (int k=0; k<BF_coeffs->getNumberBasisFunctions(); k++)
	    if (BF_coeffs->Type(k) == "s")
	      for (int l=0; l<BF_coeffs->N_Gauss(k); l++)
		{
		  sTypeCoeffs(counter,0) = BF_coeffs->Coeffs(k,l,0);
		  sTypeCoeffs(counter,1) = BF_coeffs->Coeffs(k,l,1) * 
		    ORWF_coeffs(m,BFindex+k);
		  counter++;
		}

	  // If all the expansion coefficients for the s Gaussians centered on
	  // this nucleus for this orbital are zero, there is no need for the 
	  // replacement orbital.
	  repOrbitalNeeded = false;
	  for (int k=0; k<nSGaussians; k++)
	    if ( fabs(sTypeCoeffs(k,1)) > 1e-15 )
	      {
		repOrbitalNeeded = true;
		break;
	      }

	  if (repOrbitalNeeded == false)
	    {
	      QMCElectronNucleusCuspParameters TRASH;
	      TRASH.set_rc(0.0);
	      ORParams(i,m) = TRASH;
	      continue;
	    }

	  // This function determines rc for this nucleus and orbital and fits
	  // the ideal curve.
	  determineRc();

	  // We calculate n0, which is the value of the orbital due to basis
	  // functions centered on other nuclei at this nucleus.
	  n0 = 0.0;
	  counter = 0;
	  for (int j=0; j<natoms; j++)
	    {
	      if (j != i)
		{
		  xrel = xnuc - Molecule->Atom_Positions(j,0);
		  yrel = ynuc - Molecule->Atom_Positions(j,1);
		  zrel = znuc - Molecule->Atom_Positions(j,2);
		  rsqrel = xrel*xrel + yrel*yrel + zrel*zrel;
		  BF_coeffs = BF->getBFCoeffs(j);
		  for (int k=0; k<BF_coeffs->getNumberBasisFunctions(); k++)
		    {
		      a = BF_coeffs->xyz_powers(k,0);
		      b = BF_coeffs->xyz_powers(k,1);
		      c = BF_coeffs->xyz_powers(k,2);
		      xyz = pow(xrel,a)*pow(yrel,b)*pow(zrel,c);
		      temp = 0.0;
		      for (int l=0; l<BF_coeffs->N_Gauss(k); l++)
			{
			  p0 = BF_coeffs->Coeffs(k,l,0);
			  p1 = BF_coeffs->Coeffs(k,l,1);
			  temp += p1*exp(-p0*rsqrel);
			}
		      n0 += temp*xyz*ORWF_coeffs(m,counter);
		      counter++;
		    }
		}
	      else if (j == i)
		counter += BF->getNumberBasisFunctions(i);
	    }

	  // This function fits the exponential function to the original
	  // orbital.  The deviation of the local energy of the replacement 
	  // orbital is optimized with respect to its value at the nucleus.
	  ORParams(i,m) = fitOrbitalParameters();
	  if (globalInput.flags.print_replacement_orbitals == 1)
	    {
	      clog << "Orbital replacement parameters for atom " << i;
	      clog << " orbital " << m << endl;
	      ORParams(i,m).printParameters();
	    }
	}
    }
}

void QMCElectronNucleusCusp::determineRc()
{
  // The maximum radius of correction is 1/Z.
  const double rcmax = 1.0/Z;

  // We fit the ideal curve at rcmax.
  if (Z > 1)
    {
      idealCurveParams.allocate(9);

      idealCurveParams(8) = 1.94692;
      idealCurveParams(7) = -12.1316;
      idealCurveParams(6) = 31.2276;
      idealCurveParams(5) = -42.8705;
      idealCurveParams(4) = 33.7308;
      idealCurveParams(3) = -15.0126;
      idealCurveParams(2) = 3.25819;
      idealCurveParams(1) = 0.0;
      idealCurveParams(0) = 0.0;

      idealCurve.initialize(idealCurveParams);
      idealCurve.evaluate(rcmax);
      idealCurveParams(0) = getOrigLocalEnergy(rcmax)/(Z*Z) - idealCurve.getFunctionValue();

      idealCurve.initialize(idealCurveParams);
    }

  else if (Z == 1)
    {
      idealCurveParams.allocate(1);
      idealCurveParams(0) = getOrigLocalEnergy(rcmax);
      idealCurve.initialize(idealCurveParams);
    }
  else
    {
      cerr << "ERROR in QMCElectronNucleusCusp: Z = " << Z << endl;
      exit(0);
    }

  // Now we find the roots of the orbital.  We evaluate the orbital at 100 
  // points between rcmax and 0 and see if it changes sign.
  const double dr = rcmax/100;
  double old_value = evaluateOrigOrbital(rcmax);
  double new_value = 0.0;
  int nRoots = 0;
  list<int> rtIndex;
  Array1D<double> roots;

  for (int i=0; i<101; i++)
    {
      new_value = evaluateOrigOrbital(rcmax-i*dr);
      if (new_value*old_value < 0)
	{
	  rtIndex.push_back(i);
	  nRoots++;
	}
      old_value = new_value;
    }

  // If there are any roots we find them in their intervals.  We use a
  // bisection method.
  if (nRoots > 0)
    {
      roots.allocate(nRoots);
      int counter = 0;

      double rl = 0.0;
      double rh = 0.0;
      double fl = 0.0;
      double fh = 0.0;
      double rm = 0.0;
      double fm = 0.0;

      const double tol = 1e-8;

      for (list<int>::iterator rp=rtIndex.begin(); rp!=rtIndex.end(); ++rp)
	{
	  rl = rcmax-(*rp)*dr;
	  fl = evaluateOrigOrbital(rl);
	  rh = rl+dr;
	  fh = evaluateOrigOrbital(rh);

	  if (fl*fh > 0)
	    {
	      cerr << "ERROR in QMCElectronNucleusCusp::determineRc(): ";
	      cerr << rl << " and " << rh << " do not bracket a root" << endl;
	      exit(0);
	    }
	  if (fabs(fl)<tol)
	    {
	      roots(counter) = rl;
	      counter++;
	      continue;
	    }
	  else if (fabs(fh)<tol)
	    {
	      roots(counter) = rh;
	      counter++;
	      continue;
	    }
	  else
	    {
	      for (int i=0; i<10000; i++)
		{
		  if (i == 9999)
		    {
		      cerr << "ERROR in ";
		      cerr << "QMCElectronNucleusCusp::determineRc(): ";
		      cerr << " too many iterations in finding root" << endl;
		      exit(0);
		    }
		  rm = (rl+rh)/2;
		  fm = evaluateOrigOrbital(rm);

		  if (fabs(fm)<tol)
		    {
		      roots(counter) = rm;
		      counter++;
		      break;
		    }
		  else if (fm*fl<0)
		    {
		      // rl and rm now bracket the root.
		      rh = rm;
		      fh = fm;
		    }
		  else
		    {
		      // rm and rh bracket the root.
		      rl = rm;
		      fl = fm;
		    }
		}
	      continue;
	    }
	}
    }

  // Now we have an array of all the roots, in order of decreasing r.
  // We want to look starting at rcmax and going closer, excluding regions 
  // around the nodes, for the point where the deviation from the ideal
  // curve gets large.

  double origGradientNearNucleus = getOrigGradient(rcmax/1000);
  double origValueAtNucleus = evaluateOrigOrbital(0.0);

  // The gradient of the original orbital at the nucleus is zero.  We find the 
  // value of the gradient near the nucleus.  If it is positive, we place rc in
  // a region where the original value has a greater value than it does at the
  // nucleus.  If the gradient is negative, we place rc in a region where the
  // original orbital has a value less than its value at the nucleus.

  // This is to fix a problem that came up in fitting some orbitals that were
  // decreasing from rc to the origin, but had a positive gradient at the
  // origin.  The fit function ended up having the wrong cusp condition.

  rc = rcmax;
  double r_trial = rcmax;
  int rootIndex = -1;
  double r_bound = 0.0;
  const double ddr = dr/10;
  const double maxdev = Z*Z/50.0;
  
  if (nRoots>0)
    {
      rootIndex = nRoots-1;
      r_bound = roots(rootIndex) + rcmax/5;
    }

  while (r_trial > 0.0)
    {
      if (origGradientNearNucleus*(evaluateOrigOrbital(r_trial)-origValueAtNucleus) > 0.0)
	{
	  if (r_trial > r_bound)
	    {
	      idealCurve.evaluate(r_trial);
	      old_value = getOrigLocalEnergy(r_trial);
	      if (fabs(Z*Z*idealCurve.getFunctionValue()-old_value) > maxdev)
		{
		  r_trial += dr;
		  while (r_trial > 0.0)
		    {
		      r_trial -= ddr;
		      if (origGradientNearNucleus*(evaluateOrigOrbital(r_trial)-origValueAtNucleus) > 0.0)
			{
			  idealCurve.evaluate(r_trial);
			  old_value = getOrigLocalEnergy(r_trial);
			  if (fabs(Z*Z*idealCurve.getFunctionValue()-old_value)>maxdev)
			    {
			      rc = r_trial;
			      break;
			    }
			}
		    }
		  break;
		}
	    }
	  else
	    {
	      idealCurve.evaluate(r_bound);
	      old_value = getOrigLocalEnergy(r_bound);
	      if ( (origGradientNearNucleus*(evaluateOrigOrbital(r_trial)-origValueAtNucleus) > 0.0) && (fabs(Z*Z*idealCurve.getFunctionValue()-old_value) > maxdev) )
		{
		  while (r_trial > 0.0)
		    {
		      r_trial -= ddr;
		      idealCurve.evaluate(r_trial);
		      old_value = getOrigLocalEnergy(r_trial);
		      if (fabs(Z*Z*idealCurve.getFunctionValue()-old_value)>maxdev)
			{
			  rc = r_trial;
			  break;
			}
		    }
		  break;
		}	      
	      else
		{
		  r_trial = roots(rootIndex) - rcmax/5 + dr;
		  if (rootIndex > 0)
		    {
		      rootIndex--;
		      r_bound = roots(rootIndex) + rcmax/5;
		    }
		  else if (rootIndex == 0)
		    r_bound = 0.0;
		}
	    }
	}
      r_trial -= dr;
    }

  // We fit the ideal curve at rc.
  if (Z > 1)
    {
      idealCurveParams(0) = 0.0;

      idealCurve.initialize(idealCurveParams);
      idealCurve.evaluate(rc);
      idealCurveParams(0) = getOrigLocalEnergy(rc)/(Z*Z) - idealCurve.getFunctionValue();

      idealCurve.initialize(idealCurveParams);
    }

  else if (Z == 1)
    {
      idealCurveParams.allocate(1);
      idealCurveParams(0) = getOrigLocalEnergy(rc);
      idealCurve.initialize(idealCurveParams);
    }
}

QMCElectronNucleusCuspParameters QMCElectronNucleusCusp::fitOrbitalParameters()
{
  // With rc, this function fits the orbital parameters.  phi0 is adjustable, 
  // and is optimized to minimize the max deviation of the replacement orbital
  // from the ideal curve.

  QMCElectronNucleusCuspParameters loParams;
  QMCElectronNucleusCuspParameters midParams;
  QMCElectronNucleusCuspParameters hiParams;

  loParams.initialize(sTypeCoeffs,n0,rc,idealCurve,Z);
  midParams.initialize(sTypeCoeffs,n0,rc,idealCurve,Z);
  hiParams.initialize(sTypeCoeffs,n0,rc,idealCurve,Z);

  // The value of phi0 will always have the same sign as the original orbital,
  // with equal or greater magnitude.
  double origAtOrigin = evaluateOrigOrbital(0.0);
  double origAtRc = evaluateOrigOrbital(rc);

  double lo_phi = 0.0;
  double mid_phi = 0.0;
  double hi_phi = 0.0;

  if (origAtOrigin-origAtRc < 0.0)
    {
      hi_phi = origAtOrigin;
      mid_phi = origAtOrigin - 0.05*fabs(origAtOrigin-origAtRc);
      lo_phi = origAtOrigin - 0.1*fabs(origAtOrigin-origAtRc);
    }
  else
    {
      lo_phi = origAtOrigin;
      mid_phi = origAtOrigin + 0.05*fabs(origAtOrigin-origAtRc);
      hi_phi = origAtOrigin + 0.1*fabs(origAtOrigin-origAtRc);
    }
  
  loParams.fitReplacementOrbital(lo_phi);
  midParams.fitReplacementOrbital(mid_phi);
  hiParams.fitReplacementOrbital(hi_phi);

  double loSigmaSq = loParams.getSigmaSq();
  double midSigmaSq = midParams.getSigmaSq();
  double hiSigmaSq = hiParams.getSigmaSq();
  
#define GOLD 0.618034

  // We will look for the minimum between hi_phi and lo_phi.
  // We are not going to search beyond 10% of the original value.
  // We don't want the cusp replacement to drastically change the original orbital.  
 
  while (hi_phi-lo_phi > fabs(mid_phi/10000))
    {
      if (loSigmaSq < midSigmaSq)
        {
          // Low is less than mid
          if (midSigmaSq < hiSigmaSq)
            {
              // Mid is less than high, we set high to mid and reset mid
              hi_phi = mid_phi;
              hiSigmaSq = midSigmaSq;
              hiParams = midParams;
 
              mid_phi = lo_phi+GOLD*(hi_phi-lo_phi);
              midParams.fitReplacementOrbital(mid_phi);
              midSigmaSq = midParams.getSigmaSq();
              continue;
            }
          else
            {
              // Mid is greater than high.
              if (loSigmaSq < hiSigmaSq)
                {
                  // Low is less than high, we set high to mid and reset mid
                  hi_phi = mid_phi;
                  hiSigmaSq = midSigmaSq;
                  hiParams = midParams;
 
                  mid_phi = lo_phi+GOLD*(hi_phi-lo_phi);
                  midParams.fitReplacementOrbital(mid_phi);
                  midSigmaSq = midParams.getSigmaSq();
                  continue;
                }
	      else
		{
		  // High is less than low, we set low to mid and reset mid
		  lo_phi = mid_phi;
		  loSigmaSq = midSigmaSq;
		  loParams = midParams;
 
		  mid_phi = hi_phi-GOLD*(hi_phi-lo_phi);
		  midParams.fitReplacementOrbital(mid_phi);
		  midSigmaSq = midParams.getSigmaSq();
		  continue;
		}
	    }
	}
      else
	{
	  // low is higher than mid
	  if (hiSigmaSq < midSigmaSq)
	    {
	      // mid is higher than hi, set lo to mid and reset mid
	      lo_phi = mid_phi;
	      loSigmaSq = midSigmaSq;
	      loParams = midParams;
 
	      mid_phi = hi_phi-GOLD*(hi_phi-lo_phi);
	      midParams.fitReplacementOrbital(mid_phi);
	      midSigmaSq = midParams.getSigmaSq();
	      continue;
	    }
	  else
	    {
	      // mid is lower than hi
	      if (loSigmaSq < hiSigmaSq)
		{
		  // set hi to mid, reset mid
		  hi_phi = mid_phi;
		  hiParams = midParams;
		  hiSigmaSq = midSigmaSq;
 
		  mid_phi = hi_phi - GOLD*(hi_phi-lo_phi);
		  midParams.fitReplacementOrbital(mid_phi);
		  midSigmaSq = midParams.getSigmaSq();
		  continue;
		}
	      else
		{
		  // set lo to mid, reset mid
		  lo_phi = mid_phi;
		  loSigmaSq = midSigmaSq;
		  loParams = midParams;

                  mid_phi = lo_phi + GOLD*(hi_phi-lo_phi);
                  midParams.fitReplacementOrbital(mid_phi);
                  midSigmaSq = midParams.getSigmaSq();
                  continue;
		}
	    }
	}
    }
#undef GOLD

  // The best set of parameters is returned.
  return midParams;
}

void QMCElectronNucleusCusp::replaceCusps(Array2D<double> & X, 
					  int Start, int Stop,
					  Array2D<qmcfloat> & D,
					  Array2D<qmcfloat> & GradX,
					  Array2D<qmcfloat> & GradY,
					  Array2D<qmcfloat> & GradZ,
					  Array2D<qmcfloat> & Laplacian_D)
{
  // For each nucleus, for each electron, calculates x,y,z,r.
  // For each orbital, checks if r<rc.  If so, replaces the appropriate element
  // of the Slater matrices.
  double xnuc,ynuc,znuc;
  double xrel,yrel,zrel,rrel;
  double temp_value,temp_gradx,temp_grady,temp_gradz,temp_lap;
  int counter = 0;

  if(Stop - Start + 1 != D.dim1())
    {
      clog << "Warning: dimensions don't match in QMCElectronNucleusCusp." << endl;  
      clog << " Start = " << Start << endl;
      clog << "  Stop = " << Stop << endl;
      clog << "D.dim1 = " << D.dim1() << endl;
    }

  for (int i=0; i<natoms; i++)
    {
      if(globalInput.Molecule.usesPseudo(i))
	continue;

      xnuc = Molecule->Atom_Positions(i,0);
      ynuc = Molecule->Atom_Positions(i,1);
      znuc = Molecule->Atom_Positions(i,2);
      
      counter = 0;
      for (int j=Start; j<=Stop; j++)
	{
	  xrel = X(j,0) - xnuc;
	  yrel = X(j,1) - ynuc;
	  zrel = X(j,2) - znuc;
	  rrel = sqrt(xrel*xrel + yrel*yrel + zrel*zrel);

	  for (int k=0; k<norbitals; k++)
	    {
	      if (rrel < ORParams(i,k).get_rc())
		{
		  temp_value = 0.0;
		  temp_gradx = 0.0;
		  temp_grady = 0.0;
		  temp_gradz = 0.0;
		  temp_lap = 0.0;

		  ORParams(i,k).replaceOrbitalValues(xrel,yrel,zrel,rrel,
			 temp_value,temp_gradx,temp_grady,temp_gradz,temp_lap);

		  D(counter,k) += temp_value;
		  GradX(counter,k) += temp_gradx;
		  GradY(counter,k) += temp_grady;
		  GradZ(counter,k) += temp_gradz;
		  Laplacian_D(counter,k) += temp_lap;
		}
	    }
	  counter++;
	}
    }
}

void QMCElectronNucleusCusp::replaceCusps(Array2D<double> & X, 
					  int Start, int Stop,
					  Array2D<qmcfloat> & D)
{
  // For each nucleus, for each electron, calculates x,y,z,r.
  // For each orbital, checks if r<rc.  If so, replaces the appropriate element
  // of the Slater matrices.
  double xnuc,ynuc,znuc;
  double xrel,yrel,zrel,rrel;
  double temp_value;
  int counter = 0;

  if(Stop - Start + 1 != D.dim1())
    {
      clog << "Warning: dimensions don't match in QMCElectronNucleusCusp." << endl;  
      clog << " Start = " << Start << endl;
      clog << "  Stop = " << Stop << endl;
      clog << "D.dim1 = " << D.dim1() << endl;
    }

  for (int i=0; i<natoms; i++)
    {
      if(globalInput.Molecule.usesPseudo(i))
	continue;

      xnuc = Molecule->Atom_Positions(i,0);
      ynuc = Molecule->Atom_Positions(i,1);
      znuc = Molecule->Atom_Positions(i,2);
      
      counter = 0;
      for (int j=Start; j<=Stop; j++)
	{
	  xrel = X(j,0) - xnuc;
	  yrel = X(j,1) - ynuc;
	  zrel = X(j,2) - znuc;
	  rrel = sqrt(xrel*xrel + yrel*yrel + zrel*zrel);

	  for (int k=0; k<norbitals; k++)
	    {
	      if (rrel < ORParams(i,k).get_rc())
		{
		  temp_value = 0.0;

		  ORParams(i,k).replaceOrbitalValues(xrel,yrel,zrel,rrel,
						     temp_value);
		  
		  D(counter,k) += temp_value;
		}
	    }
	  counter++;
	}
    }
}

double QMCElectronNucleusCusp::getOrigLocalEnergy(double r_orig)
{
  double r_sq = r_orig*r_orig;
  int nGaussians = sTypeCoeffs.dim1();

  double term_value = 0.0;
  double p0 = 0.0;
  double p1 = 0.0;

  double orig_value = 0.0;
  double laplacian = 0.0;

  for (int i=0; i<nGaussians; i++)
    {
      p0 = sTypeCoeffs(i,0);
      p1 = sTypeCoeffs(i,1);
      term_value = p1*exp(-p0*r_sq);

      orig_value += term_value;
      laplacian += 2*p0*(2*p0*r_sq-3)*term_value; 
    }

  laplacian = laplacian/orig_value;

  return -0.5*laplacian - Z/r_orig;
}

double QMCElectronNucleusCusp::getOrigGradient(double r_orig)
{
  double r_sq = r_orig*r_orig;
  int nGaussians = sTypeCoeffs.dim1();
  
  double term_gradient = 0.0;
  double p0 = 0.0;
  double p1 = 0.0;
  
  double orig_gradient = 0.0;

  for (int i=0; i<nGaussians; i++)
    {
      p0 = sTypeCoeffs(i,0);
      p1 = sTypeCoeffs(i,1);
      term_gradient = -p0*2*r_orig*p1*exp(-p0*r_sq);

      orig_gradient += term_gradient;
    }
  return orig_gradient;
}

double QMCElectronNucleusCusp::evaluateOrigOrbital(double r_orig)
{
  double r_sq = r_orig*r_orig;
  int nGaussians = sTypeCoeffs.dim1();
  
  double term_value = 0.0;
  double p0 = 0.0;
  double p1 = 0.0;
  
  double orig_value = 0.0;

  for (int i=0; i<nGaussians; i++)
    {
      p0 = sTypeCoeffs(i,0);
      p1 = sTypeCoeffs(i,1);
      term_value = p1*exp(-p0*r_sq);

      orig_value += term_value;
    }
  return orig_value;
}

ostream& operator <<(ostream& strm, QMCElectronNucleusCusp &rhs)
{
  // This function is never used.  It was added to get rid of some compiler 
  // errors.
  return strm;
}
