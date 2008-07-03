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

#include "QMCPotential_Energy.h"
#include "Lebedev-Laikov.h"
#include "QMCBasisFunction.h"
#include "MathFunctions.h"
#include "QMCJastrowElectronElectron.h"
#include "QMCJastrowElectronNuclear.h"
#include "QMCThreeBodyJastrow.h"
#include "Stopwatch.h"

QMCPotential_Energy::QMCPotential_Energy()
{}

void QMCPotential_Energy::operator=( const QMCPotential_Energy & rhs )
{
  Energy_total = rhs.Energy_total ;
  Energy_ne    = rhs.Energy_ne;
  Energy_ee    = rhs.Energy_ee;

  Input        = rhs.Input;
  HartreeFock  = rhs.HartreeFock;
  WF           = rhs.WF;
  MOL          = rhs.MOL;
  wd           = 0;

  ElectronNucleusCuspA = rhs.ElectronNucleusCuspA;
  ElectronNucleusCuspB = rhs.ElectronNucleusCuspB;

  P_nn = rhs.P_nn;
  P_en = rhs.P_en;
  P_ee = rhs.P_ee;
}

void QMCPotential_Energy::initialize(QMCInput * Ip, QMCHartreeFock * HF)
{
  Input       = Ip;
  HartreeFock = HF;
  wd          = 0;
  WF          = &Input->WF;
  MOL         = &Input->Molecule;

  Energy_total.allocate(globalInput.flags.walkers_per_pass);
  Energy_ne.allocate(globalInput.flags.walkers_per_pass);
  Energy_ee.allocate(globalInput.flags.walkers_per_pass);

  calc_P_nn();

  if (globalInput.flags.replace_electron_nucleus_cusps == 1)
    {
      coeffs = WF->getCoeff(true);
      ElectronNucleusCuspA.initialize(Input, * coeffs);
      ElectronNucleusCuspA.fitReplacementOrbitals();
      coeffs = WF->getCoeff(false);
      ElectronNucleusCuspB.initialize(Input, * coeffs);
      ElectronNucleusCuspB.fitReplacementOrbitals();
    }

  int numN = MOL->Atom_Positions.dim1();
  angularGrids.allocate(numN);
  for(int nuc=0; nuc<numN; nuc++)
    if(MOL->usesPsuedo(nuc)) 
      {
	Array2D<double> grTemp = MOL->getGrid(nuc,1.0,false);
	globalInput.BF.angularGrid(grTemp,nuc,angularGrids(nuc));
      }
}

double QMCPotential_Energy::evaluatePsuedoPotential(Array2D<double> & R, int elec, int nuc)
{
  bool isAlpha;
  int elecIndex;
  double r = wd->riA(elec,nuc);
  int lmax = (MOL->gridLegendre(nuc)).dim1();

  Array1D<double> vnl(lmax);
  Array1D<double> integral(lmax);

  /*
    Most psuedopotentials have local and nonlocal (i.e. requiring a separate integration)
    contributions.
  */
  double Vlocal = MOL->evaluatePotential(MOL->Vlocal(nuc),r);
  Vlocal -= MOL->Zeff(nuc) / r;
  
  bool small = true;
  for(int l=0; l<lmax; l++)
    {
      vnl(l) = MOL->evaluatePotential((MOL->Vnonlocal(nuc))(l),r);
      if(fabs(vnl(l)) > globalInput.flags.psuedo_cutoff) small = false;
    }

  //If all the nonlocal components are small, then we can cut out early.
  if(small) return Vlocal;

  if(elec < WF->getNumberElectrons(true))
    {
      isAlpha = true;
      elecIndex = elec;
      Dc_inv = wd->Dc_invA;
    } else {
      isAlpha = false;
      elecIndex = elec - WF->getNumberElectrons(true);
      Dc_inv = wd->Dc_invB;
    }

  /*
    The grid is composed of points on a spherical shell, all
    the same distance from the nucleus as the electron itself.
  */
  Array2D<double> grid = MOL->getGrid(nuc,r,true);

  int grSize = grid.dim1();
  X.allocate(grSize,WF->getNumberBasisFunctions());
  D.allocate(grSize,WF->getNumberOrbitals(isAlpha));
  ciDet.allocate(grSize,WF->getNumberElectrons(isAlpha));
  integrand.allocate(grSize);
  X.allocate(grSize,WF->getNumberBasisFunctions());

  /*
    Take a look at the paper from 1991:
    Nonlocal psuedopotentials and diffusion Monte Carlo
    Mitas, Shirley, Ceperley 
    http://link.aip.org/link/?JCPSA6/95/3467/1

    In their paper, they talk about Gaussian and Chebyshev type grids. But I am going
    to use Lebedev type grids because I was able to find good, recent code online for them, and
    they appear to be the grid of choice for DFT.

    We need to evaluate the wavefunction at all the grid points. We don't need
    to calculate the local energy, so we skip all the extra procedure in QMCSlater. We
    will want to evaluate the nonlocal integral W Psi / Psi, where W is a multiplicative operator,
    as seen in formula 28. Our integrand contains:

    \frac{ \Psi(r1, r2, ... , ri', ... , rn) }{ \Psi(r1, r2, ... , ri, ... , rn) }
    Where ri' are the grid points, and ri is the initial position of the electron.

    First, evaluate the basis functions and orbitals for the grid points around this nucleus
  */

  
  //These two function calls should produce the same result.
  //globalInput.BF.evaluateBasisFunctions(grid,X); 
  globalInput.BF.basisFunctionsOnGrid(grid,nuc,angularGrids(nuc),X);
  X.gemm(*WF->getCoeff(isAlpha),D,true);

  if (globalInput.flags.replace_electron_nucleus_cusps == 1)
    if(isAlpha)
      ElectronNucleusCuspA.replaceCusps(grid,0,D.dim1()-1,D);
    else
      ElectronNucleusCuspB.replaceCusps(grid,0,D.dim1()-1,D);

  integrand    = 0.0;
  double ratio = 1.0;
  double psi   = 0.0;
  for(int ci=0; ci<WF->getNumberDeterminants(); ci++){
    double termPsi = globalInput.WF.CI_coeffs(ci) * wd->DcA(ci) * wd->DcB(ci);
    psi += termPsi;
    WF->getDataForCI(ci,isAlpha,D,ciDet);

    double lastRatio = 1.0;
    for(int gp=0; gp<grSize; gp++){
      ratio = Dc_inv(ci).inverseUpdateOneRow(elecIndex,ciDet,gp);
      integrand(gp) += ratio * lastRatio * termPsi;
      lastRatio *= ratio; // ratio * lastRatio = \Psi_prime / \Psi
    }
  }

  /*
    Now we move on to evaluating the Jastrow function at all the grid points.

    We don't need to worry about electron-nucleus terms since \Jen(ri') = \Jen(ri), and therefore
    the terms cancel in the integrand.
  */
  double denom;
  QMCJastrowElectronElectron jee;
  denom = jee.jastrowOnGrid(globalInput.JP,elec,R,grid,integrand);
  psi *= denom;
  QMCJastrowElectronNuclear jen;
  denom = jen.jastrowOnGrid(globalInput.JP,elec,R,grid,integrand);
  psi *= denom;
  if(globalInput.flags.use_three_body_jastrow == 1){
    QMCThreeBodyJastrow jeen;
    denom = jeen.jastrowOnGrid(globalInput.JP,elec,R,grid,integrand);
    psi *= denom;
  }

  /*
    This is the actual integral part. The Legendre function selects the angular component
    to integrate and the weights are a part of the Lebedev-Laikov grid. See the Appendix
    of the aforementioned paper for a description of how Gaussian quadratures work for
    integration on the unit sphere, and see the src/Lebedev-Laikov.cpp code for references.
  */
  integral = 0.0;
  for(int l=0; l<lmax; l++)
    for(int gp=0; gp<grSize; gp++)
      integral(l) += integrand(gp) * (MOL->gridLegendre(nuc))(l,gp) * (MOL->gridWeights(nuc))(gp);

  /*
    This is the summation in formula 28
  */
  double Vnonlocal = 0.0;
  for(int l=0; l<lmax; l++)
    Vnonlocal += (2.0*l + 1.0) * vnl(l) * integral(l) / psi;

  return Vlocal + Vnonlocal;  
}

void QMCPotential_Energy::evaluate(Array1D<Array2D<double>*> &R,
				   Array1D<QMCWalkerData *> &walkerData,
				   int num)
{
  for(int walker = 0; walker < num; walker++)
    {
      wd = walkerData(walker);

      calc_P_en(*R(walker));
      calc_P_ee(*R(walker));
      Energy_total(walker) = P_nn + P_ee + P_en;

      Energy_ne(walker) = P_en;
      Energy_ee(walker) = P_ee;
    }

  wd = 0;
}

void QMCPotential_Energy::calc_P_en(Array2D<double> &R)
{
  P_en = 0.0;
  double Zeff;

  for(int nuc=0; nuc<MOL->Atom_Positions.dim1(); nuc++)
    {
      if(MOL->usesPsuedo(nuc))
	{
	  for(int el=0; el<R.dim1(); el++)
	    P_en += evaluatePsuedoPotential(R,el,nuc);
	} else {
	  Zeff = MOL->Zeff(nuc);	  
	  for(int el=0; el<R.dim1(); el++)
	    P_en -= Zeff/wd->riA(el,nuc);
        }
    }
}

void QMCPotential_Energy::calc_P_nn()
{
  P_nn = 0.0;
  double r;
  double chargei, chargej;

  //loop over every atom

  for(int i=0;i<MOL->Atom_Positions.dim1();i++)
    {
      chargei = MOL->Zeff(i);

      //loop over every other atom starting one beyond the current atom.
      //No double counting!

      for(int j=i+1;j<MOL->Atom_Positions.dim1();j++)
        {
          chargej = MOL->Zeff(j);
          r = MathFunctions::rij(MOL->Atom_Positions,i,j);
          P_nn += (chargei*chargej)/r;
        }
    }
}

void QMCPotential_Energy::calc_P_ee(Array2D<double> &R)
{
  P_ee = 0.0;

  double chargei, chargej;

  if (globalInput.flags.use_hf_potential == 1)
    {
      for (int i=0; i<R.dim1(); i++)
	P_ee += HartreeFock->GetVEff(i, R(i, 0), R(i, 1), R(i, 2));
    }
  else
    {
      // Loop over each unique pair of electrons.
      for (int i=0; i<R.dim1(); i++)
	{
	  chargei = -1.0;
	  for (int j=0; j<i; j++)
	    {
	      chargej = -1.0;
	      P_ee += (chargei*chargej)/wd->rij(i,j);
	    }
	}
    }
}

double QMCPotential_Energy::getEnergy(int which)
{
  return Energy_total(which);
}

double QMCPotential_Energy::getEnergyNE(int which)
{
  return Energy_ne(which);
}

double QMCPotential_Energy::getEnergyEE(int which)
{
  return Energy_ee(which);
}

