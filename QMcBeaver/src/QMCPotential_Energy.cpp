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

QMCPotential_Energy::QMCPotential_Energy()
{
}

void QMCPotential_Energy::operator=( const QMCPotential_Energy & rhs )
{
  Energy_total = rhs.Energy_total ;
  Input  = rhs.Input;

  P_nn = rhs.P_nn;
  P_en = rhs.P_en;
  P_ee = rhs.P_ee;
}

void QMCPotential_Energy::initialize(QMCInput *Ip)
{
  Input = Ip;
  calc_P_nn();
}

void QMCPotential_Energy::evaluate(Array2D<double> &R)
{
  calc_P_en(R);
  calc_P_ee(R);
  Energy_total = P_nn + P_ee + P_en;
}

void QMCPotential_Energy::calc_P_nn()
{
  P_nn = 0.0;
  double r;
  double chargei, chargej;

  //loop over every atom

  for(int i=0;i<Input->Molecule.Atom_Positions.dim1();i++)
    {
      chargei = Input->Molecule.Z(i);

      //loop over every other atom starting one beyond the current atom. 
      //No double counting!

      for(int j=i+1;j<Input->Molecule.Atom_Positions.dim1();j++)
	{
	  chargej = Input->Molecule.Z(j);
	  r = rij(Input->Molecule.Atom_Positions,i,j);
	  P_nn += (chargei*chargej)/r;
	}
    }
}

void QMCPotential_Energy::calc_P_en(Array2D<double> &R)
{
  P_en = 0.0;
  double r;
  double chargei, chargej;

  //loop over every atom

  for(int i=0;i<Input->Molecule.Atom_Positions.dim1();i++)
    {
      chargei = Input->Molecule.Z(i);

      //loop over every electron

      for(int j=0;j<R.dim1();j++)
	{
	  chargej = -1.0;

	  r = rij(Input->Molecule.Atom_Positions,R,i,j);
	  P_en += (chargei*chargej)/r;
	}
    }
}

void QMCPotential_Energy::calc_P_ee(Array2D<double> &R)
{
  P_ee = 0.0;
  double r;
  double chargei, chargej;
  
  //loop over every electron

  for(int i=0;i<R.dim1();i++)
    {
      chargei = -1.0;

      // loop over every electron that the ee energy hasn't be calculaed 
      // yet for

      for(int j=i+1;j<R.dim1();j++)
	{
	  chargej = -1.0;

	  r = rij(R,i,j);
	  P_ee += (chargei*chargej)/r;
	}
    }
}


double QMCPotential_Energy::rij(Array2D<double> &position, int i, int j)
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

double QMCPotential_Energy::rij(Array2D<double> &positioni,
				Array2D<double> &positionj, int i, int j)
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

double QMCPotential_Energy::getEnergy()
{
  return Energy_total;
}



