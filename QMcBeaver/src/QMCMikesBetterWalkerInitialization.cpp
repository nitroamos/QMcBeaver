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

#include "QMCMikesBetterWalkerInitialization.h"

#define PI 3.14159265359
#define a0 0.529177257507

QMCMikesBetterWalkerInitialization::
  QMCMikesBetterWalkerInitialization( QMCInput *INPUT )
{
  Input = INPUT;
}

Array2D<double> QMCMikesBetterWalkerInitialization::initializeWalkerPosition()
{
  int nelectrons = Input->WF.getNumberElectrons();
  Array2D<double> R(nelectrons,3);

  //Make a bunch of walkers
  int nTestWalkers = Input->flags.walker_initialization_combinations; 
  Array3D<double> severalRs(nelectrons,3,nTestWalkers);
  severalRs = initializeBunchOfWalkersPosition();

  //Pull the best walker/(combination of walkers) from this bunch
  R = FindBestWalker(severalRs);

  //Return this good walker
  return R;
}


Array3D<double> QMCMikesBetterWalkerInitialization::initializeBunchOfWalkersPosition()
{
  int nelectrons = Input->WF.getNumberElectrons();
  int nTestWalkers = Input->flags.walker_initialization_combinations;
  int natoms = Input->flags.Natoms;
  Array2D<double> atom_centers(natoms,3);
  atom_centers = Input->Molecule.Atom_Positions;
  Array3D<double> severalRs(nelectrons,3,nTestWalkers);
  double x0 = 0.0;
  double y0 = 0.0;
  double z0 = 0.0;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  int which_atom = 0;

  //For each electron
  for(int i=0;i<nelectrons;i++)
    {
      //For each TestWalker
      for(int j=0;j<nTestWalkers;j++)
	{
	  //invert most probable atomic orbital
	  //NOT DONE HERE!!!!!!!!
	  which_atom=0;//?????????
	  
	  x0=atom_centers(which_atom,0);
	  y0=atom_centers(which_atom,1);
	  z0=atom_centers(which_atom,2);
	  
	  //AOI.initialize(Xn,Yn,Zn,B,C);
	  //AOI.get_xyz(&Input->flags.iseed,x,y,z);
	  
	  x=x+x0;
	  y=y+y0;
	  z=z+z0;
	  
	  //END OF NOT DONE HERE!!!!!!!!

	  //put this electron into severalRs
	  severalRs(i,0,j)=x;
	  severalRs(i,1,j)=y;
	  severalRs(i,2,j)=z;
	}
    }

  return severalRs;
}

Array2D<double> QMCMikesBetterWalkerInitialization::
FindBestWalker(Array3D<double> &severalRs)
{
  //initialize data structures
  int nelectrons = severalRs.dim1();
  int nTestWalkers = severalRs.dim3(); 

  Array2D<double> Occupations(nelectrons,nTestWalkers);
  Array2D<double> GradOccupations(nelectrons,nTestWalkers);

  for(int i=0;i<nelectrons;i++)
    {
      for(int j=0;j<nTestWalkers;i++)
	{
	  Occupations(i,j) = 1.0/nTestWalkers;
	}
    }
  
  //We must recognize this is a constrained quadratic optimization problem
  //All the occupations are free excluding the last which is determined 
  //With the constraint that the sum of the occupations be one.

  double grad_bound=1.0;
  double step_size=0.1;
  int number_of_steps=4;

  //For number_of_steps we will move a maximum of step_size units in occupation
  for(int i=0;i<number_of_steps;i++)
    {
      //Find the gradient of the Objective function
      GradObjectiveFunctionForWalkers(severalRs,Occupations,GradOccupations);

      //Truncate gradient to bound max occupation change
      BoundGradOccupations(GradOccupations,grad_bound);

      //Take a step of size step_size units in the gradient direction
      MoveOccupations(Occupations,GradOccupations,step_size);
    }


  //Now will take the maximum Occupation for each electron and return
  //this electron to the user with the vector R[][]
  Array2D<double> R(nelectrons,3);

  //For each electron
  for(int i=0;i<nelectrons;i++)
    {
      //Find our guess that the best
      double max_occ=0.0;
      int max_occ_index=-99999;

      //For each TestWalker
      for(int j=0;j<nTestWalkers;j++)
	{
	  if(max_occ<Occupations(i,j))
	    {
	      max_occ=Occupations(i,j);
	      max_occ_index=j;
	    }
	}
      
      //Put the best into R
      R(i,0)=severalRs(i,0,max_occ_index);
      R(i,1)=severalRs(i,1,max_occ_index);
      R(i,2)=severalRs(i,2,max_occ_index);
    }
  return R;
}

void QMCMikesBetterWalkerInitialization::
FixConstraints(Array2D<double> &Occupations)
{
  double sum;
  for(int i=0;i<Occupations.dim1();i++)
    {
      sum=0.0;
      for(int j=0;j<(Occupations.dim2()-1);j++)
	{
	  //keep each electron in the right domain
	  if(Occupations(i,j)>1.0)
	    {
	      Occupations(i,j)=1.0;
	    }
	  //keep each electron in the right domain
	  if(Occupations(i,j)<0.0)
	    {
	      Occupations(i,j)=0.0;
	    }
	  sum = sum + Occupations(i,j);
	}
      if(sum>1.0)
	{
	  //pull the sum back down to 1.0
	  for(int j=0;j<(Occupations.dim2()-1);j++)
	    {
	      Occupations(i,j)=Occupations(i,j)/sum;
	    }
	  sum=1.0;
	}
      //fill in the last test electron's occupation to match the constraint
      Occupations(i,Occupations.dim2()-1) = 1.0-sum;
    }
}

double QMCMikesBetterWalkerInitialization::
ObjectiveFunctionForWalkers(Array3D<double> &severalRs,Array2D<double> &Occupations)
{
  double ee_energy = 0.0;
  double en_energy = 0.0;
  double temp=0.0;
  double x1,y1,z1,x2,y2,z2;
  double occ1,occ2;
  double r;

  //int nbeta  = Input->WF.getNumberBetaElectrons();  //number of beta electrons
  int nalpha = Input->WF.getNumberAlphaElectrons(); //number of alpha electrons

  //Get the electron-electron part of the energy
  //for each electron
  for(int i=0;i<severalRs.dim1();i++)
    {
      //for each combination of this electron
      for(int j=0;j<severalRs.dim3();j++)
	{
	  //for each other electron
	  for(int m=(i+1);m<severalRs.dim1();m++)
	    {
	      //for each combination of this electron
	      for(int n=0;n<severalRs.dim3();n++)
		{
		  x1=severalRs(i,0,j);
		  y1=severalRs(i,1,j);
		  z1=severalRs(i,2,j);
		  x2=severalRs(m,0,n);
		  y2=severalRs(m,1,n);
		  z2=severalRs(m,2,n);
		  occ1=Occupations(i,j);
		  occ2=Occupations(m,n);
		  
		  r=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));

		  //DETERMINE IF THEY ARE PARALLEL OR OPPOSITE SPIN ELECTRONS
		  //the alphas come first followed by the beta (right?)
		  //NOT DONE CHECK THIS WORK
		  if((i<nalpha && m<nalpha) || (i>=nalpha && m>=nalpha))
		    {
		      temp=Energy_parallel(r);
		    }
		  else
		    {
		      temp=Energy_opposite(r);
		    }
		  //END OF NOT DONE
		  ee_energy=ee_energy+temp*occ1*occ2;
		}
	    }
	}
    }

  //Get the electron-nuclear part of the energy
  int natoms = Input->flags.Natoms;
  Array2D<double> atom_centers(natoms,3);
  Array1D<int> atom_charges(natoms);
  atom_centers = Input->Molecule.Atom_Positions;
  atom_charges = Input->Molecule.Z;
  int charge;

  //for each electron
  for(int i=0;i<severalRs.dim1();i++)
    {
      //for each combination of this electron
      for(int j=0;j<severalRs.dim3();j++)
	{
	  //for each nuclei
	  for(int k=0;k<natoms;k++)
	    {
	      x1=severalRs(i,0,j);
	      y1=severalRs(i,1,j);
	      z1=severalRs(i,2,j);
	      x2=atom_centers(k,0);
	      y2=atom_centers(k,1);
	      z2=atom_centers(k,2);
	      occ1=Occupations(i,j);
	      
	      r=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
	      
	      charge=atom_charges(k);

	      temp=Energy_el_nuclr(r,charge);
	      
	      en_energy=en_energy+temp*occ1;
	    }
	}
    }

  return (ee_energy+en_energy);
}

void QMCMikesBetterWalkerInitialization::
GradObjectiveFunctionForWalkers(Array3D<double> &severalRs,
				Array2D<double> &Occupations,
				Array2D<double> &GradOccupations)
{
  double initial_energy = ObjectiveFunctionForWalkers(severalRs,Occupations);
  double dr = 0.1;
  double perturb_energy;
  double r_orig;

  //for each electron
  for(int i=0;i<severalRs.dim1();i++)
    {
      //for each combination of this electron
      for(int j=0;j<severalRs.dim3();j++)
	{
	  r_orig = Occupations(i,j);

	  Occupations(i,j) = Occupations(i,j)+dr;
	  FixConstraints(Occupations);

	  perturb_energy = ObjectiveFunctionForWalkers(severalRs,Occupations);
      
	  GradOccupations(i,j) = (perturb_energy-initial_energy)/dr;

	  Occupations(i,j) = r_orig;
	  FixConstraints(Occupations);
	}
    }
}

void QMCMikesBetterWalkerInitialization::
BoundGradOccupations(Array2D<double> &GradOccupations,double bound)
{
  double mag,temp;
  //for each electron
  for(int i=0;i<GradOccupations.dim1();i++)
    {
      mag = 0.0;
      //for each combination of this electron
      for(int j=0;j<GradOccupations.dim2();j++)
	{
	  temp = GradOccupations(i,j);
	  mag = mag+temp*temp;
	}

      mag = sqrt(mag);

      for(int j=0;j<GradOccupations.dim2();j++)
	{
	  GradOccupations(i,j) = GradOccupations(i,j)/mag*bound;
	}
    }
}

void QMCMikesBetterWalkerInitialization::
MoveOccupations(Array2D<double> &Occupations,
		Array2D<double> &GradOccupations,
		double dr)
{
  //for each electron
  for(int i=0;i<GradOccupations.dim1();i++)
    {
      //for each combination of this electron
      for(int j=0;j<(GradOccupations.dim2()-1);j++)
	{
	  Occupations(i,j)=Occupations(i,j)
	    +dr*GradOccupations(i,j);
	}
      FixConstraints(Occupations);
    }
}

double QMCMikesBetterWalkerInitialization::Energy_parallel(double r)
{
  //can make this slightly worse since we want the parallels to stay apart more
  return 2.0/r;
}

double QMCMikesBetterWalkerInitialization::Energy_opposite(double r)
{
  return 1.0/r;
}

double QMCMikesBetterWalkerInitialization::Energy_el_nuclr(double r, 
							   int charge)
{
  return -1.0*charge/r;
}










