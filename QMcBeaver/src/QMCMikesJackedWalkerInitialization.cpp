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

#include "QMCMikesJackedWalkerInitialization.h"

#define PI 3.14159265359
#define a0 0.529177257507

QMCMikesJackedWalkerInitialization::
  QMCMikesJackedWalkerInitialization( QMCInput *INPUT )
{
  Input = INPUT;
}

Array2D<double> QMCMikesJackedWalkerInitialization::initializeWalkerPosition()
{
  int natoms = Input->flags.Natoms;
  Array2D<double> atom_centers(natoms,3);
  atom_centers = Input->Molecule.Atom_Positions;

  int nelectrons = Input->WF.getNumberElectrons();
  int nelectrons_left = nelectrons;
  int nbeta  = Input->WF.getNumberBetaElectrons();  //number of beta electrons
  int nalpha = Input->WF.getNumberAlphaElectrons(); //number of alpha electrons
  int nbeta_left  = nbeta;
  int nalpha_left = nalpha;

  Array2D <double> atom_e_count(natoms,2);

  atom_e_count = electrons_and_radii(); 

  //returns the number of electrons on each atom and the relative 
  //size of the atom

  Array2D <double> atom_ab_count(natoms,3);  //this will hold the number 
                                             //of total,alpha, and beta 
                                             //electrons on each atom

  int n_e_pairs = nalpha;
  if(nalpha>nbeta)
    {
      n_e_pairs = nbeta;
    }
  int n_e_pairs_left = n_e_pairs;

  int n_lone_e = nelectrons-2*n_e_pairs;
  int n_lone_e_left = n_lone_e;

  for(int np=0; np<natoms;np++)
    {
      atom_ab_count(np,0) = atom_e_count(np,0);  
      //total number of electrons on this atom
      
      atom_ab_count(np,1) = 0;  //zero alphas
      atom_ab_count(np,2) = 0;  //zero betas
    }

  //fill in unpaired electrons on atoms to make the number of 
  //electrons even on each atom

  int flag1, flag2;
  for(int np=0;np<natoms;np++)
    {  
      flag1 = 0;
      flag2 = 0;
      if( ( (int) (0.1+atom_ab_count(np,0)) ) %2 == 1)  //odd number
	{
	  //find out which electron type we have more of
	  if(nalpha_left>nbeta_left)      flag1 = 1;
	  else if(nalpha_left<nbeta_left) flag1 = -1;
	  else                            flag2 = 1;    //they are equal

	  //what should we do with this situation
	  if(flag1==1)  //add an alpha
	    {
	      nalpha_left--;
	      nelectrons_left--;
	      n_lone_e_left--;
	      atom_ab_count(np,0)--;
	      atom_ab_count(np,1)++;
	    }
	  else if(flag1==-1)  //add a beta
	    {
	      nbeta_left--;
	      nelectrons_left--;
	      n_lone_e_left--;
	      atom_ab_count(np,0)--;
	      atom_ab_count(np,2)++; 
	    }
	  else if(flag2==1)   //there are no lone electrons yet 
	                      //this atom has an odd number of electrons
	    {
	      //add an alpha
	      nalpha_left--;
	      nelectrons_left--;
	      n_lone_e_left++;
	      n_e_pairs_left--;
	      atom_ab_count(np,0)--;
	      atom_ab_count(np,1)++;
	    }
	  else  //can't happen
	    {
	      cerr<< "ERROR: Incorrect number of electrons left A" << endl;
	      exit(1);
	    }
	}
    }
  //now there are only even numbers to put on each atom

  //we will now fill up each atom with electrons two at a time 
  //until we run out of pairs

  while(n_e_pairs_left!=0)
    {
      for(int np=0;np<natoms;np++)
	{
	  if(atom_ab_count(np,0) >= 2)
	    //we still can put at least 2 electrons on this atom
	    {
	      atom_ab_count(np,0) -= 2; 
	      //putting 2 electrons on this atom
	      atom_ab_count(np,1)++;  //adding an alpha
	      atom_ab_count(np,2)++;  //adding a beta
	      nelectrons_left -= 2;
	      n_e_pairs_left--;
	      nalpha_left--;
	      nbeta_left--;
	    }
	}
    }

  //now we have only lone electrons to fill up the remaining atoms with
  while(n_lone_e_left!=0)
    {
      for(int np=0;np<natoms;np++)
	{
	  if(atom_ab_count(np,0)>0)
	    //we can still put at least one electron on this atom
	    {
	      if(nalpha_left>0)  //first we will use all the alphas
		{
		  nalpha_left--;
		  nelectrons_left--;
		  n_lone_e_left--;
		  atom_ab_count(np,0)--;
		  atom_ab_count(np,1)++; 
		}
	      else if(nbeta_left>0)  //then we will use all the betas
		{
		  nbeta_left--;
		  nelectrons_left--;
		  n_lone_e_left--;
		  atom_ab_count(np,0)--;
		  atom_ab_count(np,2)++;  
		}
	      else
		{
		  cerr << "Error: Alpha and beta counting has gone bad!"
		       << endl;
		  exit(1);
		}
	    }
	}
    }

  //we now have no n_e_pairs left and no n_lone_pairs_left
  //check our work

  if(nalpha_left!=0)
    {
      cerr<<"Error: In nalpha_left!=0" << endl;
      exit(1);
    }
  if(nbeta_left!=0)
    {
      cerr << "Error: In nbeta_left!=0" << endl;
      exit(1);
    }
  if(n_e_pairs_left!=0)
    {
      cerr << "Error: In n_e_pairs_left!=0" << endl;
      exit(1);
    }
  if(n_lone_e_left!=0)
    {
      cerr << "Error: In n_lone_e_left!=0" << endl;
      exit(1);
    }
  if(nelectrons_left!=0)
    {
      cerr << "Error: In nelectrons_left!=0" << endl;
      exit(1);
    }

  //now we will actually start filling in the R with respect to
  //the number of alpha and beta that we have on each atom

  double x,y,z,xx,yy,zz,SD,atom_size,AA;

  //these will keep track of where to put the next alpha electron
  int alpha_index = 0;
  int beta_index  = nalpha;

  Array2D<double> R(nelectrons,3);
  for(int np=0;np<natoms;np++)
    { 
      x = atom_centers(np,0);
      y = atom_centers(np,1);
      z = atom_centers(np,2);
      nalpha_left = (int)(0.1+atom_ab_count(np,1));
      nbeta_left = (int)(0.1+atom_ab_count(np,2));
      atom_size = atom_e_count(np,1);
      while(nalpha_left!=0)
	{
	  SD = atom_size;  //gaussian with SD as the 3D standard deviation
	  AA = SD*sqrt(2*PI/(3*PI-8.0));
	  //AA is what makes the 3D gaussian use SD as 
	  //the standard deviation

	  xx = AA*ran.gasdev();
	  yy = AA*ran.gasdev();
	  zz = AA*ran.gasdev();

	  R(alpha_index,0) = xx+x;
	  R(alpha_index,1) = yy+y;
	  R(alpha_index,2) = zz+z;
	  alpha_index++;
	  nalpha_left--;
	}
      while(nbeta_left!=0)
	{
	  SD = atom_size; //gaussian with SD as the 3D standard deviation
	  AA = SD*sqrt(2*PI/(3*PI-8.0));
	  //AA is what makes the 3D gaussian use SD as 
	  //the standard deviation

	  xx = AA*ran.gasdev();
	  yy = AA*ran.gasdev();
	  zz = AA*ran.gasdev();

	  R(beta_index,0) = xx+x;
	  R(beta_index,1) = yy+y;
	  R(beta_index,2) = zz+z;
	  beta_index++;
	  nbeta_left--;
	}  
    } 

  return R;
}

Array2D <double> QMCMikesJackedWalkerInitialization::electrons_and_radii()
{
    int Norbitals  = Input->WF.getNumberOrbitals();
    int Nbasisfunc = Input->WF.getNumberBasisFunctions();

    Array2D<qmcfloat> OrbitalScratch(Input->WF.OrbitalCoeffs);
    
    //Square MO coefficients
    for(int i=0; i<Norbitals; i++)
      for(int j=0; j<Nbasisfunc; j++)
        {
	  OrbitalScratch(i,j) = OrbitalScratch(i,j)*OrbitalScratch(i,j);
        }
    
    //"Normalize" MO by summing the squared coeffs and dividing

    double OrbitalSum;

    for(int i=0; i<Norbitals; i++)
      {
        OrbitalSum = 0;
        for(int j=0; j<Nbasisfunc; j++) 
	  {
	    OrbitalSum += OrbitalScratch(i,j);
	  }
        
        for(int j=0; j<Nbasisfunc; j++)
	  {
            OrbitalScratch(i,j) = OrbitalScratch(i,j)/OrbitalSum;
	  }
      }
    
    // Create an Array that links the index of an orbital to the atom it 
    // is associated with
    Array1D <int> OrbPos(Input->WF.getNumberBasisFunctions());

    int Natoms = Input->flags.Natoms;
    
    int i=0;
    for(int j=0; j<Natoms; j++)
      for(int k=0; k<Input->BF.getNumberBasisFunctions(j); k++)
	{
	  OrbPos(i) = j; 
	  i++;
	}
    
    //Initialize the array to be returned.  The first column has the
    //number of electrons for the atom and the second column has the
    //radi of the atom
    
    Array2D <double> EandR(Natoms,2);
    for(int i=0; i<Natoms; i++) 
      for(int j=0; j<2; j++) 
	EandR(i,j) = 0.0;
    
    double rv;  //uniform [0,1) random variable
    double sum;
    for(int i=0; i<Norbitals; i++)
      {
	if (Input->WF.AlphaOccupation(0,i) != Input->WF.getUnusedIndicator())
	  {
	    sum = 0.0;
            rv = ran.unidev();
            
            for(int j=0; j<Nbasisfunc; j++)
	      {
                sum += OrbitalScratch(i,j);
                if (sum >= rv) 
		  { 
		    EandR(OrbPos(j),0) += 1; 
		    break; 
		  }
	      }
	  }
	
	if (Input->WF.BetaOccupation(0,i) != Input->WF.getUnusedIndicator())
	  {
	    sum = 0.0;
	    rv = ran.unidev();

	    for (int k=0; k<Nbasisfunc; k++)
	      {
		sum += OrbitalScratch(i,k);
		if (sum >= rv)
		  {
		    EandR(OrbPos(k),0) += 1;
		    break;
		  }
	      }
	  }
      }

    for(int i=0; i<Natoms; i++)
      {
        EandR(i,1) = covalent_radi(Input->Molecule.Z(i));
      }
    
    return EandR;
}


double QMCMikesJackedWalkerInitialization::covalent_radi(int ZZ)
{

  double value = 0;

  if(ZZ < 1) 
    {
      cerr << "ERROR: Z < 1" << endl;
      cerr << "Z = " << ZZ << endl;
      exit(1);
    }

  //Table of covalent radi in angstroms
  switch(ZZ)
    {
    case 1:  value = 0.15; break;//
    case 2:  value = 0.25; break;//
    case 3:  value = 0.40; break;//
    case 4:  value = 0.30; break;//
    case 5:  value = 0.30; break;//
    case 6:  value = 0.30; break;//
    case 7:  value = 0.50; break;//
    case 8:  value = 0.73; break;
    case 9:  value = 0.72; break;
    case 10:  value = 0.71; break;
    case 11:  value = 1.54; break;
    case 12:  value = 1.36; break;
    case 13:  value = 1.18; break;
    case 14:  value = 1.11; break;
    case 15:  value = 1.06; break;
    case 16:  value = 1.02; break;
    case 17:  value = 0.99; break;
    case 18:  value = 0.98; break;
    case 19:  value = 2.03; break;
    case 20:  value = 1.74; break;
    case 21:  value = 1.44; break;
    case 22:  value = 1.32; break;
    case 23:  value = 1.22; break;
    case 24:  value = 1.18; break;
    case 25:  value = 1.17; break;
    case 26:  value = 1.17; break;
    case 27:  value = 1.16; break;
    case 28:  value = 1.15; break;
    case 29:  value = 1.17; break;
    case 30:  value = 1.25; break;
    case 31:  value = 1.26; break;
    case 32:  value = 1.22; break;
    case 33:  value = 1.20; break;
    case 34:  value = 1.16; break;
    case 35:  value = 1.14; break;
    case 36:  value = 1.89; break;
    case 37:  value = 2.16; break;
    case 38:  value = 1.91; break;
    case 39:  value = 1.62; break;
    case 40:  value = 1.45; break;
    case 41:  value = 1.34; break;
    case 42:  value = 1.30; break;
    case 43:  value = 1.27; break;
    case 44:  value = 1.25; break;
    case 45:  value = 1.25; break;
    case 46:  value = 1.28; break;
    case 47:  value = 1.34; break;
    case 48:  value = 1.41; break;
    case 49:  value = 1.44; break;
    case 50:  value = 1.41; break;
    case 51:  value = 1.40; break;
    case 52:  value = 1.36; break;
    case 53:  value = 1.33; break;
    case 54:  value = 1.31; break;
    case 55:  value = 2.35; break;
    case 56:  value = 1.98; break;
    case 57:  value = 1.25; break;
    case 58:  value = 1.65; break;
    case 59:  value = 1.65; break;
    case 60:  value = 1.64; break;
    case 61:  value = 1.63; break;
    case 62:  value = 1.62; break;
    case 63:  value = 1.85; break;
    case 64:  value = 1.61; break;
    case 65:  value = 1.59; break;
    case 66:  value = 1.59; break;
    case 67:  value = 1.58; break;
    case 68:  value = 1.57; break;
    case 69:  value = 1.56; break;
    case 70:  value = 1.70; break;
    case 71:  value = 1.56; break;
    case 72:  value = 1.44; break;
    case 73:  value = 1.34; break;
    case 74:  value = 1.30; break;
    case 75:  value = 1.28; break;
    case 76:  value = 1.26; break;
    case 77:  value = 1.27; break;
    case 78:  value = 1.30; break;
    case 79:  value = 1.34; break;
    case 80:  value = 1.49; break;
    case 81:  value = 1.48; break;
    case 82:  value = 1.47; break;
    case 83:  value = 1.46; break;
    case 84:  value = 1.53; break;
    case 85:  value = 1.47; break;
    case 90:  value = 1.65; break;
    case 92:  value = 1.42; break;
    }

  if(value == 0)
    {
      cerr << "ERROR: Do not have a covalent radi for Z = " << ZZ << endl;
      exit(1);
    }

  return value/a0;
}
















