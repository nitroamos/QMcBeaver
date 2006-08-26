// This is dans_walker_initialization.  It creates initial configurations that 
// require very little time for the run to equilibrate.

#include "QMCDansWalkerInitialization.h"

#define PI 3.14159265359

QMCDansWalkerInitialization::QMCDansWalkerInitialization(QMCInput * INPUT)
{
  Input = INPUT;
  if (!arraysInitialized)
    {
      arraysInitialized = true;
      initializeArrays();
    }
}

bool QMCDansWalkerInitialization::arraysInitialized = false;

Array1D<CubicSpline> QMCDansWalkerInitialization::phiSplines;
Array1D<int> QMCDansWalkerInitialization::phiSplinesMade;

Array1D<CubicSpline> QMCDansWalkerInitialization::thetaSplines;
Array1D<int> QMCDansWalkerInitialization::thetaSplinesMade;

Array2D<CubicSpline> QMCDansWalkerInitialization::radialSplines;
Array2D<int> QMCDansWalkerInitialization::radialSplinesMade;

Array1D<double> QMCDansWalkerInitialization::x_array;

void QMCDansWalkerInitialization::initializeArrays()
{
  phiSplines.allocate(19);
  phiSplinesMade.allocate(19);
  phiSplinesMade = 0;

  thetaSplines.allocate(19);
  thetaSplinesMade.allocate(19);
  thetaSplinesMade = 0;

  radialSplines.allocate(18,3);
  radialSplinesMade.allocate(18,3);
  radialSplinesMade = 0;

  x_array.allocate(21);
  for (int i=0; i<21; i++) x_array(i) = i/20.0;
}

Array2D<double> QMCDansWalkerInitialization::initializeWalkerPosition()
{
  int natoms = Input->flags.Natoms;
  Array2D<double> atom_centers(natoms,3);
  atom_centers = Input->Molecule.Atom_Positions;

  int nelectrons = Input->WF.getNumberElectrons();
  int nbeta = Input->WF.getNumberBetaElectrons();
  int nalpha = Input->WF.getNumberAlphaElectrons();

  // This array holds the numbers of total, alpha, and beta electrons on each
  // center.

  Array2D<int> ab_count(natoms,3);
  ab_count = assign_electrons_to_nuclei();

  // Redistribute electrons if there are charged centers.

  int icharge_init, icharge, kcharge, partner, partnercharge, give, get;

  for (int i=0; i<natoms; i++)
    if (Input->Molecule.Z(i) != ab_count(i,0)) // If this atom is charged
      {
	icharge_init = Input->Molecule.Z(i) - ab_count(i,0);
        for (int j=0; j<abs(icharge_init); j++)
	  {
	    icharge = Input->Molecule.Z(i) - ab_count(i,0);
	    partner = -1;
	    give = -1;
	    get = -1;
	    partnercharge = 0;
	    for (int k=0; k<natoms; k++)  // Find an atom to exchange with
	      {
		if (k != i)
		  {
		    kcharge = Input->Molecule.Z(k) - ab_count(k,0);
		    if (icharge < 0 && kcharge >= icharge+2)
		      {
			if (partner == -1)
			  {
			    partner = k;
			    partnercharge = kcharge;
			  }
			else if (kcharge > partnercharge)
			  {
			    partner = k;
			    partnercharge = kcharge;
			  }
		      }
		    else if (icharge > 0 && kcharge <= icharge-2)
		      {
			if (partner == -1)
			  {
			    partner = k;
			    partnercharge = kcharge;
			  }
			else if (kcharge < partnercharge)
			  {
			    partner = k;
			    partnercharge = kcharge;
			  }
		      }
		  }
	      }  // Now we know which atom we are trading electrons with
	    if (partner != -1)
	      {
		if (icharge < partnercharge)
		  // We are giving an electron to the partner
		  {
		    give = i;
		    get = partner;
		  }
		else if (icharge > partnercharge)
		  // We are getting an electron from the partner
		  {
		    get = i;
		    give = partner;
		  }
		if (ab_count(give,1) == ab_count(give,2))  // if na=nb on -
		  {
		    if (ab_count(get,1) <= ab_count(get,2))
		      {  // transfer an alpha
			ab_count(give,0) -= 1;
			ab_count(give,1) -= 1;
			ab_count(get,0) += 1;
			ab_count(get,1) += 1;
		      }
		    else if (ab_count(get,1) > ab_count(get,2))
		      {  // transfer a beta
			ab_count(give,0) -= 1;
			ab_count(give,2) -= 1;
			ab_count(get,0) += 1;
			ab_count(get,2) += 1;
		      }
		  }
		else if (ab_count(give,1) > ab_count(give,2)) // if na>nb on -
		  {
		    if (ab_count(get,1) <= ab_count(get,2))
		      {  // transfer an alpha
			ab_count(give,0) -= 1;
			ab_count(give,1) -= 1;
			ab_count(get,0) += 1;
			ab_count(get,1) += 1;
		      }
		    else if (ab_count(get,1) > ab_count(get,2))
		      // find a center where nb>na, and trade an alpha from - 
		      // for a beta for +.
		      {
			for (int l=0; l<natoms; l++) 
			  if ( (l != give && l != get) && 
			       (ab_count(l,1) < ab_count(l,2)) )
			    { 
			      ab_count(give,0) -= 1;
			      ab_count(give,1) -= 1;
			      ab_count(l,1) += 1;
			      ab_count(l,2) -= 1;
			      ab_count(get,0) += 1;
			      ab_count(get,2) += 1;
			      break;
			    }
		      }
		  }
		else if (ab_count(give,1) < ab_count(give,2))
		  {
		    if (ab_count(get,1) >= ab_count(get,2))
		      { // transfer a beta
			ab_count(give,0) -= 1;
			ab_count(give,2) -= 1;
			ab_count(get,0) += 1;
			ab_count(get,2) += 1;
		      }
		    else if (ab_count(get,1) < ab_count(get,2))
		      // find a center where na>nb, and trade a beta from -
		      // for an alpha for +
		      {
			for (int m=0; m<natoms; m++)
			  if ( (m != give && m != get) &&  
			       (ab_count(m,1) > ab_count(m,2)) )
			    {
			      ab_count(give,0) -= 1;
			      ab_count(give,2) -= 1;
			      ab_count(m,2) += 1;
			      ab_count(m,1) -= 1;
			      ab_count(get,0) += 1;
			      ab_count(get,1) += 1;
			      break;
			    }
		      }
		  }
	      }
	  }
      }

  // Now we check to make sure no center has too many electrons of the same 
  // type.

  for (int i=0; i<natoms; i++)
    {
      if (Input->Molecule.Z(i) <= 18 && ab_count(i,1) > 9)
	{
	  int extra_alphas = ab_count(i,1)-9;
	  for (int j=0; j<extra_alphas; j++)
	    for (int k=0; k<natoms; k++)
	      if (k != i && ab_count(k,1) < ab_count(k,2))
		{
		  ab_count(i,1) -= 1;
		  ab_count(k,1) += 1;
		  ab_count(i,2) += 1;
		  ab_count(k,2) -= 1;
		  break;
		}
	}
      if (Input->Molecule.Z(i) <= 18 && ab_count(i,2) > 9)
	{
	  int extra_betas = ab_count(i,2)-9;
	  for (int m=0; m<extra_betas; m++)
	    for (int n=0; n<natoms; n++)
	      if (n != i && ab_count(n,1) > ab_count(n,2))
		{
		  ab_count(i,2) -= 1;
		  ab_count(n,2) += 1;
		  ab_count(i,1) += 1;
		  ab_count(n,1) -= 1;
		  break;
		}
	}
      if (Input->Molecule.Z(i) <= 10 && ab_count(i,1) > 5)
	{
	  int extra_alphas = ab_count(i,1)-5;
	  for (int j=0; j<extra_alphas; j++)
	    for (int k=0; k<natoms; k++)
	      if (k != i && ab_count(k,1) < ab_count(k,2))
		{ 
		  ab_count(i,1) -= 1;
	          ab_count(k,1) += 1;
	          ab_count(i,2) += 1;
	          ab_count(k,2) -= 1;
	          break;
	        }
	}
      if (Input->Molecule.Z(i) <= 10 && ab_count(i,2) > 5)
	{
	  int extra_betas = ab_count(i,2)-5;
	  for (int m=0; m<extra_betas; m++)
	    for (int n=0; n<natoms; n++)
	      if (n != i && ab_count(n,1) > ab_count(n,2))
		{
		  ab_count(i,2) -= 1;
		  ab_count(n,2) += 1;
		  ab_count(i,1) += 1;
		  ab_count(n,1) -= 1;
		  break;
		}
	}
      if (Input->Molecule.Z(i) <= 2 && ab_count(i,1) > 1)
	{
	  int extra_alphas = ab_count(i,1)-1;
	  for (int p=0; p<extra_alphas; p++)
	    for (int q=0; q<natoms; q++)
              if (q != i && ab_count(q,1) < ab_count(q,2))
                {
                  ab_count(i,1) -= 1;
                  ab_count(q,1) += 1;
                  ab_count(i,2) += 1;
                  ab_count(q,2) -= 1;
                  break;
                }
	}
      if (Input->Molecule.Z(i) <= 2 && ab_count(i,2) > 1)
	{
	  int extra_betas = ab_count(i,2)-1;
	  for (int r=0; r<extra_betas; r++)
	    for (int s=0; s<natoms; s++)
              if (s != i && ab_count(s,1) > ab_count(s,2))
                {
                  ab_count(i,2) -= 1;
                  ab_count(i,1) += 1;
                  ab_count(s,2) += 1;
                  ab_count(s,1) -= 1;
                  break;
                }
	}
    }

  // Check to see that all electrons have been assigned.

  int check_elec = 0;
  int check_alpha = 0;
  int check_beta = 0;

  for (int i=0; i<natoms; i++)
    {
      if ( (ab_count(i,1)+ab_count(i,2)) != ab_count(i,0) )
	{
	  cerr << "Alpha and beta electrons do not add up for atom " << i;
	  cerr << endl;
	  exit(1);
	}
      check_elec += ab_count(i,0);
      check_alpha += ab_count(i,1);
      check_beta += ab_count(i,2);
    }

  if (check_elec != nelectrons)
    { 
      cerr << "Electron counting has gone bad." << endl; 
      exit(1); 
    }

  if (check_alpha != nalpha)
    { 
      cerr << "Alpha counting has gone bad." << endl; 
      exit(1); 
    }
  
  if (check_beta != nbeta)
    { 
      cerr << "Beta counting has gone bad." << endl; 
      exit(1); 
    }

  // We put the coordinates for the alpha electrons into one array, and the 
  // coordinates for the beta electrons in another one.

  Array2D<double> temp_alpha(nalpha,3);
  int alpha_index = 0;

  Array2D<double> temp_beta(nbeta,3);
  int beta_index = 0;

  Array2D<double> temp_coords;
  int n_e,n_a,n_b;

  for (int i=0; i<natoms; i++)
    {
      n_e = ab_count(i,0);
      n_a = ab_count(i,1);
      n_b = ab_count(i,2);

      if (n_e > 0)
	{
	  temp_coords.allocate(n_e,3);
	  temp_coords = dist_center(Input->Molecule.Z(i),n_e,n_a,n_b);
	
	  for (int j=0; j<n_a; j++)
	    {
	      for (int k=0; k<3; k++)
		temp_alpha(alpha_index,k)=temp_coords(j,k)+atom_centers(i,k);
	      alpha_index++;
	    }
      
	  for (int m=0; m<n_b; m++)
	    {
	      for (int n=0; n<3; n++)
		temp_beta(beta_index,n)=temp_coords(n_a+m,n)+atom_centers(i,n);
	      beta_index++;
	    }
	}
    }

  // The electronic coordinates are stored in this array, with the alpha
  // electrons first, then the betas.

  Array2D<double> R(nelectrons,3);

  for (int i=0; i<nalpha; i++)
    for (int j=0; j<3; j++)
      R(i,j) = temp_alpha(i,j);

  for (int i=0; i<nbeta; i++)
    for (int j=0; j<3; j++)
      R(i+nalpha,j) = temp_beta(i,j);

  // We shuffle the order of the alpha electrons and the order of the beta 
  // electrons.

  /*
  Array1D<double> rnd;
  Array1D<int> index;

  double temp_rnd;
  int temp_index;

  if (nalpha > 1)
    {
      rnd.allocate(nalpha);
      index.allocate(nalpha);

      for (int i=0; i<nalpha; i++)
	{
	  rnd(i) = ran1(&Input->flags.iseed);
	  index(i) = i;
	}

      for (int i=1; i<nalpha; i++)
	for (int j=0; j<nalpha-i; j++)
	  if (rnd(j) > rnd(j+1))
	    {
	      temp_rnd = rnd(j+1);
	      temp_index = index(j+1);
	      
	      rnd(j+1) = rnd(j);
	      rnd(j) = temp_rnd;
	      
	      index(j+1) = index(j);
	      index(j) = temp_index;
	    }

      for (int i=0; i<nalpha; i++)
	for (int j=0; j<3; j++)
	  R(i,j) = temp_alpha(index(i),j);
    }

  if (nbeta > 1)
    {
      rnd.allocate(nbeta);
      index.allocate(nbeta);

      for (int i=0; i<nbeta; i++)
	{
	  rnd(i) = ran1(&Input->flags.iseed);
	  index(i) = i;
	}
      
      for (int i=1; i<nbeta; i++)
	for (int j=0; j<nbeta-i; j++)
	  if (rnd(j) > rnd(j+1))
	    {
	      temp_rnd = rnd(j+1);
	      temp_index = index(j+1);
	      
	      rnd(j+1) = rnd(j);
	      rnd(j) = temp_rnd;

	      index(j+1) = index(j);
	      index(j) = temp_index;
	    }

      for (int i=0; i<nbeta; i++)
	for (int j=0; j<3; j++)
	  R(i+nalpha,j) = temp_beta(index(i),j);
    }
  */

  return R;
}

Array2D<int> QMCDansWalkerInitialization::assign_electrons_to_nuclei()
{
  int Norbitals  = Input->flags.Norbitals;
  int Nbasisfunc = Input->flags.Nbasisfunc;

  Array2D<qmcfloat> AlphaScratch(Input->WF.AlphaCoeffs);
  Array2D<qmcfloat> BetaScratch(Input->WF.BetaCoeffs);

  // Square MO coefficients.

  for (int i=0; i<Norbitals; i++)
    for (int j=0; j<Nbasisfunc; j++)
      {
	AlphaScratch(i,j) = AlphaScratch(i,j)*AlphaScratch(i,j);
	BetaScratch(i,j) = BetaScratch(i,j)*BetaScratch(i,j);
      }
  //"Normalize" MO by summing the squared coeffs and dividing

  double AlphaSum;
  double BetaSum;

  for (int i=0; i<Norbitals; i++)
    {
      AlphaSum = 0.0;
      BetaSum = 0.0;
      for(int j=0; j<Nbasisfunc; j++)
	{ 
	  AlphaSum += AlphaScratch(i,j);
	  BetaSum += BetaScratch(i,j);
	}
      for(int j=0; j<Nbasisfunc; j++)
	{
	  AlphaScratch(i,j) = AlphaScratch(i,j)/AlphaSum;
	  BetaScratch(i,j) = BetaScratch(i,j)/BetaSum;
	}
    }

  // Create an Array that links the index of an basis function to the atom it
  // is associated with

  int Natoms = Input->flags.Natoms;

  Array1D<int> OrbPos(Nbasisfunc);

  int i = 0;
  for (int j=0; j<Natoms; j++)
    for(int k=0; k<(Input->BF.getNumberBasisFunctions(j)); k++)
      {
	OrbPos(i) = j;
	i++;
      }

  //Initialize the array to be returned.

  Array2D <int> atom_occ(Natoms,3);
  atom_occ = 0;

  double sum = 0.0;
  double rv;  //uniform [0,1] random variable

  // First we decide which determinant we will use to distribute the elecs.
  rv = ran1(&Input->flags.iseed);

  int determinant_index = 0;
  for (int i=0; i<Input->WF.getNumberDeterminants(); i++)
    {
      sum += Input->WF.CI_coeffs(i)*Input->WF.CI_coeffs(i);
      if (sum >= rv)
	{
	  determinant_index = i;
	  break;
	}
    }

  for (int i=0; i<Norbitals; i++)
    {
      if (Input->WF.AlphaOccupation(determinant_index,i) == 1)
	{
          sum = 0.0;
          rv = ran1(&Input->flags.iseed);

          for (int j=0; j<Nbasisfunc; j++)
            {
              sum += AlphaScratch(i,j);
              if (sum >= rv) 
		{
		  atom_occ(OrbPos(j),0) += 1;
                  atom_occ(OrbPos(j),1) += 1; 
                  break;
		}
	    }
	}

      if (Input->WF.BetaOccupation(determinant_index,i) == 1)
	{
	  sum = 0.0;
          rv = ran1(&Input->flags.iseed);

	  for (int k=0; k<Nbasisfunc; k++)
	    {
	      sum += BetaScratch(i,k);
	      if (sum >= rv)
		{
		  atom_occ(OrbPos(k),0) += 1;
		  atom_occ(OrbPos(k),2) += 1; 
		  break;
		}
	    }
        }
    }
  return atom_occ;
}

Array2D<double> QMCDansWalkerInitialization::
dist_center(int atomic_charge, int n_e, int n_a, int n_b)
{
  Array2D<double> e_locs(n_e,3);

  if (n_e <= 2) 
    e_locs = dist_energy_level(atomic_charge,1,n_a,n_b);
      
  else if (n_e > 2 && n_e <= 10)
    {
      Array2D<double> level_1_of_2(2,3);
      level_1_of_2 = dist_energy_level(atomic_charge,1,1,1);
      for (int i=0; i<3; i++)
        {
	  e_locs(0,i) = level_1_of_2(0,i);
	  e_locs(n_a,i) = level_1_of_2(1,i);
	}
      Array2D<double> level_2_of_2(n_e-2,3);
      level_2_of_2 = dist_energy_level(atomic_charge,2,n_a-1,n_b-1); 

      for (int j=0; j<n_a-1; j++)
	for (int k=0; k<3; k++)
	  e_locs(1+j,k) = level_2_of_2(j,k);

      for (int m=0; m<n_b-1; m++)
	for (int n=0; n<3; n++)
	  e_locs(1+n_a+m,n) = level_2_of_2(n_a-1+m,n);
    }

  else if (n_e > 10 && n_e <=18)
    {
      Array2D<double> level_1_of_3(2,3); 
      level_1_of_3 = dist_energy_level(atomic_charge,1,1,1);
      
      for (int i=0; i<3; i++)
        {
	  e_locs(0,i) = level_1_of_3(0,i);
	  e_locs(n_a,i) = level_1_of_3(1,i);
	}

      Array2D<double> level_2_of_3(8,3);
      level_2_of_3 = dist_energy_level(atomic_charge,2,4,4); 

      for (int j=0; j<4; j++)
	for (int k=0; k<3; k++)
	  {
	    e_locs(1+j,k) = level_2_of_3(j,k);
	    e_locs(n_a+1+j,k) = level_2_of_3(4+j,k);
	  }

      Array2D<double> level_3_of_3(n_e-10,3);
      level_3_of_3 = dist_energy_level(atomic_charge,3,n_a-5,n_b-5);

      for (int j=0; j<n_a-5; j++)
	for (int k=0; k<3; k++)
	  e_locs(5+j,k) = level_3_of_3(j,k);

      for (int m=0; m<n_b-5; m++)
	for (int n=0; n<3; n++)
	  e_locs(n_a+5+m,n) = level_3_of_3(n_a-5+m,n);
    }
  return e_locs;
}

Array2D<double> QMCDansWalkerInitialization::
dist_energy_level(int Z, int n, int nalpha, int nbeta)
{
  // Distribute radial distances.

  Array1D<double> r_locs;
  r_locs = generateRadialDistances(Z,n,nalpha+nbeta);

  // The interval array determines which distribution will be retrieved from
  // the angle_dist_arrays file for each electron.  The intervals are assigned
  // based on the number of each type of electron.

  int a_interval = 0;
  int b_interval = 0;

  if (nalpha == 1) 
    {
      if (nbeta == 0) a_interval = 0;
      else a_interval = 1;  
    }
  else if (nalpha == 2) a_interval = 1;
  else if (nalpha == 3) a_interval = 5;
  else if (nalpha == 4) a_interval = 11;
  if (nbeta == 1) 
    {
      if (nalpha == 0) b_interval = 0;
      else if (nalpha == 1) b_interval = 2;
      else b_interval =  3;
    }
  else if (nbeta == 2) b_interval = 3;
  else if (nbeta == 3) b_interval = 8;
  else if (nbeta == 4) b_interval = 15;
  
  if (a_interval > 11 || a_interval < 0) 
    {
      cerr << "Error in alpha interval." << endl;
      exit (1);
    }
  if (b_interval > 15 || b_interval < 0)
    {
      cerr << "Error in beta interval." << endl;
      exit (1);
    }

  Array1D<double> a_phi(nalpha);
  Array1D<double> b_phi(nbeta);
  Array1D<double> a_theta(nalpha);
  Array1D<double> b_theta(nbeta);

  // Distributing the theta and phi points based on their distributions.

  for (int i=0; i<nalpha; i++)
    {
      a_phi(i) = generatePhiCoordinate(a_interval+i);
      a_theta(i) = generateThetaCoordinate(a_interval+i);
    }

  for (int j=0; j<nbeta; j++)
    {
      b_phi(j) = generatePhiCoordinate(b_interval+j);
      b_theta(j) = generateThetaCoordinate(b_interval+j);
    }

  // Now we convert the spherical coordinates into cartesians.

  Array2D<double> crds(nalpha+nbeta,3);

  for (int i=0; i<nalpha; i++)
    {
      crds(i,0)=r_locs(i)*sin(a_theta(i))*cos(a_phi(i));
      crds(i,1)=r_locs(i)*sin(a_theta(i))*sin(a_phi(i));
      crds(i,2)=r_locs(i)*cos(a_theta(i));
    }

  for (int j=0; j<nbeta; j++)
    {
      crds(nalpha+j,0)=r_locs(nalpha+j)*sin(b_theta(j))*cos(b_phi(j));
      crds(nalpha+j,1)=r_locs(nalpha+j)*sin(b_theta(j))*sin(b_phi(j));
      crds(nalpha+j,2)=r_locs(nalpha+j)*cos(b_theta(j));
    }

  // We use a random vector of unit length as our axis, and rotate all the 
  // points in the energy level about it by a random angle.

  double phi = ran1(&Input->flags.iseed)*2*PI;
  double theta = sindev(&Input->flags.iseed);
  double angle = ran1(&Input->flags.iseed)*2*PI;

  Array1D<double> axis(3);
  axis(0) = sin(theta)*cos(phi);
  axis(1) = sin(theta)*sin(phi);
  axis(2) = cos(theta);

  for (int i=0; i<(nalpha+nbeta); i++)
    {
      Array1D<double> temp_coords(3);
      for (int j=0; j<3; j++) temp_coords(j) = crds(i,j);
      temp_coords.rotate(axis,angle);
      for (int k=0; k<3; k++) crds(i,k) = temp_coords(k);
    }
  return crds;
}

double QMCDansWalkerInitialization::generatePhiCoordinate(int index)
{
  // These distributions are uniform with respect to phi.
  if (index == 0 || index == 1 || index == 2 || index == 5) 
    {
      return ran1(&Input->flags.iseed)*2*PI;
    }
  else
    {
      // If no spline has been made for the distribution, we make one and then
      // use it.
      if (phiSplinesMade(index) == 0)
	{
	  // Some of the distributions are identical in phi, so we generate
          // both splines at once.
	  if (index == 15 || index == 12)
	    {
	      Array1D<double> phi_array; 
	      phi_array = AngleDistributions::getPhiArray(12);
	      double derivative=(phi_array(1)+phi_array(20)-phi_array(19))*10;
	      phiSplines(12).initializeWithFunctionValues(x_array,phi_array,\
							derivative,derivative);
              phiSplines(15) = phiSplines(12);
              phiSplinesMade(12) = 1;
	      phiSplinesMade(15) = 1;
	    }
	  else if (index == 16 || index == 11)
	    {
	      Array1D<double> phi_array; 
              phi_array = AngleDistributions::getPhiArray(11);
	      double derivative=(phi_array(1)+phi_array(20)-phi_array(19))*10;
	      phiSplines(11).initializeWithFunctionValues(x_array,phi_array,\
							derivative,derivative);
              phiSplines(16) = phiSplines(11);
	      phiSplinesMade(11) = 1;
              phiSplinesMade(16) = 1;
	    }
	  else if (index == 17 || index == 14)
	    {
	      Array1D<double> phi_array; 
              phi_array = AngleDistributions::getPhiArray(14);
	      double derivative=(phi_array(1)+phi_array(20)-phi_array(19))*10;
	      phiSplines(14).initializeWithFunctionValues(x_array,phi_array,\
							derivative,derivative);
              phiSplines(17) = phiSplines(14);
	      phiSplinesMade(14) = 1;
              phiSplinesMade(17) = 1;
	    }
	  else if (index == 18 || index == 13)
	    {
	      Array1D<double> phi_array;
              phi_array = AngleDistributions::getPhiArray(13);
	      double derivative=(phi_array(1)+phi_array(20)-phi_array(19))*10;
	      phiSplines(13).initializeWithFunctionValues(x_array,phi_array,\
							derivative,derivative);
              phiSplines(18) = phiSplines(13);
	      phiSplinesMade(13) = 1;
	      phiSplinesMade(18) = 1;
	    }	      
	  else
	    {
	      Array1D<double> phi_array;
              phi_array = AngleDistributions::getPhiArray(index);
	      double derivative=(phi_array(1)+phi_array(20)-phi_array(19))*10;
	      phiSplines(index).initializeWithFunctionValues(x_array,\
                                              phi_array,derivative,derivative);
              phiSplinesMade(index) = 1;
	    }
	}
      phiSplines(index).evaluate(ran1(&Input->flags.iseed));
      return phiSplines(index).getFunctionValue();
    }
}

double QMCDansWalkerInitialization::generateThetaCoordinate(int index)
{
  // This distribution is uniform in theta.
  if (index == 0)
    {
      return sindev(&Input->flags.iseed);
    }
  else
    {
      // If the spline has not been made yet, we make it and use it.
      if (thetaSplinesMade(index) == 0)
	{
	  // Some distributions are identical in theta.  We generate both
          // splines at the same time.
	  if (index == 3 || index == 4)
	    {
	      Array1D<double> th_array;
              th_array = AngleDistributions::getThetaArray(3);
              double derivative = (th_array(1)+th_array(20)-th_array(19))*10;
              thetaSplines(3).initializeWithFunctionValues(x_array,th_array,\
							derivative,derivative);
	      thetaSplines(4) = thetaSplines(3);
	      thetaSplinesMade(3) = 1;
	      thetaSplinesMade(4) = 1;
	    }
	  else if (index == 6 || index == 7)
	    {
	      Array1D<double> th_array;
              th_array = AngleDistributions::getThetaArray(6);
              double derivative = (th_array(1)+th_array(20)-th_array(19))*10;
              thetaSplines(6).initializeWithFunctionValues(x_array,th_array,\
							derivative,derivative);
	      thetaSplines(7) = thetaSplines(6);
	      thetaSplinesMade(6) = 1;
	      thetaSplinesMade(7) = 1;
	    }
	  else if (index == 8 || index == 9 || index == 10)
	    {
	      Array1D<double> th_array = AngleDistributions::getThetaArray(8);
              double derivative = (th_array(1)+th_array(20)-th_array(19))*10;
              thetaSplines(8).initializeWithFunctionValues(x_array,th_array,\
							derivative,derivative);
	      thetaSplines(9) = thetaSplines(8);
              thetaSplines(10) = thetaSplines(8);
	      thetaSplinesMade(8) = 1;
	      thetaSplinesMade(9) = 1;
              thetaSplinesMade(10) = 1;
	    }
	  else if (index == 11 || index == 12 || index == 17 || index == 18)
	    {
	      Array1D<double> th_array; 
	      th_array = AngleDistributions::getThetaArray(11);
              double derivative = (th_array(1)+th_array(20)-th_array(19))*10;
              thetaSplines(11).initializeWithFunctionValues(x_array,th_array,\
							derivative,derivative);
	      thetaSplines(12) = thetaSplines(11);
	      thetaSplines(17) = thetaSplines(11);
	      thetaSplines(18) = thetaSplines(11);
	      thetaSplinesMade(11) = 1;
	      thetaSplinesMade(12) = 1;
	      thetaSplinesMade(17) = 1;
	      thetaSplinesMade(18) = 1;
	    }
	  else if (index == 13 || index == 14 || index == 15 || index == 16)
	    {
	      Array1D<double> th_array;
	      th_array = AngleDistributions::getThetaArray(13);
              double derivative = (th_array(1)+th_array(20)-th_array(19))*10;
              thetaSplines(13).initializeWithFunctionValues(x_array,th_array,\
							derivative,derivative);
	      thetaSplines(14) = thetaSplines(13);
	      thetaSplines(15) = thetaSplines(13);
	      thetaSplines(16) = thetaSplines(13);
	      thetaSplinesMade(13) = 1;
	      thetaSplinesMade(14) = 1;
	      thetaSplinesMade(15) = 1;
	      thetaSplinesMade(16) = 1;
	    }
          else
	    {
              Array1D<double> th_array;
	      th_array = AngleDistributions::getThetaArray(index);
	      double derivative = (th_array(1)+th_array(20)-th_array(19))*10;
	      thetaSplines(index).initializeWithFunctionValues(x_array,
					       th_array,derivative,derivative);
	      thetaSplinesMade(index) = 1;
	    }
	}
      thetaSplines(index).evaluate(ran1(&Input->flags.iseed));
      return thetaSplines(index).getFunctionValue();
    }
}

Array1D<double> QMCDansWalkerInitialization::\
                              generateRadialDistances(int Z, int n, int nelecs)
{
  int atomicNumberIndex = Z-1;
  int energyLevelIndex = n-1;
  if (radialSplinesMade(atomicNumberIndex,energyLevelIndex) == 0)
    {
      Array1D<double> r_array;
      r_array = RadialDistributions::getRadialArray(Z,n);
      double derivativeAtZero = r_array(1)*20;
      double derivativeAtOne = (r_array(20) - r_array(19))*20;
      radialSplines(atomicNumberIndex,energyLevelIndex).\
	initializeWithFunctionValues\
	                    (x_array,r_array,derivativeAtZero,derivativeAtOne);
      radialSplinesMade(atomicNumberIndex,energyLevelIndex) = 1;
    }
  Array1D<double> r_locs(nelecs);
  for (int i=0; i<nelecs; i++)
    {
      radialSplines(atomicNumberIndex,energyLevelIndex).evaluate\
                                                   (ran1(&Input->flags.iseed));
      r_locs(i) = radialSplines(atomicNumberIndex,energyLevelIndex).\
                                                            getFunctionValue();
    }
  return r_locs;
}

# undef PI
