// This is dans_walker_initialization.  It creates initial configurations that 
// require very little time for the run to equilibrate.

#include "QMCDansWalkerInitialization.h"

#define PI 3.14159265359

QMCDansWalkerInitialization::QMCDansWalkerInitialization(QMCInput * INPUT)
{
  Input = INPUT;
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

  // Redistribute electrons if there are charged centers or if there are more
  // than 5 of one type around a center with Z <= 10.

  int needed_e, extra_e, temp_na, temp_nb;

  for (int i=0; i<natoms; i++)
    if (Input->Molecule.Z(i) < ab_count(i,0))  // If the center is negative.
      {
	extra_e = ab_count(i,0) - Input->Molecule.Z(i);
	temp_na = ab_count(i,1);
	temp_nb = ab_count(i,2);
        for (int j=0; j<extra_e; j++)
	  for (int k=0; k<natoms; k++)  // Find a center that is positive.
	    if (Input->Molecule.Z(k) > ab_count(k,0))
      	      {
                if (temp_na == temp_nb)  // if na=nb on -
		  {
		    if (ab_count(k,1) <= ab_count(k,2)) // and na <= nb on +
		      {  // transfer an alpha
			ab_count(i,0) -= 1;
			ab_count(i,1) -= 1;
			ab_count(k,0) += 1;
			ab_count(k,1) += 1;
			temp_na--;
		        break;
		      }
		    else if (ab_count(k,1) > ab_count(k,2))  // if na>nb on +
		      {  // transfer a beta
			ab_count(i,0) -= 1;
			ab_count(i,2) -= 1;
			ab_count(k,0) += 1;
			ab_count(k,2) += 1;
			temp_nb--;
			break;
		      }
		  }
		else if (temp_na > temp_nb) // if na>nb on -
		  {
		    if (ab_count(k,1) <= ab_count(k,2)) // and na<=nb on +
		      {  // transfer an alpha
			ab_count(i,0) -= 1;
			ab_count(i,1) -= 1;
			ab_count(k,0) += 1;
			ab_count(k,1) += 1;
			temp_na--;
			break;
		      }
		    else if (ab_count(k,1) > ab_count(k,2)) // if na>nb on +
		      // find a neutral center where nb>na, and trade an alpha
		      // from - for a beta for +.
		      {
		        for (int l=0; l<natoms; l++) 
			  if (Input->Molecule.Z(l) == ab_count(l,0) &&
			                        ab_count(l,1) < ab_count(l,2))
			    {  
			      ab_count(i,0) -= 1;
			      ab_count(i,1) -= 1;
			      ab_count(l,1) += 1;
			      ab_count(l,2) -= 1;
			      ab_count(k,0) += 1;
			      ab_count(k,2) += 1;
			      temp_na--;
			      break;
			    }
		        break;
		      }
		  }
		else if (temp_na < temp_nb)  // if nb>na on -
		  {
		    if (ab_count(k,1) >= ab_count(k,2)) // and na>=nb on +
		      {  // transfer a beta
			ab_count(i,0) -= 1;
			ab_count(i,2) -= 1;
			ab_count(k,0) += 1;
			ab_count(k,2) += 1;
			temp_nb--;
			break;
		      }
		    else if (ab_count(k,1) < ab_count(k,2)) // if nb>na on +
		      // find a neutral center where na>nb, and trade a beta
		      // from - for an alpha for +
		      {
			for (int m=0; m<natoms; m++)
			  if (Input->Molecule.Z(m) == ab_count(m,0) &&
			                         ab_count(m,1) > ab_count(m,2))
			    {
			      ab_count(i,0) -= 1;
			      ab_count(i,2) -= 1;
			      ab_count(m,2) += 1;
			      ab_count(m,1) -= 1;
			      ab_count(k,0) += 1;
			      ab_count(k,1) += 1;
			      temp_nb--;
			      break;
			    }
			break;
		      }
		  }
	      }
      }

  for (int i=0; i<natoms; i++)
    if (Input->Molecule.Z(i) > ab_count(i,0))  // If the center is positive.
      {
	needed_e = Input->Molecule.Z(i) - ab_count(i,0);
	temp_na = ab_count(i,1);
	temp_nb = ab_count(i,2);
        for (int j=0; j<needed_e; j++)
	  for (int k=0; k<natoms; k++)  // Find a center that is negative.
	    if (Input->Molecule.Z(k) < ab_count(k,0))
      	      {
                if (temp_na == temp_nb)  // if na=nb on +
		  {
		    if (ab_count(k,1) <= ab_count(k,2)) // and na <= nb on -
		      {  // transfer a beta
			ab_count(i,0) += 1;
			ab_count(i,2) += 1;
			ab_count(k,0) -= 1;
			ab_count(k,2) -= 1;
			temp_nb++;
		        break;
		      }
		    else if (ab_count(k,1) > ab_count(k,2))  // if na>nb on -
		      {  // transfer an alpha
			ab_count(i,0) += 1;
			ab_count(i,1) += 1;
			ab_count(k,0) -= 1;
			ab_count(k,1) -= 1;
			temp_na++;
			break;
		      }
		  }
		else if (temp_na > temp_nb) // if na>nb on +
		  {
		    if (ab_count(k,1) <= ab_count(k,2)) // and na<=nb on -
		      {  // transfer a beta
			ab_count(i,0) += 1;
			ab_count(i,2) += 1;
			ab_count(k,0) -= 1;
			ab_count(k,2) -= 1;
			temp_na++;
			break;
		      }
		    else if (ab_count(k,1) > ab_count(k,2)) // if na>nb on -
		      // find a neutral center where nb>na, and trade an alpha
		      // from - for a beta for +.
		      {
		        for (int l=0; l<natoms; l++) 
			  if (Input->Molecule.Z(l) == ab_count(l,0) &&
			                        ab_count(l,1) < ab_count(l,2))
			    {  
			      ab_count(i,0) += 1;
			      ab_count(i,2) += 1;
			      ab_count(l,1) += 1;
			      ab_count(l,2) -= 1;
			      ab_count(k,0) -= 1;
			      ab_count(k,1) -= 1;
			      temp_nb++;
			      break;
			    }
		        break;
		      }
		  }
		else if (temp_na < temp_nb)  // if nb>na on +
		  {
		    if (ab_count(k,1) >= ab_count(k,2)) // and na>=nb on -
		      {  // transfer an alpha
			ab_count(i,0) += 1;
			ab_count(i,1) += 1;
			ab_count(k,0) -= 1;
			ab_count(k,1) -= 1;
			temp_na++;
			break;
		      }
		    else if (ab_count(k,1) < ab_count(k,2)) // if nb>na on -
		      // find a neutral center where na>nb, and trade a beta
		      // from - for an alpha for +
		      {
			for (int m=0; m<natoms; m++)
			  if (Input->Molecule.Z(m) == ab_count(m,0) &&
			                         ab_count(m,1) > ab_count(m,2))
			    {
			      ab_count(i,0) += 1;
			      ab_count(i,1) += 1;
			      ab_count(m,2) += 1;
			      ab_count(m,1) -= 1;
			      ab_count(k,0) -= 1;
			      ab_count(k,2) -= 1;
			      temp_na++;
			      break;
			    }
			break;
		      }
		  }
	      }
      }

  // Now we check to make sure no second row atom has too many electrons of the
  // same type.

  for (int i=0; i<natoms; i++)
    {
      if (Input->Molecule.Z(i) <= 18 && ab_count(i,1) > 9)
	{
	  int extra_alphas = ab_count(i,1)-9;
	  for (int j=0; j<extra_alphas; j++)
	    for (int k=0; k<natoms; k++)
	      if (ab_count(k,1) < ab_count(k,2))
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
	      {
		if (ab_count(n,1) > ab_count(n,2))
		  {
		    ab_count(i,2) -= 1;
		    ab_count(n,2) += 1;
		    ab_count(i,1) += 1;
		    ab_count(n,1) -= 1;
		    break;
		  }
	      }
	}
      if (Input->Molecule.Z(i) <= 10 && ab_count(i,1) > 5)
	{
	  int extra_alphas = ab_count(i,1)-5;
	  for (int j=0; j<extra_alphas; j++)
	    for (int k=0; k<natoms; k++)
	      if (ab_count(k,1) < ab_count(k,2))
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
	      {
		if (ab_count(n,1) > ab_count(n,2))
		  {
		    ab_count(i,2) -= 1;
		    ab_count(n,2) += 1;
		    ab_count(i,1) += 1;
		    ab_count(n,1) -= 1;
		    break;
		  }
	      }
	}
      if (Input->Molecule.Z(i) <= 2 && ab_count(i,1) > 1)
	{
	  int extra_alphas = ab_count(i,1)-1;
	  for (int p=0; p<extra_alphas; p++)
	    for (int q=0; q<natoms; q++)
	      {
		if ( Input->Molecule.Z(q) > 2 && ab_count(q,1) < 5)
		  {
		    ab_count(i,0) -= 1;
		    ab_count(i,1) -= 1;
		    ab_count(q,0) += 1;
		    ab_count(q,1) += 1;
		    break;
		  }
	      }
	}
      if (Input->Molecule.Z(i) <= 2 && ab_count(i,2) > 1)
	{
	  int extra_betas = ab_count(i,2)-1;
	  for (int r=0; r<extra_betas; r++)
	    for (int s=0; s<natoms; s++)
	      {
		if (Input->Molecule.Z(s) > 2 && ab_count(s,2) < 5)
		  {
		    ab_count(i,0) -= 1;
		    ab_count(i,2) -= 1;
		    ab_count(s,0) += 1;
		    ab_count(s,2) += 1;
		    break;
		  }
	      }
	}
    }

  /*
  cout << "Final ab_count array:" << endl;

  for (int i=0; i<natoms; i++)
    {
      cout << "\t" << "Z = " << Input->Molecule.Z(i);
      cout << "\t" << "alpha = " << ab_count(i,1);
      cout << "\t" << "beta = " << ab_count(i,2) << endl;
    }   	      
  */  

  // Check to see that all electrons have been assigned.

  int check_elec = 0;
  int check_alpha = 0;
  int check_beta = 0;

  for (int i=0; i<natoms; i++)
    {
      check_elec += ab_count(i,0);
      check_alpha += ab_count(i,1);
      check_beta += ab_count(i,2);
    }

  if (check_elec != nelectrons)
    { cerr << "Electron counting has gone bad." << endl; exit(1); }

  if (check_alpha != nalpha)
    { cerr << "Alpha counting has gone bad." << endl; exit(1); }
  
  if (check_beta != nbeta)
    { cerr << "Beta counting has gone bad." << endl; exit(1); }

  // The electronic coordinates are stored in this array, with the alpha
  // electrons first, then the betas.

  Array2D<double> R(nelectrons,3);

  int alpha_index = 0;
  int beta_index = nalpha;

  Array2D<double> temp_coords;
  int n_e, n_a, n_b;

  for (int i=0; i<natoms; i++)
    {
      n_e = ab_count(i,0);
      n_a = ab_count(i,1);
      n_b = ab_count(i,2);

      temp_coords.allocate(n_e,3);
      temp_coords = dist_center(Input->Molecule.Z(i),n_e,n_a,n_b);

      for (int j=0; j<n_a; j++)
	{
      	  for (int k=0; k<3; k++)
	    {
	      R(alpha_index,k) = temp_coords(j,k) + atom_centers(i,k);
	    }
	  alpha_index++;
	}
      
      for (int m=0; m<n_b; m++)
	{
	  for (int n=0; n<3; n++)
	    {
	      R(beta_index,n) = temp_coords(n_a+m,n) + atom_centers(i,n);
	    }
	  beta_index++;
	}
    }

  return R;
}

Array2D<int> QMCDansWalkerInitialization::assign_electrons_to_nuclei()
{
  int Norbitals  = Input->flags.Norbitals;
  int Nbasisfunc = Input->flags.Nbasisfunc;

  Array2D<double> Scratch(Input->WF.Coeffs);

  // Square MO coefficients.

  for (int i=0; i<Nbasisfunc; i++)
    {
      for (int j=0; j<Norbitals; j++)
        {
          Scratch(i,j) = Scratch(i,j)*Scratch(i,j);
        }
    }

  //"Normalize" MO by summing the squared coeffs and dividing

  double sum;

  for (int i=0; i<Norbitals; i++)
    {
      sum = 0.0;
      for(int j=0; j<Nbasisfunc; j++) sum += Scratch(j,i);

      for(int j=0; j<Nbasisfunc; j++)
        {
          Scratch(j,i) = Scratch(j,i)/sum;
        }
    }

  // Create an Array that links the index of an basis function to the atom it
  // is associated with

  int Natoms = Input->flags.Natoms;

  Array1D<int> OrbPos(Nbasisfunc);

  int i = 0;
  for (int j=0; j<Natoms; j++)
    {
      for(int k=0; k<(Input->BF.getNumberBasisFunctions(j)); k++)
        {
          OrbPos(i) = j;
          i++;
        }
    }

  //Initialize the array to be returned.

  Array2D <int> atom_occ(Natoms,3);

  for (int i=0; i<Natoms; i++)
    for (int j=0; j<3; j++) atom_occ(i,j) = 0;

  double rv;  //uniform [0,1] random variable
  for (int i=0; i<Norbitals; i++)
    {
      if (Input->WF.AlphaOccupation(i) == 1)
	{
          sum = 0.0;
          rv = ran1(&Input->flags.iseed);

          for (int j=0; j<Nbasisfunc; j++)
            {
              sum += Scratch(j,i);
              if (sum >= rv) 
		{
		  atom_occ(OrbPos(j),0) += 1;
                  atom_occ(OrbPos(j),1) += 1; 
                  break;
		}
	    }
	}

      if (Input->WF.BetaOccupation(i) == 1)
	{
	  sum = 0.0;
          rv = ran1(&Input->flags.iseed);

	  for (int k=0; k<Nbasisfunc; k++)
	    {
	      sum += Scratch(k,i);
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
    {
      e_locs = dist_energy_level(atomic_charge,1,n_a,n_b);
    }
      
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
  // Get the radial array.
  ifstream r_file ("radial_dist_arrays");
  r_file.clear();
  if (!r_file)
    { cerr << "Error opening radial array file." << endl; exit(1); }
  string temp_string;
  int length = 0, skip = 0;
 find_ampersand:
  do r_file >> temp_string; while (temp_string != "&");
  r_file >> temp_string;
  if (atoi(temp_string.c_str()) != Z) goto find_ampersand;
  else if (atoi(temp_string.c_str()) == Z) 
    {
     find_n:
      r_file >> temp_string;
      if (atoi(temp_string.c_str()) != n) 
	{
          r_file >> temp_string;
          skip = atoi(temp_string.c_str());
          for (int i=0; i<2*skip; i++) r_file >> temp_string;
          goto find_n; 
	}     
      else if (atoi(temp_string.c_str()) == n) 
	{
          r_file >> temp_string;
          length = atoi(temp_string.c_str()); 
	}
      else
	{ 
	  cerr << "Error- no energy level " << n << "for Z = " << Z << endl;
	  exit(1);
	}
    }

  // Read the array from the file.

  Array2D<double> r_array(length,2);
  for (int i=0; i<length; i++) 
    {
      r_file >> temp_string;
      r_array(i,0) = atof(temp_string.c_str());
      r_file >> temp_string;
      r_array(i,1) = atof(temp_string.c_str()); 
    }

  r_file.close();	

  // Distribute the radial points with respect to the distribution.

  Array1D<double> r_locs(nalpha+nbeta); 
  r_locs  = dist_wrt_array(nalpha+nbeta, r_array);

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

  // The angle distributions are stored in 2D arrays.

  Array2D<double> a_phi_array(1001,nalpha+1);
  Array2D<double> a_theta_array(501,nalpha+1);
  Array2D<double> b_phi_array(1001,nbeta+1);
  Array2D<double> b_theta_array(501,nbeta+1);

  ifstream a_file("angle_dist_arrays");
  a_file.clear();
  if (!a_file)
    { cerr << "Error opening angle array file." << endl; exit(1); }
  do a_file >> temp_string; while (temp_string != "&");
  for (int i=0; i<1001; i++)
    {
      a_file >> temp_string;
      a_phi_array(i,0) = atof(temp_string.c_str());
      for (int j=0; j<a_interval; j++) a_file >> temp_string;
      for (int k=0; k<nalpha; k++)
	{
	  a_file >> temp_string;
	  a_phi_array(i,k+1) = atof(temp_string.c_str());
	}
      for (int l=0; l<(19-a_interval-nalpha); l++) a_file >> temp_string;
    }
  do a_file >> temp_string; while (temp_string != "&");
  for (int l=0; l<501; l++)
    {
      a_file >> temp_string;
      a_theta_array(l,0) = atof(temp_string.c_str());
      for (int m=0; m<a_interval; m++) a_file >> temp_string;
      for (int n=0; n<nalpha; n++)
	{
	  a_file >> temp_string;
	  a_theta_array(l,n+1) = atof(temp_string.c_str());
	}
      for (int p=0; p<(19-a_interval-nalpha); p++) a_file >> temp_string;
    }
  a_file.close();

  ifstream b_file("angle_dist_arrays");
  b_file.clear();
  if (!b_file)
    { cerr << "Error opening angle array file." << endl; exit(1); }
  do b_file >> temp_string; while (temp_string != "&");
  for (int i=0; i<1001; i++)
    {
      b_file >> temp_string;
      b_phi_array(i,0) = atof(temp_string.c_str());
      for (int j=0; j<b_interval; j++) b_file >> temp_string;
      for (int k=0; k<nbeta; k++)
	{
	  b_file >> temp_string;
	  b_phi_array(i,k+1) = atof(temp_string.c_str());
	}
      for (int l=0; l<(19-b_interval-nbeta); l++) b_file >> temp_string;
    }
  do b_file >> temp_string; while (temp_string != "&");
  for (int l=0; l<501; l++)
    {
      b_file >> temp_string;
      b_theta_array(l,0) = atof(temp_string.c_str());
      for (int m=0; m<b_interval; m++) b_file >> temp_string;
      for (int n=0; n<nbeta; n++)
	{
	  b_file >> temp_string;
	  b_theta_array(l,n+1) = atof(temp_string.c_str());
	}
      for (int p=0; p<(19-b_interval-nbeta); p++) b_file >> temp_string;
    }
  b_file.close();

  Array1D<double> a_phi(nalpha);
  Array1D<double> b_phi(nbeta);
  Array1D<double> a_theta(nalpha);
  Array1D<double> b_theta(nbeta);

  // Distributing the theta and phi points based on their distributions.

  for (int i=0; i<nalpha; i++)
    {
      a_phi(i) = dist_wrt_array(a_phi_array,i+1); 
      a_theta(i) = dist_wrt_array(a_theta_array,i+1);
    }

  for (int j=0; j<nbeta; j++)
    {
      b_phi(j) = dist_wrt_array(b_phi_array,j+1); 
      b_theta(j) = dist_wrt_array(b_theta_array,j+1);
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
  double theta = ran1(&Input->flags.iseed)*PI;
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

Array1D<double> QMCDansWalkerInitialization::
dist_wrt_array(int npts, Array2D<double> dist)
{
  // Check to make sure this is a good distribution.
  // y values should go from 0 to 1, both x and y should monotonically 
  // increase.

  int length = dist.dim1();
  if (dist(0,1) != 0.0 || dist(length-1,1) > 1.0) 
    {
      cerr << "Bad distribution passed to dist_wrt_array." << endl;
      exit(1);
    }

  Array1D<double> x_values(npts);
  double rv=0.0, remainder;
  int lo, hi, mid;
  double dx = (dist(length-1,0)-dist(0,0))/length;
  for (int i=0; i<npts; i++)
    { 
      rv = ran1(&Input->flags.iseed);

      // If the random number is greater than the greatest y value, the point
      // is distributed exponentially.  

      if (rv > dist(length-1,1))
      {
        remainder=(rv-dist(length-1,1))/(1-dist(length-1,1));
        x_values(i) = dist(length-1,0) - dx*log(1-remainder);
      }

      // Otherwise a binary search algorithm determines which x value 
      // corresponds to the random value.

      else
      {
	lo = 0;
        hi = length-1;
	mid = (lo + hi)/2;
        while (mid != lo)
	{
	  if (rv < dist(mid,1)) { hi = mid; }
	  else if (rv >= dist(mid,1)) { lo = mid; }
          mid = (lo + hi)/2;
	}
	remainder = (rv-dist(lo,1))/(dist(hi,1)-dist(lo,1));
	x_values(i) = dist(lo,0) + dx*remainder;
      }
    }
  return x_values;
}

   
double QMCDansWalkerInitialization::
dist_wrt_array(Array2D<double> dist, int index)
{
  int length = dist.dim1();
  if (dist(0,index) != 0.0 || dist(length-1,index) > 1.0) 
    {
      cerr << "Bad distribution passed to dist_wrt_array." << endl;
      exit(1);
    }

  double x_value;
  double rv=0.0, remainder;
  int lo, hi, mid;
  double dx = (dist(length-1,0)-dist(0,0))/length;
  rv = ran1(&Input->flags.iseed);
  if (rv > dist(length-1,index))
    {
      remainder=(rv-dist(length-1,index))/(1-dist(length-1,index));
      x_value = dist(length-1,0) - dx*log(1-remainder);
    }
  else
    {
      lo = 0;
      hi = length-1;
      mid = (lo + hi)/2;
    while (mid != lo)
      {
        if (rv < dist(mid,index)) { hi = mid; }
        else if (rv >= dist(mid,index)) { lo = mid; }
        mid = (lo + hi)/2;
      }
      remainder = (rv-dist(lo,index))/(dist(hi,index)-dist(lo,index));
      x_value = dist(lo,0) + dx*remainder;
    }
  return x_value;
}
