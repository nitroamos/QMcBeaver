// This is dans_walker_initialization.  It creates initial configurations that 
// require very little time for the run to equilibrate.

#include "QMCDansWalkerInitialization.h"

#define PI 3.14159265359

QMCDansWalkerInitialization::QMCDansWalkerInitialization(QMCInput * INPUT)
{
  Input = INPUT;

  for (int i=0; i<Input->flags.Natoms; i++){
    if(Input->Molecule.Z(i) >= 18){
      clog << "Warning: QMCDansWalkerInitialization can only handle atoms with up to Z=18, but atom "
	   << i << " has Z=" << Input->Molecule.Z(i) << endl;
      clog << "         Instead of quitting, we'll proceed with our best guess..." << endl;
    }
  }

  if (!arraysInitialized)
    {
      arraysInitialized = true;
      initializeArrays();
    }
}

QMCDansWalkerInitialization::~QMCDansWalkerInitialization()
{
  arraysInitialized = false;
  phiSplines.deallocate();
  phiSplinesMade.deallocate();
  thetaSplines.deallocate();
  thetaSplinesMade.deallocate();
  radialSplines.deallocate();
  radialSplinesMade.deallocate();
  x_array.deallocate();
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

  // This will return the number of explicit electrons.  Electrons replaced by
  // pseudopotentials are not included in these numbers.
  int nelectrons = Input->WF.getNumberElectrons();
  int nbeta = Input->WF.getNumberElectrons(false);
  int nalpha = Input->WF.getNumberElectrons(true);

  // This array holds the numbers of total, alpha, and beta electrons on each
  // center.

  Array2D<int> ab_count(natoms,3);
  ab_count = assign_electrons_to_nuclei();

  // Redistribute electrons if there are charged centers.
  // We use Zeff here- the effective nuclear charge shielded by pseudo 
  // pseudo electrons because we are only placing the explicit electrons.

  int icharge_init, icharge, kcharge, partner, partnercharge, give, get;

  for (int i=0; i<natoms; i++)
    if (Input->Molecule.Zeff(i) != ab_count(i,0)) // If this atom is charged
      {
	icharge_init = Input->Molecule.Zeff(i) - ab_count(i,0);
        for (int j=0; j<abs(icharge_init); j++)
	  {
	    icharge = Input->Molecule.Zeff(i) - ab_count(i,0);
	    partner = -1;
	    give = -1;
	    get = -1;
	    partnercharge = 0;
	    for (int k=0; k<natoms; k++)  // Find an atom to exchange with
	      {
		if (k != i)
		  {
		    kcharge = Input->Molecule.Zeff(k) - ab_count(k,0);
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
			       (ab_count(l,1) <= ab_count(l,2)) )
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
			       (ab_count(m,1) >= ab_count(m,2)) )
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

  // Now we check to make sure no center has too many or too few electrons of 
  // one type.  We want to make sure each energy level is full before we start
  // filling the next one.  We assume the previous checks have made the atoms
  // neutral.  One example of what we want to prevent in this section is a Li
  // atom with three alphas.  It is neutral, but the first energy level was not
  // filled before we started occupying the second.

  // We are going to add the pseudo electrons back into ab_count for this
  // section.  Trying to tally the energy levels with some electrons missing 
  // would be too complicated.

  for (int i=0; i<natoms; i++)
    {
      if (Input->Molecule.usesPseudo(i) == true)
	{
	  int pseudo = Input->Molecule.Z(i) - Input->Molecule.Zeff(i);
	  if (pseudo == 0)
	    {
	      cerr << "ERROR: Atom " << i << " has a pseudopotential, but ";
	      cerr << "Zeff = Z = " << Input->Molecule.Z(i) << endl;
	      exit(1);
	    }
	  ab_count(i,0) += pseudo;
	  ab_count(i,1) += pseudo/2;
	  ab_count(i,2) += pseudo/2;
	}
    }

  for (int i=0; i<natoms; i++)
    {
      if (Input->Molecule.Z(i) <= 2)
        {
          // If this atom is in the first row, it should not have more than one
          // electron of one type.

          if (ab_count(i,1) > 1)
            {
              int extra_alphas = ab_count(i,1)-1;
              for (int p=0; p<extra_alphas; p++)
		{
		  for (int q=0; q<natoms; q++)
		    if (q != i)
		      if (ab_count(q,1) <= ab_count(q,2))
			{
			  ab_count(i,1) -= 1;
			  ab_count(q,1) += 1;
			  ab_count(i,2) += 1;
			  ab_count(q,2) -= 1;
			  break;
			}
		      else if (Input->Molecule.Z(q) - ab_count(q,0) > 0)
			{
			  ab_count(i,1) -= 1;
			  ab_count(i,0) -= 1;
			  ab_count(q,1) += 1;
			  ab_count(q,0) += 1;
			  break;
			}
		    }
	    }
          if (ab_count(i,2) > 1)
            {
              int extra_betas = ab_count(i,2)-1;
              for (int r=0; r<extra_betas; r++)
                for (int s=0; s<natoms; s++)
		  if (s != i)
		    {
		      if (ab_count(s,1) >= ab_count(s,2))
			{
			  ab_count(i,2) -= 1;
			  ab_count(i,1) += 1;
			  ab_count(s,2) += 1;
			  ab_count(s,1) -= 1;
			  break;
			}
		      else if (Input->Molecule.Z(s) - ab_count(s,0) > 0)
			{
			  ab_count(i,2) -= 1;
			  ab_count(i,0) -= 1;
			  ab_count(s,2) += 1;
			  ab_count(s,0) += 1;
			  break;
			}
		    }
	    }
	}
      if (Input->Molecule.Z(i) <= 10 && Input->Molecule.Z(i) > 2)
        {
          // If this atom is in the second row, it should have at least one but
          // no more than five electrons of one type.
          if (ab_count(i,1) > 5)
            {
              int extra_alphas = ab_count(i,1)-5;
              for (int j=0; j<extra_alphas; j++)
                for (int k=0; k<natoms; k++)
		  if (k != i)
		    {
		      if (ab_count(k,1) <= ab_count(k,2))
			{ 
			  ab_count(i,1) -= 1;
			  ab_count(k,1) += 1;
			  ab_count(i,2) += 1;
			  ab_count(k,2) -= 1;
			  break;
			}
		      else if (Input->Molecule.Z(k) - ab_count(k,0) > 0)
			{
			  ab_count(i,1) -= 1;
			  ab_count(i,0) -= 1;
			  ab_count(k,1) += 1;
			  ab_count(k,0) += 1;
			  break;
			}
		    }
            }
          if (ab_count(i,1) == 0)
            {
              // If there are zero alphas on this atom, the first energy level
              // is not filled, so we can't start the second one.  We need to
              // get an alpha from another atom.
              for (int j=0; j<natoms; j++)
                if (j != i)
		  {
		    if (ab_count(j,1) >= ab_count(j,2))
		      {
			ab_count(i,1) += 1;
			ab_count(j,1) -= 1;
			ab_count(i,2) -= 1;
			ab_count(j,2) += 1;
			break;
		      }
		    else if (Input->Molecule.Z(j) - ab_count(j,0) < 0)
		      {
			ab_count(i,1) += 1;
			ab_count(i,0) += 1;
			ab_count(j,1) -= 1;
			ab_count(j,0) -= 1;
			break;
		      }
		  }
            }
          if (ab_count(i,2) > 5)
            {
              int extra_betas = ab_count(i,2)-5;
              for (int m=0; m<extra_betas; m++)
                for (int n=0; n<natoms; n++)
                  if (n != i)
		    {
		      if (ab_count(n,1) >= ab_count(n,2))
			{
			  ab_count(i,2) -= 1;
			  ab_count(n,2) += 1;
			  ab_count(i,1) += 1;
			  ab_count(n,1) -= 1;
			  break;
			}
		      else if (Input->Molecule.Z(n) - ab_count(n,0) > 0)
			{
			  ab_count(i,1) -= 1;
			  ab_count(i,0) -= 1;
			  ab_count(n,1) += 1;
			  ab_count(n,0) += 1;
			  break;
			}
		    }
            }
          if (ab_count(i,2) == 0)
            {
              for (int j=0; j<natoms; j++)
                if (j != i)
                  {
		    if (ab_count(j,2) >= ab_count(j,1))
		      {
			ab_count(i,2) += 1;
			ab_count(j,2) -= 1;
			ab_count(i,1) -= 1;
			ab_count(j,1) += 1;
			break;
		      }
		    else if (Input->Molecule.Z(j) - ab_count(j,0) < 0)
		      {
			ab_count(i,2) += 1;
			ab_count(i,0) += 1;
			ab_count(j,2) -= 1;
			ab_count(j,0) -= 1;
			break;
		      }
		  }
            }
        }
      if (Input->Molecule.Z(i) <= 18 && Input->Molecule.Z(i) > 10)
        {
          // If the atom is in the third row, we want it to have at least five
          // but no more than nine electrons of each type.
          if (ab_count(i,1) > 9)
            {
              int extra_alphas = ab_count(i,1)-9;
              for (int j=0; j<extra_alphas; j++)
                for (int k=0; k<natoms; k++)
                  if (k != i)
		    {
		      if (ab_count(k,1) <= ab_count(k,2))
			{
			  ab_count(i,1) -= 1;
			  ab_count(k,1) += 1;
			  ab_count(i,2) += 1;
			  ab_count(k,2) -= 1;
			  break;
			}
		      else if (Input->Molecule.Z(k) - ab_count(k,0) > 0)
			{
			  ab_count(i,1) -= 1;
			  ab_count(i,0) -= 1;
			  ab_count(k,1) += 1;
			  ab_count(k,0) += 1;
			  break;
			}
		    }
            }
          if (ab_count(i,1) < 5)
            {
              int alphas_needed = 5-ab_count(i,1);
              for (int j=0; j<alphas_needed; j++)
                for (int k=0; k<natoms; k++)
                  if (k != i)
                    {
		      if (ab_count(k,1) >= ab_count(k,2))
			{
			  ab_count(i,1) += 1;
			  ab_count(k,1) -= 1;
			  ab_count(i,2) -= 1;
			  ab_count(k,2) += 1;
			  break;
			}
		      else if (Input->Molecule.Z(k) - ab_count(k,0) < 0)
			{
			  ab_count(i,1) += 1;
			  ab_count(i,0) += 1;
			  ab_count(k,1) -= 1;
			  ab_count(k,0) -= 1;
			  break;
			}
		    }
            }
          if (ab_count(i,2) > 9)
            {
              int extra_betas = ab_count(i,2)-9;
              for (int m=0; m<extra_betas; m++)
                for (int n=0; n<natoms; n++)
                  if (n != i)
                    {
		      if (ab_count(n,1) >= ab_count(n,2))
			{
			  ab_count(i,2) -= 1;
			  ab_count(n,2) += 1;
			  ab_count(i,1) += 1;
			  ab_count(n,1) -= 1;
			  break;
			}
		      else if (Input->Molecule.Z(n) - ab_count(n,0) > 0)
			{
			  ab_count(i,2) -= 1;
			  ab_count(i,0) -= 1;
			  ab_count(n,2) += 1;
			  ab_count(n,0) += 1;
			  break;
			}
		    }
            }
          if (ab_count(i,2) < 5)
            {
              int betas_needed = 5-ab_count(i,2);
              for (int j=0; j<betas_needed; j++)
                for (int k=0; k<natoms; k++)
                  if (k != i)
                    {
		      if (ab_count(k,2) >= ab_count(k,1))
			{
			  ab_count(i,2) += 1;
			  ab_count(k,2) -= 1;
			  ab_count(i,1) -= 1;
			  ab_count(k,1) += 1;
			  break;
			}
		      else if (Input->Molecule.Z(k) - ab_count(k,0) < 0)
			{
			  ab_count(i,2) += 1;
			  ab_count(i,0) += 1;
			  ab_count(k,2) -= 1;
			  ab_count(k,0) -= 1;
			}
		    }
            }
        }       
    }

  // Now we take the pseudo electrons back out of the ab_count.

  for (int i=0; i<natoms; i++)
    {
      if (Input->Molecule.usesPseudo(i) == true)
	{
	  int pseudo = Input->Molecule.Z(i) - Input->Molecule.Zeff(i);

	  ab_count(i,0) -= pseudo;
	  ab_count(i,1) -= pseudo/2;
	  ab_count(i,2) -= pseudo/2;
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

  int alpha_index = 0;
  int beta_index = 0;

  Array2D<double> temp_coords;
  int n_e,n_a,n_b;

  // The electronic coordinates are stored in this array, with the alpha
  // electrons first, then the betas.

  Array2D<double> R(nelectrons,3);

  for (int i=0; i<natoms; i++)
    {
      n_e = ab_count(i,0);
      n_a = ab_count(i,1);
      n_b = ab_count(i,2);

      if (n_e > 0)
	{
	  temp_coords.allocate(n_e,3);
	  temp_coords = dist_center(Input->Molecule.Z(i),Input->Molecule.Zeff(i),n_e,n_a,n_b);
	
	  for (int j=0; j<n_a; j++)
	    {
	      for (int k=0; k<3; k++)
		R(alpha_index,k) = temp_coords(j,k) + atom_centers(i,k);
	      alpha_index++;
	    }
      
	  for (int m=0; m<n_b; m++)
	    {
	      for (int n=0; n<3; n++)
		R(nalpha+beta_index,n)=temp_coords(n_a+m,n)+atom_centers(i,n);
	      beta_index++;
	    }
	}
    }
  return R;
}

Array2D<int> QMCDansWalkerInitialization::assign_electrons_to_nuclei()
{
  int Norbitals  = Input->WF.getNumberOrbitals();
  int Nbasisfunc = Input->WF.getNumberBasisFunctions();

  Array2D<qmcfloat> OrbitalScratch(Input->WF.OrbitalCoeffs);

  // Square MO coefficients.

  for (int i=0; i<Norbitals; i++)
    for (int j=0; j<Nbasisfunc; j++)
      {
	OrbitalScratch(i,j) = OrbitalScratch(i,j)*OrbitalScratch(i,j);
      }
  //"Normalize" MO by summing the squared coeffs and dividing

  double OrbitalSum;

  for (int i=0; i<Norbitals; i++)
    {
      OrbitalSum = 0.0;
      for(int j=0; j<Nbasisfunc; j++)
	{ 
	  OrbitalSum += OrbitalScratch(i,j);
	}
      for(int j=0; j<Nbasisfunc; j++)
	{
	  OrbitalScratch(i,j) = OrbitalScratch(i,j)/OrbitalSum;
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
  rv = ran.unidev();

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
      if (Input->WF.AlphaOccupation(determinant_index,i) != Input->WF.getUnusedIndicator())
	{
          sum = 0.0;
          rv = ran.unidev();

          for (int j=0; j<Nbasisfunc; j++)
            {
              sum += OrbitalScratch(i,j);
              if (sum >= rv) 
		{
		  atom_occ(OrbPos(j),0) += 1;
                  atom_occ(OrbPos(j),1) += 1; 
                  break;
		}
	    }
	}

      if (Input->WF.BetaOccupation(determinant_index,i) != Input->WF.getUnusedIndicator())
	{
	  sum = 0.0;
          rv = ran.unidev();

	  for (int k=0; k<Nbasisfunc; k++)
	    {
	      sum += OrbitalScratch(i,k);
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
dist_center(int atomic_charge, int eff_charge, int n_e, int n_a, int n_b)
{
  // The positions of the electrons are held in this array
  Array2D<double> e_locs(n_e,3);
  e_locs = 1e20;

  Array2D<double> temp;

  int alpha_index = 0;
  int beta_index  = n_a;

  // These are how many electrons we still have to place
  int elecs_left = n_e;
  int alphas_left = n_a;
  int betas_left  = n_b;

  // These are the number of electrons being placed in an energy level
  int alphas_now  = 0;
  int betas_now   = 0;

  // This keeps track of which energy level we are working on
  int energy_level = 0;

  // If a pseudopotential is used on this atom, we skip the core energy levels
  if (atomic_charge == eff_charge) // no pseudopotential
    energy_level = 1;
  else if (atomic_charge-eff_charge == 2) // He core
    energy_level = 2;
  else if (atomic_charge-eff_charge == 10) // Ne core
    energy_level = 3;
  else if (atomic_charge-eff_charge == 18) // Ar core
    energy_level = 4;
  else if (atomic_charge-eff_charge == 36) // Kr core
    energy_level = 5;
  else if (atomic_charge-eff_charge == 54) // Xe core
    energy_level = 6;
  else if (atomic_charge-eff_charge == 86) // Rn core
    energy_level = 7;

  while (elecs_left > 0)
    {
      if (energy_level == 1)
	{
	  if (alphas_left >= 1)
	    alphas_now = 1;
	  else
	    alphas_now = 0;

	  if (betas_left >= 1)
	    betas_now = 1;
	  else
	    betas_now = 0;

	}
      else if (energy_level == 2 || energy_level == 3)
	{
	  if (alphas_left >= 4)
	    alphas_now = 4;
	  else
	    alphas_now = alphas_left;

	  if (betas_left >= 4)
	    betas_now = 4;
	  else
	    betas_now = betas_left;
	}
      else if (energy_level == 4 || energy_level == 5)
	{
	  if (alphas_left >= 9)
	    alphas_now = 9;
	  else
	    alphas_now = alphas_left;
	  
	  if (betas_left >= 9)
	    betas_now = 9;
	  else
	    betas_now = betas_left;
	}
      else if (energy_level == 6 || energy_level == 7)
	{
	  if (alphas_left >= 16)
	    alphas_now = 16;
	  else
	    alphas_now = alphas_left;

	  if (betas_left >= 16)
	    betas_now = 16;
	  else
	    betas_now = betas_left;
	}

      temp = dist_energy_level(atomic_charge,energy_level,alphas_now,betas_now);

      for (int i=0; i<alphas_now; i++)
	{
	  for (int j=0; j<3; j++)
	    e_locs(alpha_index,j) = temp(i,j);
	  alpha_index++;
	}
      for (int i=0; i<betas_now; i++)
	{
	  for (int j=0; j<3; j++)
	    e_locs(beta_index,j) = temp(alphas_now+i,j);
	  beta_index++;
	}

      alphas_left -= alphas_now;
      betas_left -= betas_now;
      elecs_left -= alphas_now + betas_now;
      energy_level++;
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

  double phi = ran.unidev()*2*PI;
  double theta = ran.sindev();
  double angle = ran.unidev()*2*PI;

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
      return ran.unidev()*2*PI;
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
      phiSplines(index).evaluate(ran.unidev());
      return phiSplines(index).getFunctionValue();
    }
}

double QMCDansWalkerInitialization::generateThetaCoordinate(int index)
{
  // This distribution is uniform in theta.
  if (index == 0)
    {
      return ran.sindev();
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
      thetaSplines(index).evaluate(ran.unidev());
      return thetaSplines(index).getFunctionValue();
    }
}

Array1D<double> QMCDansWalkerInitialization::\
                              generateRadialDistances(int Z, int n, int nelecs)
{
  int atomicNumberIndex = min(17,Z-1);
  int energyLevelIndex =  min(2,n-1);

  if(atomicNumberIndex > radialSplinesMade.dim1()){
    clog << "Error: QMCDansWalkerInitialization can not initialize atom with atomic number: " << Z << endl;
    exit(0);
  }

  if (radialSplinesMade(atomicNumberIndex,energyLevelIndex) == 0)
    {
      Array1D<double> r_array;
      r_array = RadialDistributions::getRadialArray(atomicNumberIndex+1,
						    energyLevelIndex+1);
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
                                                   (ran.unidev());
      r_locs(i) = radialSplines(atomicNumberIndex,energyLevelIndex).\
                                                            getFunctionValue();
    }
  return r_locs;
}

# undef PI
