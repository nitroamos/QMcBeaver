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

#include "QMCWavefunction.h"
#include "QMCInput.h"
#include <iomanip>

QMCWavefunction::QMCWavefunction()
{
  Norbitals = 0; 
  Nbasisfunc = 0; 
  Nalpha = 0; 
  Nbeta = 0; 
  Ncharge = 0;
  Nelectrons = 0;
  Ndeterminants = 0;
  factor = 1.0;
  unusedIndicator = 0;
  trialFunctionType = "restricted";
}

int QMCWavefunction::getNumberOrbitals()
{
  return Norbitals;
}

int QMCWavefunction::getNumberOrbitals(Array2D<int> & Occupation)
{
  Array1D<int> count(Occupation.dim1());
  count = 0;
  for(int ci=0; ci<Occupation.dim1(); ci++)
    for(int o=0; o<Occupation.dim2(); o++)
      if(Occupation(ci,o) != unusedIndicator)
	count(ci) = count(ci) + 1;

  for(int ci=0; ci<count.dim1(); ci++)
    if(count(ci) != count(0))
      {
	clog << "Error: number of orbitals doesn't match\n";
	clog << " count = " << count << endl;
	exit(0);
      }

  return count(0);
}

int QMCWavefunction::getNumberActiveOrbitals()
{
  int Nactiveorbs  = 0;
  for(int o=0; o<Norbitals; o++)
    {
      bool Aused = false;
      bool Bused = false;
      for(int ci=0; ci<Ndeterminants; ci++)
	{
	  if(AlphaOccupation(ci,o) != unusedIndicator)
	    Aused = true;
	  if(BetaOccupation(ci,o) != unusedIndicator)
	    Bused = true;
	}
      if(Aused || Bused) Nactiveorbs++;
    }
  return Nactiveorbs;
}

int QMCWavefunction::getNumberActiveOrbitals(Array2D<int> & Occupation)
{
  int Nactiveorbs = 0;
  for(int o=0; o<Norbitals; o++)
    {
      bool used = false;
      for(int ci=0; ci<Ndeterminants; ci++)
	{
	  if(Occupation(ci,o) != unusedIndicator)
	    used = true;
	}
      if(used)          Nactiveorbs++;
    }
  return Nactiveorbs;
}

int QMCWavefunction::getNumberBasisFunctions()
{
  return Nbasisfunc;
}

int QMCWavefunction::getNumberElectrons(bool isAlpha)
{
  if(isAlpha)
    return Nalpha;
  else
    return Nbeta;
}

int QMCWavefunction::getNumberElectrons()
{
  return Nelectrons;
}

int QMCWavefunction::getNumberDeterminants()
{
  return Ndeterminants;
}

void QMCWavefunction::scaleCoeffs(double scaleFactor)
{
  factor *= scaleFactor;
  OrbitalCoeffs *= scaleFactor;
}

int QMCWavefunction::getUnusedIndicator()
{
  return unusedIndicator;
}

void QMCWavefunction::sortOccupations(bool ordered)
{
  int currentUnusedIndicator = getUnusedIndicator();

  unusedIndicator = 0;
  if(ordered) unusedIndicator = -1;

  for(int ci=0; ci<Ndeterminants; ci++)
    {
      int numOrbA = 0;
      int numOrbB = 0;
      for(int o=0; o<Norbitals; o++)
	{
	  if(AlphaOccupation(ci,o) == currentUnusedIndicator)
	    AlphaOccupation(ci,o) = unusedIndicator;
	  else
	    if(ordered)
	      AlphaOccupation(ci,o) = numOrbA++;
	    else
	      AlphaOccupation(ci,o) = 1;

	  if(BetaOccupation(ci,o) == currentUnusedIndicator)
	    BetaOccupation(ci,o) = unusedIndicator;
	  else
	    if(ordered)
	      BetaOccupation(ci,o) = numOrbB++;
	    else
	      BetaOccupation(ci,o) = 1;
	}
    }
}

void QMCWavefunction::unlinkOrbitals()
{
  /*
    First we need to figure out how many and which
    orbitals will need to be duplicated.
  */
  Array1D<bool> duplicateOrbital(Norbitals);
  int finalNorbitals = Norbitals;
  for(int o=0; o<Norbitals; o++)
    {
      bool alphaUses = false;
      bool betaUses  = false;
      for(int ci=0; ci<Ndeterminants; ci++)
	{
	  if(AlphaOccupation(ci,o) != unusedIndicator)
	    alphaUses = true;
	  if(BetaOccupation(ci,o) != unusedIndicator)
	    betaUses  = true;
	}
      //If they both use a particular orbital, then we need
      //to decouple them.
      if(alphaUses && betaUses)
	{
	  duplicateOrbital(o) = true;
	  finalNorbitals++;
	} else {
	  duplicateOrbital(o) = false;
	}
    }

  //The orbitals might already be fully decoupled.
  if(finalNorbitals == Norbitals)
    return;

  Array2D<qmcfloat> newOrbitalCoeffs(finalNorbitals,Nbasisfunc);
  Array2D<int>      newAlphaOccupation(Ndeterminants, finalNorbitals);
  Array2D<int>      newBetaOccupation(Ndeterminants, finalNorbitals);

  int newIndex = Norbitals;
  for(int o=0; o<Norbitals; o++)
    {
      newOrbitalCoeffs.setRows(o,o,1,OrbitalCoeffs);
      if(duplicateOrbital(o))
	{
	  newOrbitalCoeffs.setRows(newIndex,o,1,OrbitalCoeffs);
	  for(int ci=0; ci<Ndeterminants; ci++)
	    {
	      newAlphaOccupation(ci,o)        = AlphaOccupation(ci,o);
	      newBetaOccupation(ci,o)         = unusedIndicator;
	      newAlphaOccupation(ci,newIndex) = unusedIndicator;
	      newBetaOccupation(ci,newIndex)  = BetaOccupation(ci,o);
	    }
	  newIndex++;
	} else {
	  for(int ci=0; ci<Ndeterminants; ci++)
	    {
	      newAlphaOccupation(ci,o) = AlphaOccupation(ci,o);
	      newBetaOccupation(ci,o)  = BetaOccupation(ci,o);
	    } 
	}
    }

  OrbitalCoeffs   = newOrbitalCoeffs;
  AlphaOccupation = newAlphaOccupation;
  BetaOccupation  = newBetaOccupation;
  Norbitals       = finalNorbitals;
}

void QMCWavefunction::unlinkDeterminants()
{ 
  int finalNorbitals = 0;
  vector<int> unusedOrbitals;
  for(int o=0; o<Norbitals; o++)
    {
      bool used = false;
      for(int ci=0; ci<Ndeterminants; ci++)
	{
	  if(AlphaOccupation(ci,o) != unusedIndicator ||
	     BetaOccupation(ci,o) != unusedIndicator)
	    {
	      finalNorbitals++;
	      used = true;
	    }	    
	}
      
      if(!used)
	unusedOrbitals.push_back(o);
    }

  //int finalNorbitals = Nelectrons * Ndeterminants;
  finalNorbitals += unusedOrbitals.size();

  //The orbitals might already be fully decoupled.
  if(finalNorbitals == Norbitals)
    return;

  Array2D<qmcfloat> newOrbitalCoeffs(finalNorbitals,Nbasisfunc);
  Array2D<int>      newAlphaOccupation(Ndeterminants, finalNorbitals);
  Array2D<int>      newBetaOccupation(Ndeterminants, finalNorbitals);

  newAlphaOccupation = unusedIndicator;
  newBetaOccupation  = unusedIndicator;
  
  int newIndex = 0;
  for(int o=0; o<Norbitals; o++)
    for(int ci=0; ci<Ndeterminants; ci++)
      if(AlphaOccupation(ci,o) != unusedIndicator ||
	 BetaOccupation(ci,o) != unusedIndicator)
	{
	  newAlphaOccupation(ci,newIndex) = AlphaOccupation(ci,o);
	  newBetaOccupation(ci,newIndex) = BetaOccupation(ci,o);
	  newOrbitalCoeffs.setRows(newIndex,o,1,OrbitalCoeffs);
	  newIndex++;
	}

  if(finalNorbitals - newIndex != (int)(unusedOrbitals.size()))
    {
      clog << "Error: orbital counting went bad.\n";
      clog << "   finalNorbitals = " << finalNorbitals << endl;
      clog << "         newIndex = " << newIndex << endl;
      clog << "   unusedOrbitals = " << unusedOrbitals.size() << endl;
      exit(0);
    }
  for(unsigned int o=0; o<unusedOrbitals.size(); o++)
    {
      newOrbitalCoeffs.setRows(newIndex,o,1,OrbitalCoeffs);
      newIndex++;
    }

  OrbitalCoeffs   = newOrbitalCoeffs;
  AlphaOccupation = newAlphaOccupation;
  BetaOccupation  = newBetaOccupation;
  Norbitals       = finalNorbitals;
}

QMCWavefunction QMCWavefunction::operator=( const QMCWavefunction & rhs )
{
  Norbitals         = rhs.Norbitals;
  Nbasisfunc        = rhs.Nbasisfunc;
  Nalpha            = rhs.Nalpha;
  Nbeta             = rhs.Nbeta;
  Nelectrons        = rhs.Nelectrons;
  Ndeterminants     = rhs.Ndeterminants;
  trialFunctionType = rhs.trialFunctionType;
  OrbitalCoeffs     = rhs.OrbitalCoeffs;
  CI_coeffs         = rhs.CI_coeffs;
  CI_constraints    = rhs.CI_constraints;
  CI_const_coeffs   = rhs.CI_const_coeffs;
  numCIIndependent  = rhs.numCIIndependent;
  AlphaOccupation   = rhs.AlphaOccupation;
  BetaOccupation    = rhs.BetaOccupation;
  return *this;
}

istream& operator >>(istream &strm,QMCWavefunction &rhs)
{
  string temp_string;

  if (rhs.trialFunctionType == "harmonicoscillator")
    {
      rhs.Nalpha = abs(rhs.Ncharge);
      rhs.Nbeta = 0;
      rhs.Nelectrons = rhs.Nalpha + rhs.Nbeta;
      return strm;
    }

  int finalNorbitals = rhs.Norbitals;
  if (rhs.trialFunctionType == "restricted")
    {
      rhs.OrbitalCoeffs.allocate(rhs.Norbitals,
				 rhs.Nbasisfunc);

      for(int i=0; i<rhs.Norbitals; i++)
	for(int j=0; j<rhs.Nbasisfunc; j++)
	  strm >> rhs.OrbitalCoeffs(i,j);
    }
  
  else if (rhs.trialFunctionType == "unrestricted")
    { 
      finalNorbitals *= 2;
      rhs.OrbitalCoeffs.allocate(finalNorbitals,
				 rhs.Nbasisfunc);

      //skipping the words "Alpha Molecular Orbitals"
      strm >> temp_string;
      strm >> temp_string;
      strm >> temp_string;
      
      for(int i=0; i<rhs.Norbitals; i++)
	for(int j=0; j<rhs.Nbasisfunc; j++)
	  strm >> rhs.OrbitalCoeffs(i,j);
      
      //skipping the words "Beta Molecular Orbitals"      
      strm >> temp_string;
      strm >> temp_string;
      strm >> temp_string;

      for(int i=rhs.Norbitals; i < finalNorbitals; i++)
	for(int j=0; j<rhs.Nbasisfunc; j++)
	  strm >> rhs.OrbitalCoeffs(i,j);      
    }
  
  rhs.AlphaOccupation.allocate(rhs.Ndeterminants,finalNorbitals);
  rhs.AlphaOccupation = 0;
  rhs.BetaOccupation.allocate(rhs.Ndeterminants,finalNorbitals);
  rhs.BetaOccupation  = 0;

  //skipping the words "Alpha Occupation"
  strm >> temp_string; 
  strm >> temp_string;
  
  for (int i=0; i<rhs.Ndeterminants; i++)
    for(int j=0; j<rhs.Norbitals; j++)
      strm >> rhs.AlphaOccupation(i,j);
  
  //skipping the words "Beta Occupation"
  strm >> temp_string;
  strm >> temp_string;
  
  for (int i=0; i<rhs.Ndeterminants; i++)
    {
      int start = finalNorbitals - rhs.Norbitals;
      for(int j=start; j<finalNorbitals; j++)
	strm >> rhs.BetaOccupation(i,j);
    }
  
  //skipping the words "CI Coeffs"
  strm >> temp_string;
  strm >> temp_string;

  rhs.CI_coeffs.allocate(rhs.Ndeterminants);
  rhs.CI_constraints.allocate(rhs.Ndeterminants);
  rhs.CI_const_coeffs.allocate(rhs.Ndeterminants);
  rhs.numCIIndependent = rhs.Ndeterminants;
  rhs.CI_constraints = -1;
  rhs.CI_const_coeffs = 1.0;
  for(int i=0; i<rhs.Ndeterminants; i++)
    {
      string temp;
      strm >> temp;

      if(temp.find('c') != string::npos)
	{
	  int c;
	  strm >> c;
	  double d;
	  strm >> d;
	  // This coefficient is constrained
	  // when optimizing
	  if(rhs.CI_constraints(c) == -1){
	    rhs.CI_constraints(i) = c;
	  } else {
	    clog << "Warning: constraint for determinant "
		 << i << " depends on constrained determinant "
		 << c << " so we'll set it to the original determinant "
		 << rhs.CI_constraints(c) << "." << endl;
	    rhs.CI_constraints(i) = rhs.CI_constraints(c);
	  }

	  rhs.CI_const_coeffs(i) = d;
	  rhs.numCIIndependent--;
	} else {
	  rhs.CI_coeffs(i) = atof(temp.c_str());
	}
    }

  for(int i=0; i<rhs.Ndeterminants; i++)
    if(rhs.CI_constraints(i) != -1)
	rhs.CI_coeffs(i) = rhs.CI_const_coeffs(i) * rhs.CI_coeffs(rhs.CI_constraints(i));

  if(rhs.Ncharge != 0)
    {
      cerr << "Error: molecular charges have not been included in QMcBeaver\n";
    }

  rhs.Nalpha = 0;
  rhs.Nbeta = 0;
  for(int i=0; i<rhs.Norbitals; i++)
    {
      rhs.Nalpha += rhs.AlphaOccupation(0,i);
      rhs.Nbeta += rhs.BetaOccupation(0,i);
    } 
  
  rhs.Nelectrons = rhs.Nalpha + rhs.Nbeta;

  //Sanity check
  for(int ci=1; ci<rhs.Ndeterminants; ci++)
    {
      int numE = 0;
      for(int i=0; i<rhs.Norbitals; i++)
	{
	  numE += rhs.AlphaOccupation(ci,i);
	  numE += rhs.BetaOccupation(ci,i);
	}     
      if(numE != rhs.Nelectrons)
	{
	  clog << "Error: determinant " << (ci+1) << " has " << numE << " electrons.\n";
	  clog << "       It should have " << rhs.Nelectrons << " to be consistent with\n"
	       << "       the first determinant.\n";
	  exit(0);
	}
    }

  rhs.Norbitals = finalNorbitals;
  rhs.trialFunctionType = "restricted";

  strm >> temp_string;
  if(temp_string != "&")
    {
      clog << "ERROR: there was some error in reading the wavefunction. WF = " << endl;
      clog << rhs << endl;
      exit(0);
    }

  return strm;
}

void QMCWavefunction::read(int charge, int numberOrbitals, int numberBasisFunctions,
                           int numberDeterminants, string functionType, string runfile)
{
  Norbitals  = numberOrbitals;
  Nbasisfunc = numberBasisFunctions;
  Ndeterminants = numberDeterminants;
  Ncharge = charge;
  trialFunctionType = functionType;
  
  ifstream input_file(runfile.c_str());
  
  if(!input_file)
  {
    cerr << "ERROR: Could not open " << runfile << "!" << endl;
    exit(1);
  }
  
  string temp_string;
  input_file >> temp_string;
  while((temp_string != "&wavefunction") && (input_file.eof() != 1))
  {
    input_file >> temp_string;
  }
  input_file >> *this;

  if(getNumberElectrons(true) + getNumberElectrons(false) != getNumberElectrons()){
    cerr << "Error: We are expecting number alpha electrons " << getNumberElectrons(true) <<
      " + number beta electrons " << getNumberElectrons(false) << " to add up to the total " <<
      getNumberElectrons() << endl;
  }

  //rhs.scaleCoeffs(4.0);
  
  if(globalInput.flags.optimize_Psi == 1 &&
     globalInput.flags.optimize_Orbitals == 1)
    {
      if(globalInput.flags.link_Orbital_parameters == 0)
	unlinkOrbitals();
      
      if(globalInput.flags.link_Determinant_parameters == 0)
	unlinkDeterminants();

      if(globalInput.flags.Norbitals != Norbitals)
	{
	  clog << "Notice: the number of orbitals has increased from "
	       << globalInput.flags.Norbitals << " to "
	       << Norbitals << endl;
	  globalInput.flags.Norbitals = Norbitals;
	}
      sortOccupations(true);
      
      int numOR = getNumberOrbitals();
      int numORA= getNumberActiveOrbitals();
      int numBF = getNumberBasisFunctions();
      int numCI = getNumberDeterminants();      

      numORIndependent = numORA*numBF;
      OR_constraints.allocate(numORA*numBF);

      //-2 means constrained to 0, -1 means free, i means constrained to free parameter i
      OR_constraints = -1;

      int ai = 0;
      for(int o=0; o<numOR; o++)
	{
	  bool orbitalUsed = false;
	  for(int ci=0; ci<numCI; ci++)
	    {
	      if(AlphaOccupation(ci,o) != unusedIndicator ||
		 BetaOccupation(ci,o) != unusedIndicator)
		{
		  orbitalUsed = true;
		  break;
		}
	    }
	  
	  if(!orbitalUsed) continue;
	  int orbStart = ai;
	  
	  for(int bf=0; bf<numBF; bf++)
	    {
	      if(globalInput.flags.constrain_Orbital_zeros == 1 &&
		 fabs(OrbitalCoeffs(o,bf)) < 1e-8 )
		{
		  numORIndependent--;
		  OR_constraints(ai) = -2;
		} else if(globalInput.flags.constrain_Orbital_same == 1) {

		  for(int bfp=0; bfp<bf; bfp++)
		    if(fabs(OrbitalCoeffs(o,bf) - OrbitalCoeffs(o,bfp)) < 1e-8)
		      {
			OR_constraints(ai) = orbStart+bfp;
			numORIndependent--;
			break;
		      }
		}	      
	      ai++;
	    }
	}  

      if(numORIndependent != OR_constraints.dim1())
	cout << "Notice: wavefunction went from " << OR_constraints.dim1()
	     << " optimizable orbital parameters to " << numORIndependent << endl;
      
    }
  else
    {
      sortOccupations(false);
    }

  makeCoefficients(AlphaOccupation,
		   OrbitalCoeffs,
		   AlphaOrbitals,
		   AlphaCoeffs,
		   AlphaSwaps);
  makeCoefficients(BetaOccupation,
		   OrbitalCoeffs,
		   BetaOrbitals,
		   BetaCoeffs,
		   BetaSwaps);
}

ostream& operator <<(ostream& strm, QMCWavefunction& rhs)
{
  strm << "&wavefunction" << endl;

  int wf_width = 25;
  int wf_precision = 12;
  //strm.setf(ios::scientific);

  if (rhs.trialFunctionType == "restricted")
    {
      for(int i=0; i<rhs.Norbitals; i++)
	{
	  for(int j=0; j<rhs.Nbasisfunc; j++)
	    {
	      if( j%3 == 0 && j > 0 )
		strm << endl;
	      strm.unsetf(ios::scientific);
	      strm.unsetf(ios::fixed);
	      //strm.setf(ios::fixed);
	      strm.precision(wf_precision);
	      if(rhs.OrbitalCoeffs(i,j) < 0.0) strm << " -";
	      else                             strm << "  ";
	      strm << setw(wf_width) << left << fabs(rhs.OrbitalCoeffs(i,j));
	    }
	  strm << endl << endl;
	}
    }
  else if (rhs.trialFunctionType == "unrestricted")
    {
      clog << "Warning: why is the wavefunction unrestricted?" << endl;
    }
  else if (rhs.trialFunctionType == "harmonicoscillator")
    {
      return strm;
    }
  
  strm << "Alpha Occupation" << endl;
  for (int i=0; i<rhs.Ndeterminants; i++)
    {
      for (int j=0; j<rhs.Norbitals; j++)
	{
	  //since the output will only be 1 or 0
	  if(rhs.AlphaOccupation(i,j) == rhs.unusedIndicator)
	    strm << setw(4) << 0;
	  else
	    strm << setw(4) << 1;
	}
      strm << endl;
    }
  
  strm << endl << "Beta Occupation" << endl;
  for (int i=0; i<rhs.Ndeterminants; i++)
    {
      for (int j=0; j<rhs.Norbitals; j++)
	{
	  strm.width(2);
	  if(rhs.BetaOccupation(i,j) == rhs.unusedIndicator)
	    strm << setw(4) << 0;
	  else
	    strm << setw(4) << 1;
	}
      strm << endl;
    }
  
  strm << endl << "CI Coeffs" << endl;
  for (int i=0; i<rhs.Ndeterminants; i++)
    {
      if(rhs.CI_constraints(i) == -1)
	{
	  strm.precision(wf_precision);
	  if(rhs.CI_coeffs(i) < 0.0) strm << " -";
	  else                       strm << "  ";
	  strm << setw(wf_width) << left << fabs(rhs.CI_coeffs(i)) << endl;
	} else {
	  strm << " c " << rhs.CI_constraints(i) << " ";
	  strm.precision(wf_precision);
	  strm << setw(wf_width) << left << rhs.CI_const_coeffs(i) << endl;
	}
    }
  strm << endl;
  
  strm << "&" << endl << endl;
  return strm;
}

int QMCWavefunction::getNumberCIParameters()
{
  if(globalInput.flags.optimize_CI == 0)
    return 0;
  /*
    Since the wavefunction is invariant to normalization,
    there are only N-1 independent CI coefficients.
  */
  return numCIIndependent;
}

void QMCWavefunction::getCIParameters(Array1D<double> & params, int shift)
{
  if(getNumberCIParameters() <= 0)
    return;

  int ci = 0;
  for(int i=0; i<getNumberDeterminants(); i++)
    {
      if(CI_constraints(i) == -1)
	{
	  params(ci + shift) = CI_coeffs(i);
	  ci++;
	}
    }
}

void QMCWavefunction::setCIParameters(Array1D<double> & params, int shift)
{
  if(getNumberCIParameters() <= 0)
    return;
  
  int ci = 0;
  for(int i=0; i<getNumberDeterminants(); i++)
    {
      if(CI_constraints(i) == -1)
	{
	  CI_coeffs(i) = params(ci + shift);
	  ci++;
	}
    }
  for(int i=0; i<getNumberDeterminants(); i++)
    {
      if(CI_constraints(i) != -1)
	CI_coeffs(i) = CI_const_coeffs(i) * CI_coeffs(CI_constraints(i));
    }
}

double QMCWavefunction::getCINorm()
{
  return CI_coeffs.mag();
}

void QMCWavefunction::normalizeCI()
{
  double norm = CI_coeffs.mag();
  for(int ci=0; ci<Ndeterminants; ci++)
    CI_coeffs(ci) = CI_coeffs(ci)/norm;
}

int QMCWavefunction::getNumberORParameters()
{
  if(globalInput.flags.optimize_Orbitals == 0)
    return 0;

  return numORIndependent;
}

void QMCWavefunction::getORParameters(Array1D<double> & params, int shift)
{
  if(globalInput.flags.optimize_Orbitals == 0)
    return;

  int numOR = getNumberOrbitals();
  int numBF = getNumberBasisFunctions();
  int numCI = getNumberDeterminants();
  
  int ai = shift;
  int ind = 0;
  for(int o=0; o<numOR; o++)
    {
      bool orbitalUsed = false;
      for(int ci=0; ci<numCI; ci++)
	{
	  if(AlphaOccupation(ci,o) != unusedIndicator ||
	     BetaOccupation(ci,o) != unusedIndicator)
	    {
	      orbitalUsed = true;
	      break;
	    }
	}

      if(!orbitalUsed) continue;

      for(int bf=0; bf<numBF; bf++)
	{
	  if(OR_constraints(ind) == -1)
	    {
	      params(ai) = OrbitalCoeffs(o,bf);
	      ai++;
	    }
	  ind++;
	}
    }  
}

void QMCWavefunction::setORParameters(Array1D<double> & params, int shift)
{
  if(globalInput.flags.optimize_Orbitals == 0)
    return;

  int numOR = getNumberOrbitals();
  int numBF = getNumberBasisFunctions();
  int numCI = getNumberDeterminants();
  
  int ai = shift;
  Array1D<double> newParams(OR_constraints.dim1());
  int ind = 0;
  for(int o=0; o<numOR; o++)
    {
      bool orbitalUsed = false;
      for(int ci=0; ci<numCI; ci++)
	{
	  if(AlphaOccupation(ci,o) != unusedIndicator ||
	     BetaOccupation(ci,o) != unusedIndicator)
	    {
	      orbitalUsed = true;
	      break;
	    }
	}

      if(!orbitalUsed) continue;

      for(int bf=0; bf<numBF; bf++)
	{
	  int c = OR_constraints(ind);
	  if(c == -1)
	    {
	      OrbitalCoeffs(o,bf) = params(ai);
	      ai++;
	    } else if(c == -2) {
	      OrbitalCoeffs(o,bf) = 0.0;
	    } else {
	      OrbitalCoeffs(o,bf) = newParams(c);
	    }
	  newParams(ind) = OrbitalCoeffs(o,bf);
	  ind++;
	}
    }  

  makeCoefficients(AlphaOccupation,
		   OrbitalCoeffs,
		   AlphaOrbitals,
		   AlphaCoeffs,
		   AlphaSwaps);
  makeCoefficients(BetaOccupation,
		   OrbitalCoeffs,
		   BetaOrbitals,
		   BetaCoeffs,
		   BetaSwaps);
}

Array2D<qmcfloat> * QMCWavefunction::getCoeff(bool isAlpha)
{
  if(isAlpha)
    {
      return & AlphaCoeffs;
    } else {
      return & BetaCoeffs;
    }
}

Array2D<int> * QMCWavefunction::getOccupations(bool isAlpha)
{
  if(isAlpha)
    {
      return & AlphaOrbitals;
    } else {
      return & BetaOrbitals;
    }
}

Array2D<int> * QMCWavefunction::getDeterminantSwaps(bool isAlpha)
{
  if(isAlpha)
    {
      return & AlphaSwaps;
    } else {
      return & BetaSwaps;
    }
}

void QMCWavefunction::makeCoefficients(Array2D<int> & TypeOccupations,
				       Array2D<qmcfloat> & AllCoeffs,
				       Array2D<int> & TypeIndices,
				       Array2D<qmcfloat> & TypeCoeffs,
				       Array2D<int> & TypeSwaps)
{
  int numOR = getNumberOrbitals(TypeOccupations);
  int numAC = getNumberActiveOrbitals(TypeOccupations);
  int numCI = Ndeterminants;

  TypeCoeffs.allocate(numAC,Nbasisfunc);
  TypeIndices.allocate(numCI,numOR);

  int idxAct = 0;
  Array1D<int> count(numCI);
  count = 0;
  for(int o=0; o<AllCoeffs.dim1(); o++)
    {
      bool used = false;
      for(int ci=0; ci<numCI; ci++)
	if(TypeOccupations(ci,o) != unusedIndicator)
	  {
	    used = true;
	    TypeIndices(ci,count(ci)) = idxAct;
	    count(ci) = count(ci) + 1;
	  }

      if(used)
	{
	  TypeCoeffs.setRows(idxAct,o,1,OrbitalCoeffs);
	  idxAct++;
	}
    }

  bool goodCI = true;
  for(int ci=0; ci<count.dim1(); ci++)
    if(count(ci) != numOR)
      goodCI = false;

  if(idxAct != numAC || !goodCI)
    {
      clog << "Error: we failed to make the coefficient matrices\n"
	   << " numCI = " << numCI << endl
	   << " numOR = " << numOR << endl
	   << " numAC = " << numAC << endl
	   << "idxAct = " << idxAct << endl;
      if(!goodCI)
	clog << " CI count = \n" << count << endl;
      exit(0);
    }

  /*
    We want to figure out a procedure for updating determinant
    ci by swapping orbitals out of determinant ci-1. A problem
    would occur if the same orbital occurs in both determinants,
    but in different positions.

    To avoid this, we need to swap orbitals in an order such that
    no intermediate stage of updating has the same orbital twice.
  */
  TypeSwaps.allocate(numCI,numOR);
  Array2D<int> IndexUpdater = TypeIndices;
  TypeSwaps = -1;
  for(int ci=1; ci<numCI; ci++)
    {
      int swap = 0;
      bool badSwap = false;

      /*
	We scan through all the orbitals, looking for orbitals we can update. If
	we ever ran into a bad update, then we ignore it for the moment.

	Once we've completed a scan, then hopefully the bad dependencies will be
	gone.
      */
      do {
	badSwap = false;
	for(int o=0; o<numOR; o++)
	  {
	    // Check to see if the orbitals already match, and in the same
	    // position; no update necessary
	    if(IndexUpdater(ci,o) == IndexUpdater(ci-1,o))
	      continue;

	    // Figure out if swapping this orbital would introduce a
	    // linear dependency in the determinant
	    bool goodSwap = true;
	    for(int o2=0; o2<numOR; o2++)
	      if(o != o2 && IndexUpdater(ci,o) == IndexUpdater(ci-1,o2))
		{
		  badSwap = true;
		  goodSwap = false;
		  break;		
		}
	    
	    //We're allowed to update this orbital
	    if(goodSwap)
	      {
		IndexUpdater(ci-1,o) = IndexUpdater(ci,o);
		TypeSwaps(ci,swap) = o;
		swap++;
	      }    
	  }
      } while(badSwap);
    }
}
