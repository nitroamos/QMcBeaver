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

int QMCWavefunction::getNumberActiveAlphaOrbitals()
{
  int NactiveorbsA = 0;
  for(int o=0; o<Norbitals; o++)
    {
      bool Aused = false;
      for(int ci=0; ci<Ndeterminants; ci++)
	{
	  if(AlphaOccupation(ci,o) != unusedIndicator)
	    Aused = true;
	}
      if(Aused)          NactiveorbsA++;
    }
  return NactiveorbsA;
}

int QMCWavefunction::getNumberActiveBetaOrbitals()
{
  int NactiveorbsB = 0;
  for(int o=0; o<Norbitals; o++)
    {
      bool Bused = false;
      for(int ci=0; ci<Ndeterminants; ci++)
	{
	  if(BetaOccupation(ci,o) != unusedIndicator)
	    Bused = true;
	}
      if(Bused)          NactiveorbsB++;
    }
  return NactiveorbsB;
}

int QMCWavefunction::getNumberBasisFunctions()
{
  return Nbasisfunc;
}

int QMCWavefunction::getNumberAlphaElectrons()
{
  return Nalpha;
}

int QMCWavefunction::getNumberBetaElectrons()
{
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

  Array2D<double> newOrbitalCoeffs(finalNorbitals,Nbasisfunc);
  Array2D<int>  newAlphaOccupation(Ndeterminants, finalNorbitals);
  Array2D<int>   newBetaOccupation(Ndeterminants, finalNorbitals);

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

  Array2D<double> newOrbitalCoeffs(finalNorbitals,Nbasisfunc);
  Array2D<int>  newAlphaOccupation(Ndeterminants, finalNorbitals);
  Array2D<int>   newBetaOccupation(Ndeterminants, finalNorbitals);

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
  for (int i=0; i<rhs.Ndeterminants; i++)
    {
      strm >> rhs.CI_coeffs(i);
    }

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

  if(getNumberAlphaElectrons() + 
     getNumberBetaElectrons() != getNumberElectrons()){
    cerr << "Error: We are expecting number alpha electrons " << getNumberAlphaElectrons() <<
      " + number beta electrons " << getNumberBetaElectrons() << " to add up to the total " <<
      getNumberElectrons() << endl;
  }

  //rhs.scaleCoeffs(4.0);
  
  if(globalInput.flags.optimize_Psi == 1)
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
    }

  else
    sortOccupations(false);

  makeCoefficients();
  /*
  cout << "AlphaOccupation\n" << AlphaOccupation << endl;
  cout << "BetaOccupation\n" << BetaOccupation << endl;
  cout << "numT = " << Norbitals << endl;
  cout << "numO = " << getNumberActiveOrbitals() << endl;
  cout << "numA = " << getNumberActiveAlphaOrbitals() << endl;
  cout << "numB = " << getNumberActiveBetaOrbitals() << endl;
  exit(0);
  //*/
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
	  if(rhs.AlphaOccupation(i,j) == -1)
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
	  if(rhs.BetaOccupation(i,j) == -1)
	    strm << setw(4) << 0;
	  else
	    strm << setw(4) << 1;
	}
      strm << endl;
    }
  
  strm << endl << "CI Coeffs" << endl;
  for (int i=0; i<rhs.Ndeterminants; i++)
    {
      strm.precision(wf_precision);
      if(rhs.CI_coeffs(i) < 0.0) strm << " -";
      else                       strm << "  ";
      strm << setw(wf_width) << left << fabs(rhs.CI_coeffs(i)) << endl;
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
  return getNumberDeterminants();
}

void QMCWavefunction::getCIParameters(Array1D<double> & params, int shift)
{
  for(int i=0; i<getNumberCIParameters(); i++)
    params(i + shift) = CI_coeffs(i);
}

void QMCWavefunction::setCIParameters(Array1D<double> & params, int shift)
{
  for(int i=0; i<getNumberCIParameters(); i++)
    CI_coeffs(i) = params(i + shift);
}

double QMCWavefunction::getCINorm()
{
  double norm = 0;
  for(int ci=0; ci<Ndeterminants; ci++)
    norm += CI_coeffs(ci)*CI_coeffs(ci);
  return norm;
}

void QMCWavefunction::normalizeCI()
{
  double norm = getCINorm();
  for(int ci=0; ci<Ndeterminants; ci++)
    CI_coeffs(ci) = CI_coeffs(ci)/norm;
}

int QMCWavefunction::getNumberORParameters()
{
  if(globalInput.flags.optimize_Orbitals == 0)
    return 0;

  int numBF = getNumberBasisFunctions();
  int numOR = getNumberActiveOrbitals();
  return numOR * numBF;
}

void QMCWavefunction::getORParameters(Array1D<double> & params, int shift)
{
  if(globalInput.flags.optimize_Orbitals == 0)
    return;

  int numOR = getNumberOrbitals();
  int numBF = getNumberBasisFunctions();
  int numCI = getNumberDeterminants();
  
  int ai = shift;
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
	  params(ai) = OrbitalCoeffs(o,bf);
	  ai++;
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
	  OrbitalCoeffs(o,bf) = params(ai);
	  ai++;
	}
    }  

  makeCoefficients();
}

Array2D<qmcfloat> * QMCWavefunction::getCoeff(int ci, bool isAlpha)
{
  if(isAlpha)
    {
      return & AlphaCoeffs(ci);
    } else {
      return & BetaCoeffs(ci);
    }
}

void QMCWavefunction::makeCoefficients()
{
  AlphaCoeffs.allocate(Ndeterminants);
  BetaCoeffs.allocate(Ndeterminants);

  for(int ci=0; ci<Ndeterminants; ci++)
    {
      AlphaCoeffs(ci).allocate(Nalpha, Nbasisfunc);
      BetaCoeffs(ci).allocate(Nbeta, Nbasisfunc);

      int idxa = 0;
      int idxb = 0;
      for(int o=0; o<Norbitals; o++)
	{
	  if(AlphaOccupation(ci,o) != unusedIndicator)
	    {
	      AlphaCoeffs(ci).setRows(idxa,o,1,OrbitalCoeffs);
	      idxa++;
	    }
	  if(BetaOccupation(ci,o) != unusedIndicator)
	    {
	      BetaCoeffs(ci).setRows(idxb,o,1,OrbitalCoeffs);
	      idxb++;
	    }
	}

      if(idxa != Nalpha || idxb != Nbeta)
	{
	  clog << "Error: we failed to make the coefficient matrices\n"
	       << "        ci = " << ci << endl
	       << " Norbitals = " << Norbitals << endl
	       << "    Nalpha = " << Nalpha << endl
	       << "      idxa = " << idxa << endl
	       << "     Nbeta = " << Nbeta << endl
	       << "      idxb = " << idxb << endl;
	  exit(0);
	}
    }
}
