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

#include "QMCInput.h"
#include <iomanip>

QMCInput::QMCInput()
{
}

void QMCInput::setMPIParameters(int my_rank, int nprocs)
{
  flags.my_rank = my_rank;
  flags.nprocs  = nprocs;
}

void QMCInput::read(string inputfile)
{
 flags.read_flags(inputfile);

 Molecule.initialize(flags.Natoms);
 Molecule.read(inputfile);

 BF.initialize(&flags, &Molecule);
 BF.read(inputfile);

 WF.read(flags.charge, flags.Norbitals, flags.Nbasisfunc, flags.Ndeterminants, 
	 flags.trial_function_type, inputfile);

 JP.read(Molecule.NucleiTypes,flags.link_Jastrow_parameters,flags.replace_electron_nucleus_cusps,
	 WF.getNumberAlphaElectrons(),WF.getNumberBetaElectrons(),inputfile);

 outputer = QMCConfigIO(WF.getNumberElectrons());
}

void QMCInput::openConfigFile()
{
  if(flags.print_configs == 1)
    outputer.open(flags.config_file_name,true);
}

Array1D<double> QMCInput::getParameters()
{
  Array1D<double> Guess_Jastrow_parameters = JP.getParameters();
  int numJP = JP.getParameters().dim1();

  //we never optimize the first CI coefficient
  int numCI = WF.CI_coeffs.dim1() - 1;

  Array1D<double> Guess_parameters(numJP + numCI);

  for(int i=0; i<numJP; i++)
    Guess_parameters(i) = Guess_Jastrow_parameters(i);
  for(int i=0; i<numCI; i++)
    Guess_parameters(i+numJP) = WF.CI_coeffs(i+1);
  return Guess_parameters;
}

void QMCInput::setParameterVector(Array1D<double> & params)
{
  Array1D<double> Guess_Jastrow_parameters = JP.getParameters();
  int numJP = JP.getParameters().dim1();

  //we never optimize the first CI coefficient
  int numCI = WF.CI_coeffs.dim1() - 1;

  for(int i=0; i<numJP; i++)
    Guess_Jastrow_parameters(i) = params(i);
  for(int i=0; i<numCI; i++)
    WF.CI_coeffs(i+1) = params(i+numJP);

  JP.setParameterVector(Guess_Jastrow_parameters);
}

void QMCInput::printAIArray(ostream & strm, string name, int width, Array1D<double> & array)
{
  int dim = array.dim1();
  int numEE = JP.getNumberEEParameters();
  int numNE = JP.getNumberNEParameters();

  for(int ai=0; ai<dim; ai++)
    {
      if(ai == 0)
	{
	  strm << setw(width) << name << " EE = ";
	} else if(ai == numEE && numNE > 0)
	  {
	    strm << endl << setw(width) << name << " NE = ";
	  } else if(ai == (numEE + numNE))
	    {
	      strm << endl << setw(width) << name << " CI = ";
	    }
      strm.width(20);
      strm.precision(12);
      strm << array(ai);
    }
  strm << endl;
}

ostream& operator<<(ostream & strm, QMCInput & Input)
{
  strm << Input.flags;
  strm << Input.Molecule;
  strm << Input.BF;
  strm << Input.WF;
  strm << Input.JP << endl;

  return strm;
}
