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

int QMCInput::getNumberAIParameters()
{
  int numAI = 0;
  numAI += JP.getNumberJWParameters();
  numAI += WF.getNumberCIParameters();
  numAI += WF.getNumberORParameters();
  return numAI;
}

Array1D<double> QMCInput::getAIParameters()
{
  Array1D<double> params(getNumberAIParameters());

  int shift = 0;
  JP.getJWParameters(params,shift);

  shift += JP.getNumberJWParameters();
  WF.getCIParameters(params,shift);

  shift += WF.getNumberCIParameters();
  WF.getORParameters(params,shift);

  return params;
}

void QMCInput::setAIParameters(Array1D<double> & params)
{
  int shift = 0;
  JP.setJWParameters(params,shift);

  shift += JP.getNumberJWParameters();
  WF.setCIParameters(params,shift);

  shift += WF.getNumberCIParameters();
  WF.setORParameters(params,shift);
}

void QMCInput::printAIParameters(ostream & strm, 
				 string name, 
				 int margin,
				 Array1D<double> & array,
				 bool forcePrintOrbitals)
{
  /*
    This function assumes a particular ordering of the parameters
    in the array.

    There are obviously quite a few Orbital parameters, so we probably
    don't want to print all of them.
  */
  int numAI = getNumberAIParameters();

  int width = 20;
  int prec  = 12;

  int shift = 0;
  int numEE = JP.getNumberEEParameters();
  if(numEE > 0)
    {
      strm << setw(margin) << name << " EE = ";
      for(int ai=0; ai<numEE; ai++)
	{
	  strm.width(width);
	  strm.precision(prec);
	  strm << array(ai + shift);
	}
      strm << endl;
    }

  shift += numEE;
  int numNE = JP.getNumberNEParameters();
  if(numNE > 0)
    {
      strm << setw(margin) << name << " NE = ";
      for(int ai=0; ai<numNE; ai++)
	{
	  strm.width(width);
	  strm.precision(prec);
	  strm << array(ai + shift);
	}
      strm << endl;
    }

  shift += numNE;
  int numCI = WF.getNumberCIParameters();
  if(numCI > 0)
    {
      strm << setw(margin) << name << " CI = ";
      for(int ai=0; ai<numCI; ai++)
	{
	  strm.width(width);
	  strm.precision(prec);
	  strm << array(ai + shift);
	}
      strm << endl;
    }

  shift += numCI;
  int numOR = WF.getNumberORParameters();
  if(numOR > 0)
    {
      if(forcePrintOrbitals || numOR < 10)
	{
	  int numBF = WF.getNumberBasisFunctions();
	  strm << setw(margin) << name << " OR =\n";
	  for(int ai=0; ai<numOR; ai++)
	    {
	      strm.width(width);
	      strm.precision(prec);
	      strm << array(ai + shift);
	      if(ai%5 == 0 && ai%numBF != 0)
		strm << endl;
	      if(ai%numBF == 0)
		strm << endl << endl;
	    }
	  strm << endl;
	} else {
	  strm << setw(margin) << name << " OR = ";
	  strm << "<" << numOR << " parameters not printed>";
	  strm << endl;
	}
    }
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
