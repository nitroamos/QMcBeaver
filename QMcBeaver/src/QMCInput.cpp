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

 Molecule.initialize(flags.Natoms, flags.pseudo_gridLevel);
 Molecule.readGeometry(inputfile);
 flags.use_pseudopotential = Molecule.readPseudoPotential(inputfile);

 BF.initialize(&flags, &Molecule);
 BF.read(inputfile);

 WF.read(flags.charge, flags.Norbitals, flags.Nbasisfunc, flags.Ndeterminants, 
	 flags.trial_function_type, inputfile);

 if(Molecule.getNuclearCharge() != WF.getNumberElectrons() + flags.charge)
   {
     clog << "Error: incorrect number of electrons." << endl;
     clog << "    Nuclear charge = " << Molecule.getNuclearCharge() << endl;
     clog << "  Number electrons = " << WF.getNumberElectrons() << endl;
     clog << " Electronic charge = " << flags.charge << endl;
     exit(0);
   }

 JP.read(Molecule.NucleiTypes,flags.link_Jastrow_parameters,flags.replace_electron_nucleus_cusps,
	 WF.getNumberElectrons(true),WF.getNumberElectrons(false),inputfile);

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

void QMCInput::printArray(ostream & strm,
		string name, int num,
		Array1D<double> & array, int & start,
		int margin, int width, int prec, int numPerRow)
{
  if(num <= 0) return;

  strm << setw(margin) << name;

  for(int ai=0; ai<num; ai++)
    {
      if(ai%numPerRow == 0 && ai != 0)
	strm << endl << setw(margin) << " ";
      
      strm.width(width);
      strm.precision(prec);
      strm << array(ai + start);
    }
  strm << endl;
  start += num;
}

void QMCInput::printAIParameters(ostream & strm, 
				 string name, 
				 int margin,
				 Array1D<double> & array,
				 bool forcePrintOrbitals)
{
  if(array.dim1() <= 0) return;
  /*
    This function assumes a particular ordering of the parameters
    in the array.

    There are obviously quite a few Orbital parameters, so we probably
    don't want to print all of them.
  */
  int width = 20;
  int prec  = 12;
  margin += 7;

  int shift = 0;
  printArray(strm, name + "  UD = ", JP.getNumberEupEdnParameters(),
	     array, shift, margin, width, prec, 5);
  printArray(strm, name + "  UU = ", JP.getNumberEupEupParameters(),
	     array, shift, margin, width, prec, 5);
  printArray(strm, name + "  DD = ", JP.getNumberEdnEdnParameters(),
	     array, shift, margin, width, prec, 5);

  printArray(strm, name + "  NE = ", JP.getNumberNEParameters(),
	     array, shift, margin, width, prec, 5);

  printArray(strm, name + " NUD = ", JP.getNumberNEupEdnParameters(),
	     array, shift, margin, width, prec, 5);
  printArray(strm, name + " NUU = ", JP.getNumberNEupEupParameters(),
	     array, shift, margin, width, prec, 5);
  printArray(strm, name + " NDD = ", JP.getNumberNEdnEdnParameters(),
	     array, shift, margin, width, prec, 5);

  printArray(strm, name + "  CI = ", WF.getNumberCIParameters(),
	     array, shift, margin, width, prec, 5);

  int numOR = WF.getNumberORParameters();
  if(numOR > 0)
    {
      if(forcePrintOrbitals || numOR < 10)
	{
	  int numBF = WF.getNumberBasisFunctions();
	  strm << setw(margin) << name << "  OR =";
	  for(int ai=0; ai<numOR; ai++)
	    {
	      if(ai%5 == 0 && ai%numBF != 0)
		strm << endl;
	      if(ai%numBF == 0)
		strm << endl << endl;

	      strm.width(width);
	      strm.precision(prec);
	      strm << array(ai + shift);
	    }
	  strm << endl << endl << endl;
	} else {
	  strm << setw(margin) << name << " OR = ";
	  strm << "<" << numOR << " parameters not printed>";
	  strm << endl;
	}
    }
}

void QMCInput::printAISummary()
{
  int width = 10;
  int numAI = globalInput.getNumberAIParameters();

  clog << "There are currently " << numAI << " optimizable parameters:" << endl;  
  
  int num = globalInput.JP.getNumberEEParameters();
  if(num != 0) clog << setw(width) << num << "  Electron-Electron Jastrow parameters" << endl;
  
  num = globalInput.JP.getNumberNEParameters();
  if(num != 0) clog << setw(width) << num << "  Nuclear-Electron Jastrow parameters" << endl;
  
  num = globalInput.JP.getNumberNEupEdnParameters();
  if(num != 0) clog << setw(width) << num << "  NEupEdn Jastrow parameters" << endl;
  
  num = globalInput.JP.getNumberNEupEupParameters();
  if(num != 0) clog << setw(width) << num << "  NEupEup Jastrow parameters" << endl;
  
  num = globalInput.JP.getNumberNEdnEdnParameters();
  if(num != 0) clog << setw(width) << num << "  NEdnEdn Jastrow parameters" << endl;
  
  num = globalInput.WF.getNumberCIParameters();
  if(num != 0) clog << setw(width) << num << " CI parameters" << endl;
  
  num = globalInput.WF.getNumberORParameters();
  if(num != 0) clog << setw(width) << num << " orbital parameters" << endl;
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
