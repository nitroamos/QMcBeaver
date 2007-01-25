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
 if( flags.optimize_Psi || flags.print_configs == 1 )
   outputer.open(flags.config_file_name,true);
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
