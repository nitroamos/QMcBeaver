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


/**************************************************************************
This SOFTWARE has been authored or contributed to by an employee or 
employees of the University of California, operator of the Los Alamos 
National Laboratory under Contract No. W-7405-ENG-36 with the U.S. 
Department of Energy.  The U.S. Government has rights to use, reproduce, 
and distribute this SOFTWARE.  Neither the Government nor the University 
makes any warranty, express or implied, or assumes any liability or 
responsibility for the use of this SOFTWARE.  If SOFTWARE is modified 
to produce derivative works, such modified SOFTWARE should be clearly 
marked, so as not to confuse it with the version available from LANL.   

Additionally, this program is free software; you can distribute it and/or 
modify it under the terms of the GNU General Public License. Accordingly, 
this program is  distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY;  without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A  PARTICULAR PURPOSE.  See the GNU General Public License 
for more details. 
**************************************************************************/


#include "QMCJastrowParameters.h"
#include "QMCInput.h"

void QMCJastrowParameters::operator=( const QMCJastrowParameters & rhs )
{
  NumberOfEEParameters   = rhs.NumberOfEEParameters;
  NumberOfNEParameters   = rhs.NumberOfNEParameters;
  NumberOfNEupParameters = rhs.NumberOfNEupParameters;
  EupNuclear = rhs.EupNuclear;
  EdnNuclear = rhs.EdnNuclear;
  EupEdn = rhs.EupEdn;
  EupEup = rhs.EupEup;
  EdnEdn = rhs.EdnEdn;
  EquivalentElectronUpDownParams = rhs.EquivalentElectronUpDownParams;
  NucleiTypes = rhs.NucleiTypes;
  NumberOfElectronsUp = rhs.NumberOfElectronsUp;
  NumberOfElectronsDown = rhs.NumberOfElectronsDown;
}

void QMCJastrowParameters::setJWParameters(Array1D<double> & params, int shift)
{
  if( EquivalentElectronUpDownParams )
    {
      int Counter = 0;
      int CurrentNumberOfParams = 0;
      Array1D<double> temp;

      // Set eup edn
      if(getNumberEEParameters() > 0)
	{
	  if( NumberOfElectronsUp > 0 && NumberOfElectronsDown > 0 )
	    {
	      CurrentNumberOfParams = EupEdn.getTotalNumberOfParameters();
	      
	      temp.allocate( CurrentNumberOfParams );
	      
	      for(int i=0; i<temp.dim1(); i++)
		{
		  temp(i) = params(Counter + shift);
		  Counter++;
		}
	      
	      EupEdn.setParameters( temp );
	    }
	  
	  // Set eup eup
	  
	  if( NumberOfElectronsUp > 1 )
	    {
	      CurrentNumberOfParams = EupEup.getTotalNumberOfParameters();
	      
	      temp.allocate(CurrentNumberOfParams);
	      
	      for(int i=0; i<temp.dim1(); i++)
		{
		  temp(i) = params(Counter + shift);
		  Counter++;
		}
	      
	      EupEup.setParameters( temp );
	    }
	}

      // Set eup nuc
      if(getNumberNEParameters() > 0)
	{
	  for(int i=0; i<EupNuclear.dim1(); i++)
	    {
	      CurrentNumberOfParams = 
		EupNuclear(i).getTotalNumberOfParameters();
	      
	      temp.allocate(CurrentNumberOfParams);
	      
	      for(int j=0; j<temp.dim1(); j++)
		{
		  temp(j) = params(Counter + shift);
		  Counter++;
		}
	      
	      EupNuclear(i).setParameters( temp );
	    }
	}

      // Now set the things equal that need to be
      
      if( NumberOfElectronsDown > 1 )
	{
	  EdnEdn = EupEup;
	  EdnEdn.setParticle1Type("Electron_down");
	  EdnEdn.setParticle2Type("Electron_down");
	}
      
      if( NumberOfElectronsDown > 0 )
	{
	  EdnNuclear = EupNuclear;
	  
	  for(int i=0; i<EdnNuclear.dim1(); i++)
	    {
	      EdnNuclear(i).setParticle1Type("Electron_down");
	    }
	}
    }
  else
    {
      int Counter = 0;
      int CurrentNumberOfParams = 0;
      Array1D<double> temp;

      // Set eup edn
      if(getNumberEEParameters() > 0)
	{
	  if( NumberOfElectronsUp > 0 && NumberOfElectronsDown > 0 )
	    {
	      CurrentNumberOfParams = EupEdn.getTotalNumberOfParameters();
	      
	      temp.allocate( CurrentNumberOfParams );
	      
	      for(int i=0; i<temp.dim1(); i++)
		{  
		  temp(i) = params(Counter + shift);
		  Counter++;
		}
	      
	      EupEdn.setParameters( temp );
	    }
	  
	  // Set eup eup
	  
	  if( NumberOfElectronsUp > 1 )
	    {
	      CurrentNumberOfParams = EupEup.getTotalNumberOfParameters();
	      
	      temp.allocate(CurrentNumberOfParams);
	      
	      for(int i=0; i<temp.dim1(); i++)
		{
		  temp(i) = params(Counter + shift);
		  Counter++;
		}
	      
	      EupEup.setParameters( temp );
	    }
	  
	  // Set edn edn
	  
	  if( NumberOfElectronsDown > 1 )
	    {
	      CurrentNumberOfParams = EdnEdn.getTotalNumberOfParameters();
	      
	      temp.allocate(CurrentNumberOfParams);
	      
	      for(int i=0; i<temp.dim1(); i++)
		{
		  temp(i) = params(Counter + shift);
		  Counter++;
		}
	      
	      EdnEdn.setParameters( temp );
	    }
	}

      // Set eup nuc
      if(getNumberNEParameters() > 0)
	{
	  for(int i=0; i<EupNuclear.dim1(); i++)
	    {
	      CurrentNumberOfParams = 
		EupNuclear(i).getTotalNumberOfParameters();
	      
	      temp.allocate(CurrentNumberOfParams);
	      
	      for(int j=0; j<temp.dim1(); j++)
		{
		  temp(j) = params(Counter + shift);
		  Counter++;
		}
	      
	      EupNuclear(i).setParameters( temp );
	    }
	  
	  // Set edn nuc
	  
	  if( NumberOfElectronsDown > 0 )
	    {
	      for(int i=0; i<EdnNuclear.dim1(); i++)
		{
		  CurrentNumberOfParams = 
		    EdnNuclear(i).getTotalNumberOfParameters();
		  
		  temp.allocate(CurrentNumberOfParams);
		  
		  for(int j=0; j<temp.dim1(); j++)
		    {
		      temp(j) = params(Counter + shift);
		      Counter++;
		    }
		  
		  EdnNuclear(i).setParameters( temp );
		}
	    }
	}
    }
}

int QMCJastrowParameters::getNumberJWParameters()
{
  return getNumberEEParameters() + getNumberNEParameters();
}

int QMCJastrowParameters::getNumberEEParameters()
{
  if(globalInput.flags.optimize_EE_Jastrows == 0)
    return 0;
  return NumberOfEEParameters;
}

int QMCJastrowParameters::getNumberNEParameters()
{
  if(globalInput.flags.optimize_EN_Jastrows == 0)
    return 0;
  return NumberOfNEParameters;
}

int QMCJastrowParameters::getNumberNEupParameters()
{
  return NumberOfNEupParameters;
}

Array1D<double> QMCJastrowParameters::getJWParameters()
{
  Array1D<double> params(getNumberJWParameters());
  getJWParameters(params,0);
  return params;
}

void QMCJastrowParameters::getJWParameters(Array1D<double> & params, int shift)
{ 
  if( EquivalentElectronUpDownParams )
    {
      int Counter = 0;
      Array1D<double> temp;

      // Put in eup edn params
      if(getNumberEEParameters() > 0)
	{
	  if( NumberOfElectronsUp > 0 && NumberOfElectronsDown > 0 )
	    {
	      temp = EupEdn.getParameters();
	      
	      for(int i=0; i<temp.dim1(); i++)
		{
		  params(Counter + shift) = temp(i);
		  Counter++;
		}
	    }
	  
	  // Put in eup eup params
	  
	  if( NumberOfElectronsUp > 1 )
	    {
	      temp = EupEup.getParameters();
	      
	      for(int i=0; i<temp.dim1(); i++)
		{
		  params(Counter + shift) = temp(i);
		  Counter++;
		}
	    }
	}

      // Put in eup nuclear params

      if(getNumberNEParameters() > 0)
	{
	  for(int i=0; i<EupNuclear.dim1(); i++)
	    {
	      temp = EupNuclear(i).getParameters();
	      
	      for(int j=0; j<temp.dim1(); j++)
		{
		  params(Counter + shift) = temp(j);
		  Counter++;
		}
	    }
	}
    }
  else
    {      
      int Counter = 0;
      Array1D<double> temp;
      
      // Put in eup edn params
      if(getNumberEEParameters() > 0)
	{
	  if( NumberOfElectronsUp > 0 && NumberOfElectronsDown > 0 )
	    {
	      temp = EupEdn.getParameters();
	      
	      for(int i=0; i<temp.dim1(); i++)
		{
		  params(Counter + shift) = temp(i);
		  Counter++;
		}
	    }

	  // Put in eup eup params
	  
	  if( NumberOfElectronsUp > 1 )
	    {
	      temp = EupEup.getParameters();
	      
	      for(int i=0; i<temp.dim1(); i++)
		{
		  params(Counter + shift) = temp(i);
		  Counter++;
		}
	    }
	  
	  // Put in edn edn params
	  
	  if( NumberOfElectronsDown > 1 )
	    {
	      temp = EdnEdn.getParameters();
	      
	      for(int i=0; i<temp.dim1(); i++)
		{
		  params(Counter + shift) = temp(i);
		  Counter++;
		}
	    }
	}
      
      // Put in eup nuc params
      if(getNumberNEParameters() > 0)
	{
	  for(int i=0; i<EupNuclear.dim1(); i++)
	    {
	      temp = EupNuclear(i).getParameters();
	      
	      for(int j=0; j<temp.dim1(); j++)
		{
		  params(Counter + shift) = temp(j);
		  Counter++;
		}
	    }
	  
	  // Put in edn nuc params
	  
	  if( NumberOfElectronsDown > 0 )
	    {
	      for(int i=0; i<EdnNuclear.dim1(); i++)
		{
		  temp = EdnNuclear(i).getParameters();
		  
		  for(int j=0; j<temp.dim1(); j++)
		    {
		      params(Counter + shift) = temp(j);
		      Counter++;
		    }
		}
	    }
	}
    }
}

Array1D<Complex> QMCJastrowParameters::getPoles()
{
  list<Complex> allPoles;

  if( NumberOfElectronsUp > 1 )
    {
      // up up jastrow

      Array1D<Complex> poles = EupEup.getPoles();
      
      for(int i=0; i<poles.dim1(); i++)
	{
	  allPoles.push_back(poles(i));
	}
    }

  if( NumberOfElectronsDown > 1 )
    {
      // down down jastrow

      Array1D<Complex> poles = EdnEdn.getPoles();
      
      for(int i=0; i<poles.dim1(); i++)
	{
	  allPoles.push_back(poles(i));
	}
    }

  if( NumberOfElectronsUp > 0 && NumberOfElectronsDown > 0 )
    {
      // up down jastrow

      Array1D<Complex> poles = EupEdn.getPoles();
      
      for(int i=0; i<poles.dim1(); i++)
	{
	  allPoles.push_back(poles(i));
	}
    }

  for(int i=0; i<EupNuclear.dim1(); i++)
    {
      // up nuclear jastrow

      Array1D<Complex> poles = EupNuclear(i).getPoles();
      
      for(int j=0; j<poles.dim1(); j++)
	{
	  allPoles.push_back(poles(j));
	}
    }

  for(int i=0; i<EdnNuclear.dim1(); i++)
    {
      // down nuclear jastrow

      Array1D<Complex> poles = EdnNuclear(i).getPoles();
      
      for(int j=0; j<poles.dim1(); j++)
	{
	  allPoles.push_back(poles(j));
	}
    }

  // pack up the results into a vector
  Array1D<Complex> results( (int)allPoles.size() );

  int index = 0;
  for(list<Complex>::iterator it=allPoles.begin(); it!=allPoles.end(); ++it)
    {
      results(index) = *it;
      index++;
    }

  return results;
}

double QMCJastrowParameters::calculate_penalty_function()
{
  Array1D<Complex> poles = getPoles();
  return calculate_penalty_function(poles);
}

// calculates a penalty function for getting singular parameters
double QMCJastrowParameters::calculate_penalty_function(
							Array1D<Complex> & poles)
{
  double penalty = 0.0;

  for(int i=0; i<poles.dim1(); i++)
    {
      // calculate the distance of the pole from the positive real axis
      double distance = 0.0;

      if( poles(i).real() > 0 )
	{
	  distance = fabs(poles(i).imaginary());
	}
      else
	{
	  distance = poles(i).abs();
	}

      if(distance <= 0){
	penalty = 1e100; 
	cerr << "Warning: distance from real axis is " << distance << "  in OFR, can\'t take log" << endl;
	cerr << "         poles(" << i << ") = " << poles(i) << endl;
	cerr << "         penalty will be set to " << penalty << endl;
	cerr.flush();
      } else {

	/*
	  Singularities exist when the poles are on the +'ve real
	  axis. So, a configuration's score will be improved (lowered)
	  the further it is from this region.
	*/
	penalty -= log( distance );
      }
    }
  
  return penalty;
}

void QMCJastrowParameters::read(Array1D<string> & nucleitypes, 
				bool linkparams, bool nucCuspReplacement, int nelup, int neldn, 
				string runfile)
{
  // Set the internal variable telling if the electrons of different spins
  // should use the same parameters
  EquivalentElectronUpDownParams = linkparams;

  // Save an array containing all of nuclei types
  NucleiTypes = nucleitypes;

  // Save the number of up and down electrons
  NumberOfElectronsUp = nelup;
  NumberOfElectronsDown = neldn;

  if( NumberOfElectronsDown > NumberOfElectronsUp )
    {
      cerr << "ERROR: In QMCJastrowParameters::read the number of down "
	   << "electrons is greater than the number of up electrons.  The "
	   << "opposite is assumed." << endl;
      exit(0);
    }

  // Open the file to read

  ifstream file( runfile.c_str() );

  if( !file )
    {
      cerr << "ERROR: Could not open " << runfile << "!" << endl;
      exit(0);
    }

  string temp;
  file >> temp;

  while( temp != "&Jastrow" )
    {
      file >> temp;

      if( file.eof() )
	{
	  cerr << "ERROR: No Jastrow secton found in " << runfile << "!" 
	       << endl;
	  exit(0);
	}
    }
  
  // Read in the QMCCorrelationFunctionParameters

  EupNuclear.allocate( NucleiTypes.dim1() );

  if( NumberOfElectronsDown > 0 )
    EdnNuclear.allocate( NucleiTypes.dim1() );

  // Determine the number of correlation functions

  int NumberOfCorrelationFunctions = NucleiTypes.dim1();
  
  if( NumberOfElectronsDown > 0 )
    NumberOfCorrelationFunctions += NucleiTypes.dim1();

  if( NumberOfElectronsUp >0 && NumberOfElectronsDown > 0 )
    NumberOfCorrelationFunctions++;

  if( NumberOfElectronsUp > 1 )
    NumberOfCorrelationFunctions++;

  if( NumberOfElectronsDown > 1 )
    NumberOfCorrelationFunctions++;

  // Load in the correlation functions

  for(int i=0; i<NumberOfCorrelationFunctions; i++)
    {
      QMCCorrelationFunctionParameters CP;
      bool ok = CP.read( file , nucCuspReplacement);

      if(!ok)
	{
	  clog << "ERROR: there was a problem reading the correlation function." << endl;
	  clog << "   NumberOfCorrelationFunctions = " << NumberOfCorrelationFunctions << endl;
	  clog << "   NumberOfElectronsUp          = " << NumberOfElectronsUp << endl;
	  clog << "   NumberOfElectronsDown        = " << NumberOfElectronsDown << endl;
	  clog << "   CP #" << i << " = " << endl;
	  clog << CP << endl << endl;
	  clog << "It appears you haven't included enough correlation functions. Put some more in!" << endl;
	  exit(0);
	}

      if( CP.getParticle1Type() == "Electron_up" )
	{
	  if( CP.getParticle2Type() == "Electron_up" )
	    {
	      if( NumberOfElectronsUp < 2 )
		{
		  cerr << "ERROR: Electron_up-Electron_up correlation "
		       << "function loaded when one does not belong!" << endl;
		  exit(0);
		}
	      
	      EupEup = CP;
	    }
	  else if( CP.getParticle2Type() == "Electron_down" )
	    {
	      if( NumberOfElectronsUp < 1 || NumberOfElectronsDown < 1 )
		{
		  cerr << "ERROR: Electron_up-Electron_down correlation "
		       << "function loaded when one does not belong!" << endl
		       << "Maybe you are using an old ckmf file?" << endl;
		  exit(0);
		}

	      EupEdn = CP;
	    }
	  else
	    {
	      // eup - nuclear
	      for( int j=0; j<NucleiTypes.dim1(); j++ )
		if( NucleiTypes(j) == CP.getParticle2Type() )
		    {
		      EupNuclear(j) = CP;
		      break;
		    }
	    }
	}
      else
	{
	  if( CP.getParticle2Type() == "Electron_down" )
	    {
	      if( NumberOfElectronsDown < 2 )
		{
		  cerr << "ERROR: Electron_down-Electron_down correlation "
		       << "function loaded when one does not belong!" << endl;
		  exit(0);
		}
	      
	      EdnEdn = CP;
	    }
	  else
	    {
	      // edn - nuclear
	      for( int j=0; j<NucleiTypes.dim1(); j++ )
		if( NucleiTypes(j) == CP.getParticle2Type() )
		  {
		    if( NumberOfElectronsDown < 1 )
		      {
			cerr << "ERROR: Electron_down-Nuclear correlation "
			     << "function loaded when one does not belong!"
			     << endl;
			exit(0);
		      }
		    EdnNuclear(j) = CP;
		    break;
		  }
	    }
	}
    }

  // if the particles are linked make the down electron parameters equal
  // to the up electron parameters.
  NumberOfEEParameters   = 0;
  NumberOfNEParameters   = 0;
  NumberOfNEupParameters = 0;
  NumberOfEEParameters += EupEup.getTotalNumberOfParameters();
  NumberOfEEParameters += EupEdn.getTotalNumberOfParameters();

  for(int i=0; i<EupNuclear.dim1(); i++)
    {
      NumberOfNEParameters += 
	EupNuclear(i).getTotalNumberOfParameters();
    }

  NumberOfNEupParameters = NumberOfNEParameters;

  if( EquivalentElectronUpDownParams )
    {
      // Now set the things equal that need to be

      if( NumberOfElectronsUp > 1 && NumberOfElectronsDown > 1 )
	{
	  EdnEdn = EupEup;
	  EdnEdn.setParticle1Type("Electron_down");
	  EdnEdn.setParticle2Type("Electron_down");
	}

      if( NumberOfElectronsDown > 0 )
	{
	  EdnNuclear = EupNuclear;
	  for (int i=0; i<EdnNuclear.dim1(); i++)
	    EdnNuclear(i).setParticle1Type("Electron_down");
	}
    } else {

      NumberOfEEParameters += EdnEdn.getTotalNumberOfParameters();

      for(int i=0; i<EdnNuclear.dim1(); i++)
	{
	  NumberOfNEParameters += 
	    EdnNuclear(i).getTotalNumberOfParameters();
	}
    }
}

QMCCorrelationFunctionParameters * QMCJastrowParameters::
                                         getElectronUpElectronDownParameters()
{
  return &EupEdn;
}

QMCCorrelationFunctionParameters * QMCJastrowParameters::
                                           getElectronUpElectronUpParameters()
{
  return &EupEup;
}

QMCCorrelationFunctionParameters * QMCJastrowParameters::
                                        getElectronDownElectronDownParameters()
{
  return &EdnEdn;
}

Array1D<QMCCorrelationFunctionParameters> * QMCJastrowParameters::
                                               getElectronUpNuclearParameters()
{
  return &EupNuclear;
}

Array1D<QMCCorrelationFunctionParameters> * QMCJastrowParameters::
                                             getElectronDownNuclearParameters()
{
  return &EdnNuclear;
}

Array1D<string> * QMCJastrowParameters::getNucleiTypes()
{
  return &NucleiTypes;
}

QMCJastrowParameters::QMCJastrowParameters()
{
}

QMCJastrowParameters::QMCJastrowParameters(const QMCJastrowParameters & rhs)
{
  *this = rhs;
}

ostream & operator<<(ostream &strm, QMCJastrowParameters & rhs)
{
  strm << "&Jastrow" << endl;

  if( rhs.NumberOfElectronsDown > 0 && rhs.NumberOfElectronsUp > 0 )
    {
      strm << rhs.EupEdn << endl;
    }

  if( rhs.NumberOfElectronsUp > 1 )
    {
      strm << rhs.EupEup << endl;
    }

  if( rhs.NumberOfElectronsDown > 1 )
    {
      strm << rhs.EdnEdn << endl;
    }

  for(int i=0; i<rhs.EupNuclear.dim1(); i++)
    {
      strm << rhs.EupNuclear(i) << endl;
      
      if( rhs.NumberOfElectronsDown > 0 )
	{
	  strm << rhs.EdnNuclear(i) << endl;
	}
    }

  strm << "&Jastrow" << endl;
  
  return strm;
}





