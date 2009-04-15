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
  NumberOfEupEdnParameters   = rhs.NumberOfEupEdnParameters;
  NumberOfEupEupParameters   = rhs.NumberOfEupEupParameters;
  NumberOfEdnEdnParameters   = rhs.NumberOfEdnEdnParameters;
  NumberOfNEParameters       = rhs.NumberOfNEParameters;
  NumberOfNEupParameters     = rhs.NumberOfNEupParameters;

  EupNuclear = rhs.EupNuclear;
  EdnNuclear = rhs.EdnNuclear;
  EupEdn = rhs.EupEdn;
  EupEup = rhs.EupEup;
  EdnEdn = rhs.EdnEdn;

  if (globalInput.flags.use_three_body_jastrow == 1)
    {
      NumberOfNEupEdnParameters = rhs.NumberOfNEupEdnParameters;
      NumberOfNEupEupParameters = rhs.NumberOfNEupEupParameters;
      NumberOfNEdnEdnParameters = rhs.NumberOfNEdnEdnParameters;

      EupEdnNuclear = rhs.EupEdnNuclear;
      EupEupNuclear = rhs.EupEupNuclear;
      EdnEdnNuclear = rhs.EdnEdnNuclear;
    }

  EquivalentElectronUpDownParams = rhs.EquivalentElectronUpDownParams;
  NucleiTypes = rhs.NucleiTypes;
  NumberOfElectronsUp = rhs.NumberOfElectronsUp;
  NumberOfElectronsDown = rhs.NumberOfElectronsDown;
}

void QMCJastrowParameters::print(ostream & strm)
{
  if(globalInput.flags.set_debug == 1)
    return;

  if(EupEdn.isUsed())
    {
      strm << "EupEdn:" << endl;
      EupEdn.getCorrelationFunction()->print(strm);
      strm << endl;
    }

  if(EupEup.isUsed())
    {
      strm << "EupEup:" << endl;
      EupEup.getCorrelationFunction()->print(strm);
      strm << endl;
    }
  if(!EquivalentElectronUpDownParams &&
     EdnEdn.isUsed())
    {
      strm << "EdnEdn:" << endl;
      EdnEdn.getCorrelationFunction()->print(strm);    
      strm << endl;
    }

  for(int i=0; i<EupNuclear.dim1(); i++) 
    {
      if(EupNuclear(i).isUsed())
	{
	  strm << "EupNuclear(" << NucleiTypes(i) << "):"
	       << endl;
	  EupNuclear(i).getCorrelationFunction()->print(strm);    
	  strm << endl;
	}

      if(!EquivalentElectronUpDownParams &&
	 EdnNuclear(i).isUsed())
	{
	  strm << "EdnNuclear(" << NucleiTypes(i) << "):" 
	       << endl;
	  EdnNuclear(i).getCorrelationFunction()->print(strm);    
	  strm << endl;
	}
    }

  for(int i=0; i<EupEdnNuclear.dim1(); i++) 
    {
      if(EupEdnNuclear(i).isUsed())
	{
	  if(globalInput.flags.link_NEE_Jastrows == 2)
	    strm << "EE";
	  else
	    strm << "EupEdn";

	  strm << "Nuclear(" << NucleiTypes(i) << "):" << endl
	       << " Total Parameters = " << EupEdnNuclear(i).getNumberOfTotalParameters() << endl
	       << " Free Parameters = " << EupEdnNuclear(i).getNumberOfFreeParameters() 
	       << endl;
	  EupEdnNuclear(i).getThreeBodyCorrelationFunction()->print(strm);
	  strm << endl;
	}

      if(globalInput.flags.link_NEE_Jastrows < 2 &&
	 EupEupNuclear(i).isUsed())
	{
	  strm << "EupEupNuclear(" << NucleiTypes(i) << "):" << endl
	       << " Total Parameters = " << EupEupNuclear(i).getNumberOfTotalParameters() << endl
	       << " Free Parameters = " << EupEupNuclear(i).getNumberOfFreeParameters() 
	       << endl;
	  EupEupNuclear(i).getThreeBodyCorrelationFunction()->print(strm);
	  strm << endl;
	}
      
      if(globalInput.flags.link_NEE_Jastrows < 1 &&
	 EdnEdnNuclear(i).isUsed())
	{
	  strm << "EdnEdnNuclear(" << NucleiTypes(i) << "):" << endl
	       << " Total Parameters = " << EdnEdnNuclear(i).getNumberOfTotalParameters() << endl
	       << " Free Parameters = " << EdnEdnNuclear(i).getNumberOfFreeParameters() 
	       << endl;
	  EdnEdnNuclear(i).getThreeBodyCorrelationFunction()->print(strm);
	  strm << endl;	  
	}
    }
}

void QMCJastrowParameters::setJWParameters(Array1D<double> & params, int shift)
{
  int Counter = 0;
  int CurrentNumberOfParams = 0;
  Array1D<double> temp;

  // Set eup edn
  if(getNumberEupEdnParameters() > 0 &&
     NumberOfElectronsUp > 0 &&
     NumberOfElectronsDown > 0 ){
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
  if(getNumberEupEupParameters() > 0 &&
     NumberOfElectronsUp > 1 )
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
  if(getNumberEdnEdnParameters() > 0 &&
     NumberOfElectronsDown > 1 )
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

  if(EquivalentElectronUpDownParams)
    {
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
	    EdnNuclear(i).setParticle1Type("Electron_down");
	}
    }
  else
    {
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


  if (globalInput.flags.use_three_body_jastrow == 1 &&
      globalInput.flags.optimize_NEE_Jastrows == 1)
    {
      for (int nuc=0; nuc<EupEdnNuclear.dim1(); nuc++)
	{
	  if (NumberOfElectronsUp > 0 && NumberOfElectronsDown > 0)
	    {
	      CurrentNumberOfParams = 
		EupEdnNuclear(nuc).getNumberOfFreeParameters();
	      
	      temp.allocate(CurrentNumberOfParams);
	      
	      for (int j=0; j<temp.dim1(); j++)
		{
		  temp(j) = params(Counter + shift);
		  Counter++;
		}
	      
	      EupEdnNuclear(nuc).setFreeParameters( temp );
	    }

	  //up up electrons
	  if (NumberOfElectronsUp > 1)
	    {
	      if(globalInput.flags.link_NEE_Jastrows < 2)
		{
		  CurrentNumberOfParams =
		    EupEupNuclear(nuc).getNumberOfFreeParameters();
		  
		  temp.allocate(CurrentNumberOfParams);
		  
		  for (int j=0; j<temp.dim1(); j++)
		    {
		      temp(j) = params(Counter + shift);
		      Counter++;
		    }
		  
		  EupEupNuclear(nuc).setFreeParameters( temp );		
		}
	      else
		{
		  EupEupNuclear = EupEupNuclear;
		  for (int i=0; i<EdnEdnNuclear.dim1(); i++)
		    {
		      EupEupNuclear(i).setParticle2Type("Electron_up");
		      EupEupNuclear(i).setParticle3Type("Electron_up");
		    }
		}
	    }

	  //down down electrons
	  if(NumberOfElectronsDown > 1)
	    {
	      if(globalInput.flags.link_NEE_Jastrows < 1)
		{
		  CurrentNumberOfParams =
		    EdnEdnNuclear(nuc).getNumberOfFreeParameters();
		  
		  temp.allocate(CurrentNumberOfParams);
		  
		  for (int j=0; j<temp.dim1(); j++)
		    {
		      temp(j) = params(Counter + shift);
		      Counter++;
		    }
		  
		  EdnEdnNuclear(nuc).setFreeParameters( temp );
		}
	      else
		{
		  EdnEdnNuclear = EupEupNuclear;
		  for (int i=0; i<EdnEdnNuclear.dim1(); i++)
		    {
		      EdnEdnNuclear(i).setParticle2Type("Electron_down");
		      EdnEdnNuclear(i).setParticle3Type("Electron_down");
		    }
		}
	    }
	}
    }  
}

int QMCJastrowParameters::getNumberJWParameters()
{
  return
    getNumberEEParameters() +
    getNumberNEParameters() + 
    getNumberNEupEdnParameters() +
    getNumberNEupEupParameters() + 
    getNumberNEdnEdnParameters();
}

int QMCJastrowParameters::getNumberEEParameters()
{
  return getNumberEupEdnParameters() +
    getNumberEupEupParameters() + 
    getNumberEdnEdnParameters();
}

int QMCJastrowParameters::getNumberEupEdnParameters()
{
  if(globalInput.flags.optimize_UD_Jastrows == 0)
    return 0;
  return NumberOfEupEdnParameters;
}

int QMCJastrowParameters::getNumberEupEupParameters()
{
  if(globalInput.flags.optimize_UU_Jastrows == 0)
    return 0;
  return NumberOfEupEupParameters;
}

int QMCJastrowParameters::getNumberEdnEdnParameters()
{
  if(globalInput.flags.optimize_DD_Jastrows == 0 ||
     EquivalentElectronUpDownParams)
    return 0;
  return NumberOfEdnEdnParameters;
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

int QMCJastrowParameters::getNumberNEupEdnParameters()
{
  if(globalInput.flags.optimize_NEE_Jastrows == 0)
    return 0;
  return NumberOfNEupEdnParameters;
}

int QMCJastrowParameters::getNumberNEupEupParameters()
{
  if(globalInput.flags.optimize_NEE_Jastrows == 0 ||
     globalInput.flags.link_NEE_Jastrows > 1)
    return 0;
  return NumberOfNEupEupParameters;
}

int QMCJastrowParameters::getNumberNEdnEdnParameters()
{
  if(globalInput.flags.optimize_NEE_Jastrows == 0 ||
     globalInput.flags.link_NEE_Jastrows > 0)
    return 0;
  return NumberOfNEdnEdnParameters;
}

Array1D<double> QMCJastrowParameters::getJWParameters()
{
  Array1D<double> params(getNumberJWParameters());
  getJWParameters(params,0);
  return params;
}

void QMCJastrowParameters::getJWParameters(Array1D<double> & params, int shift)
{ 
  int Counter = 0;
  Array1D<double> temp;

  // Put in eup edn params
  if(getNumberEupEdnParameters() > 0 &&
     NumberOfElectronsUp > 0 &&
     NumberOfElectronsDown > 0)
    {
      temp = EupEdn.getParameters();
      
      for(int i=0; i<temp.dim1(); i++)
	{
	  params(Counter + shift) = temp(i);
	  Counter++;
	}
    }
  
  // Put in eup eup params
  if(getNumberEupEupParameters() > 0 &&
     NumberOfElectronsUp > 1)
    {
      temp = EupEup.getParameters();
      
      for(int i=0; i<temp.dim1(); i++)
	{
	  params(Counter + shift) = temp(i);
	  Counter++;
	}
    }    

  // Put in edn edn params
  if(getNumberEdnEdnParameters() > 0 &&	  
     NumberOfElectronsDown > 1 )
    {
      temp = EdnEdn.getParameters();
      
      for(int i=0; i<temp.dim1(); i++)
	{
	  params(Counter + shift) = temp(i);
	  Counter++;
	}
    }

  // Put in eup nuclear params
  if(getNumberNEParameters() > 0)
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
  if( !EquivalentElectronUpDownParams &&
      NumberOfElectronsDown > 0 )
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

  if (globalInput.flags.optimize_NEE_Jastrows == 1)
    {
      for (int i=0; i<EupEdnNuclear.dim1(); i++)
	{
	  if (NumberOfElectronsUp > 0 && NumberOfElectronsDown > 0)
	    {
	      temp = EupEdnNuclear(i).getFreeParameters();
	      
	      for (int j=0; j<temp.dim1(); j++)
		{
		  params(Counter + shift) = temp(j);
		  Counter++;
		}
	    }
	  
	  if(globalInput.flags.link_NEE_Jastrows < 2 && NumberOfElectronsUp > 1)
	    {
	      temp = EupEupNuclear(i).getFreeParameters();
	      
	      for (int j=0; j<temp.dim1(); j++)
		{
		  params(Counter + shift) = temp(j);
		  Counter++;
		}
	    }
      
	  if(globalInput.flags.link_NEE_Jastrows < 1 && NumberOfElectronsDown > 1)
	    {
	      temp = EdnEdnNuclear(i).getFreeParameters();
	      
	      for (int j=0; j<temp.dim1(); j++)
		{
		  params(Counter + shift) = temp(j);
		  Counter++;
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

      Array1D<Complex> poles;
      try {
        poles = EupEup.getPoles();
      }

      catch(Exception e){
        cout << "Exception: EupEup " << e.getMessage() << endl;
      }

      for(int i=0; i<poles.dim1(); i++)
        {
          allPoles.push_back(poles(i));
        }
    }

  if( NumberOfElectronsDown > 1 && !EquivalentElectronUpDownParams)
    {
      // down down jastrow
      Array1D<Complex> poles;
      try {
        poles = EdnEdn.getPoles();
      }

      catch(Exception e){
        cout << "Exception: EdnEdn " << e.getMessage() << endl;
      }
      
      for(int i=0; i<poles.dim1(); i++)
        {
          allPoles.push_back(poles(i));
        }
    }

  if( NumberOfElectronsUp > 0 && NumberOfElectronsDown > 0 )
    {
      // up down jastrow
      Array1D<Complex> poles;
      try {
        poles = EupEdn.getPoles();
      }

      catch(Exception e){
        cout << "Exception: EupEdn " << e.getMessage() << endl;
      }
      
      for(int i=0; i<poles.dim1(); i++)
        {
          allPoles.push_back(poles(i));
        }
    }

  for(int i=0; i<EupNuclear.dim1(); i++)
    {
      // up nuclear jastrow
      Array1D<Complex> poles;
      try {
        poles = EupNuclear(i).getPoles();
      }

      catch(Exception e){
        cout << "Exception: EupNuclear(" << i << ") " << e.getMessage() 
	     << endl;
      }

      for(int j=0; j<poles.dim1(); j++)
        {
          allPoles.push_back(poles(j));
        }
    }

  if(!EquivalentElectronUpDownParams)
    for(int i=0; i<EdnNuclear.dim1(); i++)
      {
        // down nuclear jastrow
        Array1D<Complex> poles;
        try {
          poles = EdnNuclear(i).getPoles();
        }
        
        catch(Exception e){
          cout << "Exception: EdnNuclear(" << i << ") " << e.getMessage() 
	       << endl;
        }
        
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
	distance = fabs(poles(i).imaginary());
      else
	distance = poles(i).abs();

      if(distance <= 0)
	{
	  cerr << "Warning: distance from real axis = " << distance << " in calculate_penalty_function, can\'t take log" << endl;
	  cerr << "         poles(" << i << ") = " << poles(i) << endl;
	  if(poles(i).real() > 50.0)
	    penalty = 10;
	  else
	    penalty = 1e100;
	  cerr << "         penalty will be set to " << penalty << endl;
	  cerr.flush();
	} 
      else 
	{
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
				bool linkparams, bool nucCuspReplacement, 
				int nelup, int neldn, string runfile)
{
  // Set the internal variable telling if the electrons of different spins
  // should use the same parameters
  EquivalentElectronUpDownParams = linkparams;

  // Save an array containing all of nuclei types
  NucleiTypes = nucleitypes;

  // Save the number of up and down electrons
  NumberOfElectronsUp   = nelup;
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
  
  if(globalInput.flags.link_Jastrow_parameters == 0 && NumberOfElectronsDown > 0 )
    NumberOfCorrelationFunctions += NucleiTypes.dim1();

  if( NumberOfElectronsUp >0 && NumberOfElectronsDown > 0 )
    NumberOfCorrelationFunctions++;

  if( NumberOfElectronsUp > 1 )
    NumberOfCorrelationFunctions++;

  if(globalInput.flags.link_Jastrow_parameters == 0 &&  NumberOfElectronsDown > 1 )
    NumberOfCorrelationFunctions++;

  // Load in the correlation functions
  QMCCorrelationFunctionParameters CP;
  while(CP.read( file , nucCuspReplacement))
    {
      bool foundMatch = false;
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
	      foundMatch = true;
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
	      foundMatch = true;
	      EupEdn = CP;
	    }
	  else
	    {
	      // eup - nuclear
	      for( int j=0; j<NucleiTypes.dim1(); j++ )
		if( NucleiTypes(j) == CP.getParticle2Type() )
		    {
		      foundMatch = true;
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
	      foundMatch = true;	      
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
		    foundMatch = true;
		    EdnNuclear(j) = CP;
		    break;
		  }
	    }
	}

      if(!foundMatch){
	clog << "Warning: Jastrow from input file ("
	     << CP.getParticle1Type() << "," << CP.getParticle2Type() << ") was not used:\n";
	CP.getCorrelationFunction()->print(clog);
	clog << endl << endl;
      }
    }

  // if the particles are linked make the down electron parameters equal
  // to the up electron parameters.
  NumberOfEupEdnParameters   = 0;
  NumberOfEupEupParameters   = 0;
  NumberOfEdnEdnParameters   = 0;
  NumberOfNEParameters       = 0;
  NumberOfNEupParameters     = 0;

  NumberOfEupEupParameters  += EupEup.getTotalNumberOfParameters();
  NumberOfEupEdnParameters  += EupEdn.getTotalNumberOfParameters();

  for(int i=0; i<EupNuclear.dim1(); i++)
    {
      NumberOfNEParameters += 
	EupNuclear(i).getTotalNumberOfParameters();
    }

  NumberOfNEupParameters = NumberOfNEParameters;

  if( EquivalentElectronUpDownParams )
    {
      // Now set the things equal that need to be

      if( NumberOfElectronsUp > 1 && NumberOfElectronsDown > 1)
	{
	  EdnEdn = EupEup;
	  EdnEdn.setParticle1Type("Electron_down");
	  EdnEdn.setParticle2Type("Electron_down");
	}

      if( NumberOfElectronsDown > 0 )
	{
	  for (int i=0; i<EupNuclear.dim1(); i++)
	    {
	      EdnNuclear(i) = EupNuclear(i);
	      EdnNuclear(i).setParticle1Type("Electron_down");
	    }
	}
    } 
  else 
    {
      NumberOfEdnEdnParameters  += EdnEdn.getTotalNumberOfParameters();

      for(int i=0; i<EdnNuclear.dim1(); i++)
	{
	  NumberOfNEParameters += 
	    EdnNuclear(i).getTotalNumberOfParameters();
	}
    }

  // Read in the QMCThreeBodyCorrelationFunctionParameters
  int NumberOfThreeBodyCorrelationFunctions = 0;
  EupEdnNuclear.allocate( NucleiTypes.dim1() );
  EupEupNuclear.allocate( NucleiTypes.dim1() );
  EdnEdnNuclear.allocate( NucleiTypes.dim1() );
  
  if (NumberOfElectronsUp > 0 && NumberOfElectronsDown > 0)
    NumberOfThreeBodyCorrelationFunctions += NucleiTypes.dim1();
  
  if (NumberOfElectronsUp > 1 &&
      globalInput.flags.link_NEE_Jastrows < 2)
    NumberOfThreeBodyCorrelationFunctions += NucleiTypes.dim1();
  
  if (NumberOfElectronsDown > 1 &&
      globalInput.flags.link_NEE_Jastrows < 1)
    NumberOfThreeBodyCorrelationFunctions += NucleiTypes.dim1();
  
  // Load in the three body correlation functions
  QMCThreeBodyCorrelationFunctionParameters CP3;
  while(CP3.read( file ))
    {
      bool foundMatch = false;
      if( CP3.getParticle2Type() == "Electron_up" )
	{
	  if( CP3.getParticle3Type() == "Electron_down" )
	    {
	      if( NumberOfElectronsUp < 1 && NumberOfElectronsDown < 1 )
		{
		  cerr << "ERROR: Electron_up-Electron_down nuclear "
		       << "{correlation function loaded when one does not "
		       << "belong!" << endl;
		  exit(0);
		}
	      for (int j=0; j<NucleiTypes.dim1(); j++)
		if (NucleiTypes(j) == CP3.getParticle1Type() )
		  {
		    foundMatch = true;
		    EupEdnNuclear(j) = CP3;
		    break;
		  }
	    }
	  else if( CP3.getParticle3Type() == "Electron_up" )
	    {
	      if(NumberOfElectronsUp < 2)
		{
		  cerr << "ERROR: Electron_up-Electron_up nuclear "
		       << "correlation function loaded when one does not "
		       << "belong!" << endl
		       << "Maybe you are using an old ckmf file?" << endl;
		  exit(0);
		}
	      for (int j=0; j<NucleiTypes.dim1(); j++)
		if (NucleiTypes(j) == CP3.getParticle1Type() )
		  {
		    foundMatch = true;
		    EupEupNuclear(j) = CP3;
		    break;
		  }
	    }
	}
      else if (CP3.getParticle2Type() == "Electron_down" && 
	       CP3.getParticle3Type() == "Electron_down")
	{
	  if( NumberOfElectronsDown < 2 )
	    {
	      cerr << "ERROR: Electron_down-Electron_down nuclear "
		   << "correlation function loaded when one does not "
		   << "belong!" << endl;
	      exit(0);
	    }
	  
	  for( int j=0; j<NucleiTypes.dim1(); j++ )
	    if( NucleiTypes(j) == CP3.getParticle1Type() )
	      {
		foundMatch = true;
		EdnEdnNuclear(j) = CP3;
		break;
	      }
	}
      
      if(!foundMatch){
	clog << "Warning: Jastrow from input file ("
	     << CP3.getParticle1Type() << "," << CP3.getParticle2Type() << ","
	     << CP3.getParticle3Type() << ") was not used:\n";
	CP3.getThreeBodyCorrelationFunction()->print(clog);
	clog << endl << endl;
      }
    }
  
  // if the particles are linked make the down electron parameters equal
  // to the up electron parameters.
  NumberOfNEupEdnParameters = 0;
  NumberOfNEupEupParameters = 0;
  NumberOfNEdnEdnParameters = 0;
  
  for (int i=0; i<EupEdnNuclear.dim1(); i++)
    NumberOfNEupEdnParameters += EupEdnNuclear(i).getNumberOfFreeParameters();
  
  for (int i=0; i<EupEupNuclear.dim1(); i++)
    NumberOfNEupEupParameters += EupEupNuclear(i).getNumberOfFreeParameters();
  
  for (int i=0; i<EdnEdnNuclear.dim1(); i++)
    NumberOfNEdnEdnParameters += EdnEdnNuclear(i).getNumberOfFreeParameters();
  
  if(globalInput.flags.link_NEE_Jastrows == 1)
    {
      // Now set the things equal that need to be
      if( NumberOfElectronsUp > 1 && NumberOfElectronsDown > 1 )
	for (int i=0; i<EupEupNuclear.dim1(); i++)
	  {
	    EdnEdnNuclear(i) = EupEupNuclear(i);
	    EdnEdnNuclear(i).setParticle2Type("Electron_down");
	    EdnEdnNuclear(i).setParticle3Type("Electron_down");
	  }
    }
  else if(globalInput.flags.link_NEE_Jastrows == 2)
    {
      // Now set the things equal that need to be
      if( NumberOfElectronsUp > 1 && NumberOfElectronsDown > 1 )
	for (int i=0; i<EupEupNuclear.dim1(); i++)
	  {
	    EdnEdnNuclear(i) = EupEdnNuclear(i);
	    EupEupNuclear(i) = EupEdnNuclear(i);
	    EdnEdnNuclear(i).setParticle2Type("Electron_down");
	    EdnEdnNuclear(i).setParticle3Type("Electron_down");
	    EupEupNuclear(i).setParticle2Type("Electron_up");
	    EupEupNuclear(i).setParticle3Type("Electron_up");
	  }
    }

  bool use3 = false;
  for(int i=0; i<EupEdnNuclear.dim1(); i++)
    if(EupEdnNuclear(i).isUsed()) use3 = true;
  for(int i=0; i<EdnEdnNuclear.dim1(); i++)
    if(EdnEdnNuclear(i).isUsed()) use3 = true;
  for(int i=0; i<EupEupNuclear.dim1(); i++)
    if(EupEupNuclear(i).isUsed()) use3 = true;

  if(!use3)
    {
      if(globalInput.flags.use_three_body_jastrow == 1)
	clog << "Warning: setting use_three_body_jastrow = 0 since you didn't load any Jastrows.\n";
      /*
	The code should work fine if the switch is on,
	but turning it off will help some parts save time.
      */
      globalInput.flags.use_three_body_jastrow = 0;
      globalInput.flags.optimize_NEE_Jastrows  = 0;
      
      EupEdnNuclear.deallocate();
      EupEupNuclear.deallocate();
      EdnEdnNuclear.deallocate();
    } else {
      globalInput.flags.use_three_body_jastrow = 1;
    }

  print(clog);
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

Array1D<QMCThreeBodyCorrelationFunctionParameters> * QMCJastrowParameters::
                                   getElectronUpElectronDownNuclearParameters()
{
  return &EupEdnNuclear;
}

Array1D<QMCThreeBodyCorrelationFunctionParameters> * QMCJastrowParameters::
                                     getElectronUpElectronUpNuclearParameters()
{
  return &EupEupNuclear;
}

Array1D<QMCThreeBodyCorrelationFunctionParameters> * QMCJastrowParameters::
                                 getElectronDownElectronDownNuclearParameters()
{
  return &EdnEdnNuclear;
}

Array1D<string> * QMCJastrowParameters::getNucleiTypes()
{
  return &NucleiTypes;
}

QMCJastrowParameters::QMCJastrowParameters()
{
  NumberOfEupEdnParameters   = 0;
  NumberOfEupEupParameters   = 0;
  NumberOfEdnEdnParameters   = 0;
  NumberOfNEParameters       = 0;
  NumberOfNEupParameters     = 0;
  this->NumberOfNEParameters = 0;
}

QMCJastrowParameters::QMCJastrowParameters(const QMCJastrowParameters & rhs)
{
  *this = rhs;
}

ostream & operator<<(ostream &strm, QMCJastrowParameters & rhs)
{
  strm << "&Jastrow" << endl;

  if( rhs.NumberOfElectronsDown > 0 && rhs.NumberOfElectronsUp > 0 &&
      rhs.EupEdn.isUsed())
    {
      strm << rhs.EupEdn << endl;
    }

  if( rhs.NumberOfElectronsUp > 1 &&
      rhs.EupEup.isUsed())
    {
      strm << rhs.EupEup << endl;
    }

  if( rhs.NumberOfElectronsDown > 1 &&
      !rhs.EquivalentElectronUpDownParams && 
      rhs.EdnEdn.isUsed())
    {
      strm << rhs.EdnEdn << endl;
    }

  for(int i=0; i<rhs.EupNuclear.dim1(); i++)
    {
      if(rhs.EupNuclear(i).isUsed())
	strm << rhs.EupNuclear(i) << endl;
      
      if( rhs.NumberOfElectronsDown > 0 &&
	  !rhs.EquivalentElectronUpDownParams)
	{
	  if(rhs.EdnNuclear(i).isUsed())
	    strm << rhs.EdnNuclear(i) << endl;
	}
    }

  if (globalInput.flags.use_three_body_jastrow == 1)
    for (int i=0; i<rhs.EupEdnNuclear.dim1(); i++)
      {
	if( rhs.NumberOfElectronsDown > 0 && rhs.NumberOfElectronsUp > 0 &&
	    rhs.EupEdnNuclear(i).isUsed())
	  strm << rhs.EupEdnNuclear(i) << endl;
  
	if ( rhs.NumberOfElectronsUp > 1 &&
	     globalInput.flags.link_NEE_Jastrows < 2 &&
	     rhs.EupEupNuclear(i).isUsed())
	  strm << rhs.EupEupNuclear(i) << endl;
	
	if ( rhs.NumberOfElectronsDown > 1 &&
	     globalInput.flags.link_NEE_Jastrows < 1 &&
	     rhs.EdnEdnNuclear(i).isUsed())
	  strm << rhs.EdnEdnNuclear(i) << endl;
      }

  strm << "&Jastrow" << endl;
  
  return strm;
}





