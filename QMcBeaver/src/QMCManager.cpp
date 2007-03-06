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

#include "QMCManager.h"

/**
   Random.h has an extern call referencing this.
*/
Random ran;

bool QMCManager::SIGNAL_SAYS_QUIT         = false;  
bool QMCManager::REDUCE_ALL_NOW           = false;
bool QMCManager::INCREASE_TIME            = false;
bool QMCManager::PRINT_SIG_INFO           = false;

QMCManager::QMCManager()
{
  localTimers.getTotalTimeStopwatch()->start();
}

QMCManager::~QMCManager()
{
  finalizeOutputs();
}

void QMCManager::initialize( int argc, char **argv )
{
  initializeMPI();
  string runfile = sendAllProcessorsInputFileName( argv );
  
  // Load the input data
  Input.read( runfile );
  ran.initialize(Input.flags.iseed,Input.flags.my_rank);

  initializeOutputs();
  
  if ( Input.flags.calculate_bf_density == 1 )
  {
    bool calcDensity = true;
    Properties_total.setCalcDensity( calcDensity,
                                     Input.WF.getNumberBasisFunctions() );
  }
  
  if(Input.flags.nuclear_derivatives != "none")
  {
    if(Input.flags.nuclear_derivatives != "bin_force_density")
      {
	fwProperties_total.setCalcForces( true, Input.Molecule.getNumberAtoms(),3 );
      } else {
	fwProperties_total.setCalcForces( true, QMCNuclearForces::getNumBins(), 1 );
      }
  }

  QMCnode.initialize( &Input );

  // Set equilibrating = true so that the run equilibrates before taking
  // data
  
  if(  Input.flags.equilibrate_first_opt_step == 0 )
    equilibrating = false;
  else
    equilibrating = true;
  
  done = false;
  
  iteration = 0;
  
  // Initialize the calculation state
  initializeCalculationState();
}

void QMCManager::initializeMPI()
{
#ifdef PARALLEL
  // Create MPI Community
  
  if(  MPI_Comm_size( MPI_COMM_WORLD, &Input.flags.nprocs )  )
  {
    cerr << "ERROR: MPI_Comm_size Error" << endl;
    exit( 1 );
  }
  
  // Set Processor Rank
  if(  MPI_Comm_rank( MPI_COMM_WORLD, &Input.flags.my_rank )  )
  {
    cerr << "ERROR: MPI_Comm_rank Error" << endl;
    exit( 1 );
  }
  
#else
  Input.flags.nprocs = 1;
  
  Input.flags.my_rank = 0;
  
#endif
}

void QMCManager::initializeOutputs()
{
  Input.openConfigFile();

  if(  Input.flags.my_rank == 0 )
  {
    QMCCopyright copyright;
    
    // Allocate result stream
    qmcRslts = new ofstream( Input.flags.results_file_name.c_str() );
    
    ( *qmcRslts ).setf( ios::fixed,ios::floatfield );
    
    ( *qmcRslts ).precision( 10 );
    
    if(  !( *qmcRslts )  )
    {
      cerr << "ERROR opening " << Input.flags.results_file_name << endl;
      exit( 0 );
    }
    
    *qmcRslts << copyright;
    
    // Allocate qmc output stream
    qmcOut = new ofstream( Input.flags.output_file_name.c_str() );
    
    ( *qmcOut ).setf( ios::fixed,ios::floatfield );
    
    ( *qmcOut ).precision( 10 );
    
    if(  !( *qmcOut )  )
    {
      cerr << "ERROR opening " << Input.flags.output_file_name << endl;
      exit( 0 );
    }
    
    *qmcOut << copyright;
    
    cout.setf( ios::fixed,ios::floatfield );
    
    cout.precision( 10 );
    
      if( !true )
      cout << copyright;
  }
  else
  {
    qmcRslts = 0;
    qmcOut   = 0;
  }
}

void QMCManager::finalizeOutputs()
{
  if(  Input.flags.my_rank == 0 )
  {
    // Close and deallocate the result stream
    ( *qmcRslts ).close();
    delete qmcRslts;
    qmcRslts = 0;
    
    // Close and deallocate the qmc result stream
    ( *qmcOut ).close();
    delete qmcOut;
    qmcOut = 0;
  }
}

void QMCManager::finalize()
{
  //stop timing and package up the timer data
  localTimers.stop();
  
  if ( Input.flags.use_equilibration_array == 1 )
  {
    QMCnode.stopTimers();
    *localTimers.getPropagationStopwatch()  =
      *QMCnode.getPropagationStopwatch();
    *localTimers.getEquilibrationStopwatch()  =
      *localTimers.getEquilibrationStopwatch()  +
      *QMCnode.getEquilibrationStopwatch();
  }
  
#ifdef PARALLEL
  MPI_Reduce( &localTimers,&globalTimers,1,QMCStopwatches::MPI_TYPE,
              QMCStopwatches::MPI_REDUCE,0,MPI_COMM_WORLD );

#else
  globalTimers = localTimers;
  
#endif
}

void QMCManager::sendAllProcessorsACommand( int itag )
{
#ifdef PARALLEL
  localTimers.getSendCommandStopwatch()->start();
  
  MPI_Request *request = new MPI_Request[ Input.flags.nprocs-1 ];
  
  // Send the command integer to each processor
  
  for( int i=1;i<Input.flags.nprocs;i++ )
  {
    MPI_Isend( &itag,1,MPI_INT,i,itag,MPI_COMM_WORLD,&request[ i-1 ] );
  }
  
  // Wait for all of the sends to be completed
  MPI_Status *status = new MPI_Status[ Input.flags.nprocs-1 ];
  
  if(  MPI_Waitall( Input.flags.nprocs-1, request, status )  )
  {
    cerr << "ERROR: MPI_Waitall Error (QMCManager::MPI_command_send)"
    << endl;
  }
  
  delete []  request;
  delete []  status;
  
  localTimers.getSendCommandStopwatch()->stop();
#endif
}

void QMCManager::gatherProperties()
{
  Properties_total.zeroOut();
  fwProperties_total.zeroOut();
#ifdef PARALLEL
  localTimers.getGatherPropertiesStopwatch()->start();

  MPI_Reduce( QMCnode.getProperties(), &Properties_total, 1,
              QMCProperties::MPI_TYPE, QMCProperties::MPI_REDUCE, 0, MPI_COMM_WORLD );

  //reduce over QMCFutureWalkingProperties
  //maybe this doesn't have to happen so often...
  for(int i=0; i<fwProperties_total.props.dim1(); i++)
    {
      MPI_Reduce( QMCnode.getFWProperties()->props[i].array(),
		  fwProperties_total.props[i].array(),
		  fwProperties_total.props[i].dim1(),
		  QMCProperty::MPI_TYPE, QMCProperty::MPI_REDUCE, 0, MPI_COMM_WORLD );
    }

  localTimers.getGatherPropertiesStopwatch()->stop();
#else
  Properties_total = *QMCnode.getProperties();
  fwProperties_total = *QMCnode.getFWProperties();
#endif
}

void QMCManager::synchronizeDMCEnsemble()
{
  Properties_total.zeroOut();
#ifdef PARALLEL
  localTimers.getCommunicationSynchronizationStopwatch()->start();
  
  // Collect the global properties and send to all nodes
  
  MPI_Allreduce( QMCnode.getProperties(), &Properties_total, 1,
                 QMCProperties::MPI_TYPE, QMCProperties::MPI_REDUCE, MPI_COMM_WORLD );
  
  localTimers.getCommunicationSynchronizationStopwatch()->stop();
#else
  Properties_total = *QMCnode.getProperties();
#endif
}

void QMCManager::gatherDensities()
{
#ifdef PARALLEL
  localTimers.getGatherPropertiesStopwatch()->start();
  
  Array1D<QMCProperty> localChiDensity = QMCnode.getProperties() ->chiDensity;
  
  MPI_Reduce( QMCnode.getProperties()->chiDensity.array(),
              Properties_total.chiDensity.array(),Input.WF.getNumberBasisFunctions(),
              QMCProperty::MPI_TYPE,QMCProperty::MPI_REDUCE,0,MPI_COMM_WORLD );
  
  localTimers.getGatherPropertiesStopwatch()->stop();
#else
  
  for ( int i=0; i<Input.WF.getNumberBasisFunctions(); i++ )
  {
    Properties_total.chiDensity(i)=QMCnode.getProperties()->chiDensity( i );
  }
  
#endif
}

void QMCManager::gatherForces()
{
#ifdef PARALLEL
  localTimers.getGatherPropertiesStopwatch()->start();

  int numProps = QMCnode.getFWProperties()->nuclearForces(0).dim1() *
                 QMCnode.getFWProperties()->nuclearForces(0).dim2();

  for(int fw=0; fw<QMCnode.getFWProperties()->nuclearForces.dim1(); fw++)
  {
    MPI_Reduce( QMCnode.getFWProperties()->nuclearForces(fw).array(),
                fwProperties_total.nuclearForces(fw).array(),
                numProps,
                QMCProperty::MPI_TYPE,QMCProperty::MPI_REDUCE,0,MPI_COMM_WORLD );
  }
  localTimers.getGatherPropertiesStopwatch()->stop();
#else

  fwProperties_total.nuclearForces=QMCnode.getFWProperties()->nuclearForces;  

#endif
}

void QMCManager::gatherHistograms()
{
  int nalpha = Input.WF.getNumberAlphaElectrons();
  int nbeta = Input.WF.getNumberBetaElectrons();
  int nucleiTypes = Input.Molecule.NucleiTypes.dim1();
  
#ifdef PARALLEL
  localTimers.getGatherPropertiesStopwatch()->start();
  
  if (nalpha > 1 || nbeta > 1)
    {
      pllSpinHistogram_total.allocate(5000);
      pllSpinHistogram_total = 0.0;
    
      MPI_Reduce( QMCnode.getPllSpinHistogram()->array(),
     pllSpinHistogram_total.array(),5000,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD );
    
      if (Input.flags.my_rank != 0)
	pllSpinHistogram_total.deallocate();
    }
  
  if (nalpha > 0 && nbeta > 0)
    {
      oppSpinHistogram_total.allocate(5000);
      oppSpinHistogram_total = 0.0;
    
      MPI_Reduce( QMCnode.getOppSpinHistogram()->array(),
     oppSpinHistogram_total.array(),5000,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD );
    
      if (Input.flags.my_rank != 0)
	oppSpinHistogram_total.deallocate();
    }
  
  if (nalpha > 0)
    {
      alphaHistograms_total.allocate(nucleiTypes);
      for (int k=0; k<nucleiTypes; k++)
	{
	  alphaHistograms_total(k).allocate(5000);
	  alphaHistograms_total(k) = 0.0;
	  MPI_Reduce( (*QMCnode.getAlphaHistograms())(k).array(),
    alphaHistograms_total(k).array(),5000,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	  if (Input.flags.my_rank != 0)
	    alphaHistograms_total(k).deallocate();
	}
      if (Input.flags.my_rank != 0)
	alphaHistograms_total.deallocate();
    }
  
  if (nbeta > 0)
    {
      betaHistograms_total.allocate(nucleiTypes);
      for (int k=0; k<nucleiTypes; k++)
	{
	  betaHistograms_total(k).allocate(5000);
	  betaHistograms_total(k) = 0.0;
	  MPI_Reduce( (*QMCnode.getBetaHistograms())(k).array(),
     betaHistograms_total(k).array(),5000,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	  if (Input.flags.my_rank != 0)
	    betaHistograms_total(k).deallocate();
	}
      if (Input.flags.my_rank != 0)
	betaHistograms_total.deallocate();
    }
  
  if (nalpha > 1 || nbeta > 1)
    {
      if (Input.flags.writePllxCorrelationDiagram == 1)
	{
	  pllxCorrelationDiagram_total.allocate(1000);
	  for (int i=0; i<1000; i++)
	    {
	      pllxCorrelationDiagram_total(i).allocate(1000);
	      pllxCorrelationDiagram_total(i) = 0.0;
	      MPI_Reduce( (*QMCnode.getPllxCorrelationDiagram())(i).array(),
             pllxCorrelationDiagram_total(i).array(),1000,MPI_DOUBLE,MPI_SUM,0,
                                                               MPI_COMM_WORLD);
	      if (Input.flags.my_rank != 0)
		pllxCorrelationDiagram_total(i).deallocate();
	    }
	  if (Input.flags.my_rank != 0)
	    pllxCorrelationDiagram_total.deallocate();
	}
      if (Input.flags.writePllyCorrelationDiagram == 1)
	{
	  pllyCorrelationDiagram_total.allocate(1000);
	  for (int i=0; i<1000; i++)
	    {
	      pllyCorrelationDiagram_total(i).allocate(1000);
	      pllyCorrelationDiagram_total(i) = 0.0;
	      MPI_Reduce( (*QMCnode.getPllyCorrelationDiagram())(i).array(),
             pllyCorrelationDiagram_total(i).array(),1000,MPI_DOUBLE,MPI_SUM,0,
                                                               MPI_COMM_WORLD);
	      if (Input.flags.my_rank != 0)
		pllyCorrelationDiagram_total(i).deallocate();
	    }
	  if (Input.flags.my_rank != 0)
	    pllyCorrelationDiagram_total.deallocate();
	}
      if (Input.flags.writePllzCorrelationDiagram == 1)
	{
	  pllzCorrelationDiagram_total.allocate(1000);
	  for (int i=0; i<1000; i++)
	    {
	      pllzCorrelationDiagram_total(i).allocate(1000);
	      pllzCorrelationDiagram_total(i) = 0.0;
	      MPI_Reduce( (*QMCnode.getPllzCorrelationDiagram())(i).array(),
             pllzCorrelationDiagram_total(i).array(),1000,MPI_DOUBLE,MPI_SUM,0,
                                                               MPI_COMM_WORLD);
	      if (Input.flags.my_rank != 0)
		pllzCorrelationDiagram_total(i).deallocate();
	    }
	  if (Input.flags.my_rank != 0)
	    pllzCorrelationDiagram_total.deallocate();
	}
    }

  if (nalpha > 0 && nbeta > 0)
    {
      if (Input.flags.writeOppxCorrelationDiagram == 1)
	{
	  oppxCorrelationDiagram_total.allocate(1000);
	  for (int i=0; i<1000; i++)
	    {
	      oppxCorrelationDiagram_total(i).allocate(1000);
	      oppxCorrelationDiagram_total(i) = 0.0;
	      MPI_Reduce( (*QMCnode.getOppxCorrelationDiagram())(i).array(),
             oppxCorrelationDiagram_total(i).array(),1000,MPI_DOUBLE,MPI_SUM,0,
                                                               MPI_COMM_WORLD);
	      if (Input.flags.my_rank != 0)
		oppxCorrelationDiagram_total(i).deallocate();
	    }
	  if (Input.flags.my_rank != 0)
	    oppxCorrelationDiagram_total.deallocate();
	}
      if (Input.flags.writeOppyCorrelationDiagram == 1)
	{
	  oppyCorrelationDiagram_total.allocate(1000);
	  for (int i=0; i<1000; i++)
	    {
	      oppyCorrelationDiagram_total(i).allocate(1000);
	      oppyCorrelationDiagram_total(i) = 0.0;
	      MPI_Reduce( (*QMCnode.getOppyCorrelationDiagram())(i).array(),
             oppyCorrelationDiagram_total(i).array(),1000,MPI_DOUBLE,MPI_SUM,0,
                                                               MPI_COMM_WORLD);
	      if (Input.flags.my_rank != 0)
		oppyCorrelationDiagram_total(i).deallocate();
	    }
	  if (Input.flags.my_rank != 0)
	    oppyCorrelationDiagram_total.deallocate();
	}
      if (Input.flags.writeOppzCorrelationDiagram == 1)
	{
	  oppzCorrelationDiagram_total.allocate(1000);
	  for (int i=0; i<1000; i++)
	    {
	      oppzCorrelationDiagram_total(i).allocate(1000);
	      oppzCorrelationDiagram_total(i) = 0.0;
	      MPI_Reduce( (*QMCnode.getOppzCorrelationDiagram())(i).array(),
             oppzCorrelationDiagram_total(i).array(),1000,MPI_DOUBLE,MPI_SUM,0,
                                                               MPI_COMM_WORLD);
	      if (Input.flags.my_rank != 0)
		oppzCorrelationDiagram_total(i).deallocate();
	    }
	  if (Input.flags.my_rank != 0)
	    oppzCorrelationDiagram_total.deallocate();
	}
    }

  localTimers.getGatherPropertiesStopwatch()->stop();
#else
  
  if (nalpha > 1 || nbeta > 1)
    pllSpinHistogram_total = *QMCnode.getPllSpinHistogram();
  
  if (nalpha > 0 && nbeta > 0)
    oppSpinHistogram_total = *QMCnode.getOppSpinHistogram();
  
  if (nalpha > 0)
    {
      alphaHistograms_total.allocate(nucleiTypes);
      for (int k=0; k<nucleiTypes; k++)
	alphaHistograms_total(k) = (*QMCnode.getAlphaHistograms())(k);
    }
  
  if (nbeta > 0)
    {
      betaHistograms_total.allocate(nucleiTypes);
      for (int k=0; k<nucleiTypes; k++)
	betaHistograms_total(k) = (*QMCnode.getBetaHistograms())(k);
    }
  
  if (nalpha > 1 || nbeta > 1)
    {
      if (Input.flags.writePllxCorrelationDiagram == 1)
	{
	  pllxCorrelationDiagram_total.allocate(1000);
	  for (int i=0; i<1000; i++)
	    pllxCorrelationDiagram_total(i) = 
	      (*QMCnode.getPllxCorrelationDiagram())(i);
	}
      if (Input.flags.writePllyCorrelationDiagram == 1)
	{
	  pllyCorrelationDiagram_total.allocate(1000);
	  for (int i=0; i<1000; i++)
	    pllyCorrelationDiagram_total(i) = 
	      (*QMCnode.getPllyCorrelationDiagram())(i);
	}
      if (Input.flags.writePllzCorrelationDiagram == 1)
	{
	  pllzCorrelationDiagram_total.allocate(1000);
	  for (int i=0; i<1000; i++)
	    pllzCorrelationDiagram_total(i) = 
	      (*QMCnode.getPllzCorrelationDiagram())(i);
	}
    }

  if (nalpha > 0 && nbeta > 0)
    {
      if (Input.flags.writeOppxCorrelationDiagram == 1)
	{
	  oppxCorrelationDiagram_total.allocate(1000);
	  for (int i=0; i<1000; i++)
	    oppxCorrelationDiagram_total(i) = 
	      (*QMCnode.getOppxCorrelationDiagram())(i);
	}
      if (Input.flags.writeOppyCorrelationDiagram == 1)
	{
	  oppyCorrelationDiagram_total.allocate(1000);
	  for (int i=0; i<1000; i++)
	    oppyCorrelationDiagram_total(i) = 
	      (*QMCnode.getOppyCorrelationDiagram())(i);
	}
      if (Input.flags.writeOppzCorrelationDiagram == 1)
	{
	  oppzCorrelationDiagram_total.allocate(1000);
	  for (int i=0; i<1000; i++)
	    oppzCorrelationDiagram_total(i) = 
	      (*QMCnode.getOppzCorrelationDiagram())(i);
	}
    }
#endif
}

int QMCManager::pollForACommand()
{
#ifdef PARALLEL
  // returns the tag if there is one otherwise returns -1.
  localTimers.getCommandPollingStopwatch() ->start();
  
  int poll_result = 0;
  int itag;
  MPI_Status status;
  MPI_Request request;
  
  MPI_Iprobe(  MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
               &poll_result,&status );
  
  if( poll_result==1 )
  {
    MPI_Irecv( &itag,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&request );
    
    // Wait to get the message!
    
    if(  MPI_Waitall( 1, &request, &status )  )
    {
      cerr << "ERROR: MPI_Waitall ERROR (QMCManager::MPI_poll)" << endl;
    }
  }
  else
  {
    itag = -1;
  }
  
  localTimers.getCommandPollingStopwatch() ->stop();
  return itag;
#else
  return -1;
#endif
}

void QMCManager::run()
{
  equilibrationProperties.zeroOut();
  
  done      = false;
  iteration = 0;
  
  
  if( Input.flags.equilibrate_every_opt_step > 0 )
  {
    equilibrating = true;
    Input.flags.dt = Input.flags.dt_equilibration;
  }
  
  // Create the file to write all energies out to
  ofstream *energy_strm_out = 0;
  
  if( Input.flags.write_all_energies_out == 1 )
  {
    // Create the file name
    char my_rank_str[ 32 ];
#if defined(_WIN32) && !defined(__CYGWIN__)
    _snprintf(  my_rank_str, 32, "%d", Input.flags.my_rank );
#else
    snprintf( my_rank_str, 32, "%d", Input.flags.my_rank );
#endif
    string filename = Input.flags.base_file_name + ".energy." +
      my_rank_str;
    
    // Initilize and set the output formatting
    energy_strm_out = new ofstream( filename.c_str() );
    energy_strm_out->precision( 15 );
  }
  
  // Check to make sure a valid parallelization algorithm is chosen
  if( Input.flags.parallelization_method == "pure_iterative" )
  {
    cerr << "FIX Add back in PI parallel" << endl;
  }
  else if( Input.flags.parallelization_method == "manager_worker" )
  {}
  else
  {
    cerr << "ERROR: Not a valid form of parallelization (";
    cerr << Input.flags.parallelization_method << ") in " << endl;
    cerr << "QMCManager::run()" << endl << "exiting now" << endl;
    exit( 1 );
  }
  
  while( !done )
  {
    //each process keeps track of an iteration counter
    iteration++;

    if( equilibrating )  localTimers.getEquilibrationStopwatch()->start();
    else
    {
      if ( Input.flags.use_equilibration_array == 1 )  QMCnode.startTimers();
      else localTimers.getPropagationStopwatch()->start();
    }
    
    bool writeConfigs = !equilibrating &&
      ( Input.flags.optimize_Psi || Input.flags.print_configs == 1 )  &&
      iteration%Input.flags.print_config_frequency == 0;

    while( !QMCnode.step( writeConfigs, iteration - Input.flags.equilibration_steps ) )
      {
	//Let's assume the worst; that all our data is trash.
	cerr << "Error on node " << Input.flags.my_rank << 
	  ": QMCManager is trashing data and starting calculation all over!" << endl;


	/*
	  We're starting over here, so we need to reset everything
	*/

	//Turn off whatever stopwatches were running
	if( equilibrating )  localTimers.getEquilibrationStopwatch()->stop();
	else
	  {
	    if ( Input.flags.use_equilibration_array == 1 )  QMCnode.stopTimers();
	    else localTimers.getPropagationStopwatch()->stop();
	  }

	if(  Input.flags.equilibrate_first_opt_step == 0 )
	  {
	    equilibrating = false;
	  } else {
	    equilibrating = true;
	    Input.flags.dt = Input.flags.dt_equilibration;
	  }

	//Turn on whatever stopwatches need to be running
	if( equilibrating )  localTimers.getEquilibrationStopwatch()->start();
	else
	  {
	    if ( Input.flags.use_equilibration_array == 1 )  QMCnode.startTimers();
	    else localTimers.getPropagationStopwatch()->start();
	  }
	
	done = false;
	iteration = 0;
	zeroOut();
      }
    
    if( !equilibrating )
      {
	//this is to ensure that Properties_total represents
	//the best information that we have short of a Reduction
	//Properties_total is zeroed before a reduction so this data
	//is only useful for synchronize_dmc_ensemble and the runtime output
	double weight = QMCnode.getWeights();
	weight *= QMCnode.getPopulationSizeBiasCorrectionFactor();
	Properties_total.newSample( QMCnode.getTimeStepProperties(),
				    weight,QMCnode.getNumberOfWalkers() );
      }

    if( equilibrating )
    {
      equilibrationProperties = equilibrationProperties +
      *QMCnode.getProperties();
      updateEffectiveTimeStep( &equilibrationProperties );
      
      if( Input.flags.run_type == "diffusion" )
        updateTrialEnergy( QMCnode.getWeights(),
                           Input.flags.number_of_walkers_initial );
    }
    else if( !equilibrating && Input.flags.run_type == "diffusion" )
    {
      if( Input.flags.synchronize_dmc_ensemble == 1 )
      {
        updateEstimatedEnergy( &Properties_total );
        updateEffectiveTimeStep( &Properties_total );
      }
      else
      {
        updateEstimatedEnergy( QMCnode.getProperties() );
        updateEffectiveTimeStep( QMCnode.getProperties() );
      }
      
      updateTrialEnergy( QMCnode.getWeights(),
                         Input.flags.number_of_walkers_initial );
    }
    else
      // If this is a variational calculation, the estimated energy and the
      // trial energy aren't updated.
      updateEffectiveTimeStep( QMCnode.getProperties() );
    
    if( Input.flags.checkpoint == 1 &&
        iteration%Input.flags.checkpoint_interval == 0 )
    {
      writeCheckpoint();
    }
    
    if( !equilibrating && Input.flags.write_all_energies_out == 1 )
    {
      QMCnode.writeEnergies( *energy_strm_out );
    }
    
    if( !equilibrating && Input.flags.write_electron_densities == 1 )
    {
      QMCnode.calculateElectronDensities();
    }
    
    if( Input.flags.use_hf_potential == 1 )
    {
      QMCnode.updateHFPotential();
    }
    
    if( equilibrating )
    {
      localTimers.getEquilibrationStopwatch() ->stop();
      QMCnode.zeroOut();
    }
    else
    {
      if ( Input.flags.use_equilibration_array == 1 )  QMCnode.stopTimers();
      else localTimers.getPropagationStopwatch() ->stop();
    }
    
    //--------------------------------------------------

    if(QMCManager::SIGNAL_SAYS_QUIT)
      {
	//Since we never bother turning SIGNAL_SAYS_QUIT off, we
	//need to make sure we don't write the tag more than once
	static bool written = false;
	if(!written)
	  {
	    if(QMCManager::PRINT_SIG_INFO)
	      cout << "Rank " << Input.flags.my_rank << " received SIGNAL_SAYS_QUIT."<< endl;
	    written = true;
	  }
      }
    
      if(QMCManager::INCREASE_TIME)
	{
	  QMCManager::INCREASE_TIME = false; 
	  if(QMCManager::PRINT_SIG_INFO)
	    {
	      cout << "Rank " << Input.flags.my_rank << " received INCREASE_TIME."<< endl;
	      cout << "Max time changed from " << Input.flags.max_time_steps << " to ";
	    }
	  Input.flags.max_time_steps = (long)(Input.flags.max_time_steps*1.1);
	  if(QMCManager::PRINT_SIG_INFO)
	    cout << Input.flags.max_time_steps << endl;
	}
    
    if( Input.flags.my_rank == 0 )
    {
      
      if( QMCManager::REDUCE_ALL_NOW )
	{
	  if(QMCManager::PRINT_SIG_INFO)
	    cout << "Root received REDUCE_ALL_NOW."<< endl;
	  QMCManager::REDUCE_ALL_NOW = false;

	  sendAllProcessorsACommand( QMC_REDUCE_ALL );

	  writeEnergyResultsSummary( cout );
	  writeEnergyResultsSummary( *qmcOut );

	  synchronizationBarrier();
	  gatherProperties();
	      
	  if (Input.flags.write_electron_densities == 1)
	    {
	      gatherHistograms();
	      if( Input.flags.my_rank == 0 )
		writeElectronDensityHistograms();
	    }
        
	  if ( Input.flags.calculate_bf_density == 1 )
	    gatherDensities();
	  
	  if(Input.flags.nuclear_derivatives != "none")
	    gatherForces();

	  cout << *this;
	  *getResultsOutputStream() << *this;
	  writeTimingData( cout );

	  if(Input.flags.checkpoint == 1)
	    writeCheckpoint();

	  writeRestart();
	}

      if( iteration % Input.flags.mpireduce_interval == 0 )
	{
	  sendAllProcessorsACommand( QMC_REDUCE );
	  gatherProperties();
	  if ( Input.flags.write_electron_densities == 1)
	    {
	      gatherHistograms();
	      writeElectronDensityHistograms();
	    }
	  checkTerminationCriteria();
	}
      
      if( !equilibrating && Input.flags.run_type == "diffusion" &&
          Input.flags.synchronize_dmc_ensemble == 1 &&
          iteration % Input.flags.synchronize_dmc_ensemble_interval == 0 )
      {
        sendAllProcessorsACommand( QMC_SYNCHRONIZE );
        synchronizeDMCEnsemble();
        checkTerminationCriteria();
      }
      
      if( done || QMCManager::SIGNAL_SAYS_QUIT )
      {
	done = true;
	sendAllProcessorsACommand( QMC_TERMINATE );
        
	gatherProperties();
	      
	if (Input.flags.write_electron_densities == 1)
        {
          gatherHistograms();
          writeElectronDensityHistograms();
        }
        
        if ( Input.flags.calculate_bf_density == 1 )
          gatherDensities();

        if(Input.flags.nuclear_derivatives != "none")
          gatherForces();
      }
      
      if( Input.flags.print_transient_properties &&
          iteration%Input.flags.print_transient_properties_interval == 0 )
      {
        writeTransientProperties( iteration );
      }
      
      if( iteration%Input.flags.output_interval == 0 )
      {
        writeEnergyResultsSummary( cout );
        writeEnergyResultsSummary( *qmcOut );
      }
    }
    else//else this is a worker node
    {

      int poll_result = QMC_WORK_STEP;
      do {
      
	if( iteration%Input.flags.mpipoll_interval == 0
	    || QMCManager::REDUCE_ALL_NOW
	    || QMCManager::SIGNAL_SAYS_QUIT )
	  {
	    poll_result = pollForACommand();
	  }

	switch(poll_result)
	  {
	  case QMC_REDUCE:
	    gatherProperties();
	    if ( Input.flags.write_electron_densities == 1)
	      gatherHistograms();

	    break;
	    
	  case QMC_SYNCHRONIZE:
	    synchronizeDMCEnsemble();
	    break;
	    
	  case QMC_TERMINATE:
	    done = true;
	    
	    gatherProperties();
	    
	    if (Input.flags.write_electron_densities == 1)
	      gatherHistograms();
	    
	    if ( Input.flags.calculate_bf_density == 1 )
	      gatherDensities();
	    
	    if(Input.flags.nuclear_derivatives != "none")
	      gatherForces();

	    break;	
	    
	  case QMC_REDUCE_ALL:
	    /*
	      I only want this message printed once for each
	      time the signal is sent. It seems that the rank!=0
	      processors are falling through the MPI_Reduce
	      if the root node is not at MPI_Reduce...
	    */
	    if(QMCManager::REDUCE_ALL_NOW && QMCManager::PRINT_SIG_INFO)
	      cout << "Rank " << Input.flags.my_rank << " received REDUCE_ALL_NOW."<< endl;
	    QMCManager::REDUCE_ALL_NOW = false;
	    
	    synchronizationBarrier();
	    gatherProperties();
	    
	    if (Input.flags.write_electron_densities == 1)
	      gatherHistograms();
	    
	    if ( Input.flags.calculate_bf_density == 1 )
	      gatherDensities();
	    
	    if(Input.flags.nuclear_derivatives != "none")
	      gatherForces();

	    if(Input.flags.checkpoint == 1)
	      writeCheckpoint();

	    break;
	  }	

      } while(poll_result != -1 && poll_result != QMC_WORK_STEP);
    }

    if( equilibrating )
      {
	equilibration_step();
      }
  }

  /*
    When we've gotten to this stage in the code,
    the run has completed, so let's write the most
    up to date info we have.
   */
  if( Input.flags.my_rank == 0 )
    {
      writeEnergyResultsSummary( cout );
      writeEnergyResultsSummary( *qmcOut );
    }

  if( Input.flags.checkpoint == 1 )
    {
      writeCheckpoint();
    }
  
  if(  Input.flags.write_all_energies_out == 1 )
  {
    ( *energy_strm_out ).close();
    delete energy_strm_out;
    energy_strm_out = 0;
  }
}

void QMCManager::optimize()
{
  localTimers.getOptimizationStopwatch() ->start();
  
  int configsToSkip = 0;
  
  if ( Input.flags.use_equilibration_array == 1 )
  {
    configsToSkip = 1 + (  iteration - Input.flags.equilibration_steps -
                             QMCnode.getProperties() ->energy.getNumberSamples()  )  /
    Input.flags.print_config_frequency;
  }
  
  if(  Input.flags.optimize_Psi )
  {
    QMCCorrelatedSamplingVMCOptimization::optimize( &Input,configsToSkip );
    
    // Print out the optimized parameters
    
    if(  Input.flags.my_rank == 0 )
    {
      *qmcRslts << Input.JP << endl;
    }

    //reopen the config file
    Input.openConfigFile();
  }
  
  localTimers.getOptimizationStopwatch() ->stop();
}

void QMCManager::equilibration_step()
{
  // adjust the time steps and zeros everything out during equilibration
  
  if(  Input.flags.equilibration_function == "step" )
  {
    // step function
    
    if(  iteration >= Input.flags.equilibration_steps )
    { 
      QMCnode.zeroOut();
      equilibrating = false;
      Input.flags.dt = Input.flags.dt_run;
    }
  }
  else if(  Input.flags.equilibration_function == "ramp" )
  {
    double ddt = ( Input.flags.dt_equilibration-Input.flags.dt_run ) /
    ( Input.flags.equilibration_steps-1 );
    // ramp function
    Input.flags.dt -= ddt;

    if(  iteration >= Input.flags.equilibration_steps )
    {
      Input.flags.dt = Input.flags.dt_run;
      
      QMCnode.zeroOut();
      equilibrating = false;
    }
  }
  else if(  Input.flags.equilibration_function == "CKAnnealingEquilibration1" )
  {
    if(  Input.flags.CKAnnealingEquilibration1_parameter >=
         Input.flags.equilibration_steps )
    {
      cerr << "ERROR: For CKAnnealingEquilibration1, "
      << "CKAnnealingEquilbration1_parameter must be less than "
      << "equilibration_steps!" << endl;
      exit( 0 );
    }
    
    // step function
    if(  iteration >= Input.flags.equilibration_steps )
    {
      QMCnode.zeroOut();
      equilibrating = false;
      Input.flags.dt = Input.flags.dt_run;
    }
    else
    {
      if(  iteration == 1 )
      {
        Input.flags.old_walker_acceptance_parameter += -1000;
      }
      else if(  iteration ==
                Input.flags.CKAnnealingEquilibration1_parameter )
      {
        Input.flags.old_walker_acceptance_parameter -= -1000;
        Input.flags.dt = Input.flags.dt_run;
      }
    }
  }
  else
  {
    cerr << "ERROR: Incorrect Equilibration Method Selected!" << endl;
    exit( 1 );
  }
}

void QMCManager::writeElectronDensityHistograms()
{
  // Write out the electron density histograms.  This function will only be
  // executed by the root processor.
  
#define PI 3.14159265359
  
  int nalpha = Input.WF.getNumberAlphaElectrons();
  int nbeta = Input.WF.getNumberBetaElectrons();
  
  string baseFileName = Input.flags.base_file_name;
  
  double rValue;
  double normalHist;
  double dividedHist;
  double orbital;
  double totalWeight;
  double dr = QMCnode.getdr();
  
  if (nalpha > 1 || nbeta > 1)
  {
    string pll_spin_filename = baseFileName + ".pll_pair_density";
    ofstream * pll_spin_strm = new ofstream(pll_spin_filename.c_str());
    pll_spin_strm->precision(15);
    
    // We add up the total weight of all the samples.
    totalWeight = 0.0;
    for (int i=0; i<pllSpinHistogram_total.dim1(); i++)
      totalWeight += pllSpinHistogram_total(i);
    
    *pll_spin_strm << "#\t" <<  totalWeight << endl;
    
    for (int i=0; i<pllSpinHistogram_total.dim1(); i++)
    {
      rValue = (i+0.5)*dr;
      normalHist = pllSpinHistogram_total(i)/totalWeight;
      dividedHist = normalHist/(4*PI*rValue*rValue*dr);
      orbital = sqrt(dividedHist);
      *pll_spin_strm << rValue << "\t" << pllSpinHistogram_total(i);
      *pll_spin_strm << "\t" << normalHist << "\t" << dividedHist << "\t";
      *pll_spin_strm << orbital << endl;
    }
    delete pll_spin_strm;
    pll_spin_strm = 0;
    pllSpinHistogram_total.deallocate();
  }
  
  if (nalpha > 0 && nbeta > 0)
  {
    string opp_spin_filename = baseFileName + ".opp_pair_density";
    ofstream * opp_spin_strm = new ofstream(opp_spin_filename.c_str());
    opp_spin_strm->precision(15);
    
    // We add up the total weight of all the samples.
    totalWeight = 0.0;
    for (int i=0; i<oppSpinHistogram_total.dim1(); i++)
      totalWeight += oppSpinHistogram_total(i);
    
    *opp_spin_strm << "#\t" << totalWeight << endl;
    
    for (int i=0; i<oppSpinHistogram_total.dim1(); i++)
    {
      rValue = (i+0.5)*dr;
      normalHist = oppSpinHistogram_total(i)/totalWeight;
      dividedHist = normalHist/(4*PI*rValue*rValue*dr);
      orbital = sqrt(dividedHist);
      *opp_spin_strm << rValue << "\t" << oppSpinHistogram_total(i);
      *opp_spin_strm << "\t" << normalHist << "\t" << dividedHist << "\t";
      *opp_spin_strm << orbital << endl;
    }
    delete opp_spin_strm;
    opp_spin_strm = 0;
    oppSpinHistogram_total.deallocate();
  }
  
  if (nalpha > 0 || nbeta > 0)
  {
    int nucleiTypes = Input.Molecule.NucleiTypes.dim1();
    string nucleusType;
    
    // Write out one electron densities.
    for (int i=0; i<nucleiTypes; i++)
    {
      nucleusType = Input.Molecule.NucleiTypes(i);
      
      if (nalpha > 0)
      {
        string alpha_filename = baseFileName + "." + nucleusType + 
        "-alpha.density";
        
        ofstream * alpha_strm = new ofstream(alpha_filename.c_str());
        alpha_strm->precision(15);
        
        // We add up the total weight of all the samples.
        totalWeight = 0.0;
        for (int j=0; j<alphaHistograms_total(i).dim1(); j++)
          totalWeight += (alphaHistograms_total(i))(j);
        
        *alpha_strm << "#\t" << totalWeight << endl;
        
        for (int j=0; j<alphaHistograms_total(i).dim1(); j++)
        {
          rValue = (j+0.5)*dr;
          normalHist = (alphaHistograms_total(i))(j)/totalWeight;
          dividedHist = normalHist/(4*PI*rValue*rValue*dr);
          orbital = sqrt(dividedHist);
          *alpha_strm << rValue << "\t";
          *alpha_strm << (alphaHistograms_total(i))(j) << "\t";
          *alpha_strm << normalHist << "\t" << dividedHist << "\t";
          *alpha_strm << orbital << endl;
        }
        delete alpha_strm;
        alpha_strm = 0;
	      alphaHistograms_total(i).deallocate();
      }
      
      if (nbeta > 0)
      {
        string beta_filename = baseFileName + "." + nucleusType + 
        "-beta.density";
        
        ofstream * beta_strm = new ofstream(beta_filename.c_str());
        beta_strm->precision(15);
        
        // We add up the total weight of all the samples.
        totalWeight = 0.0;
        for (int j=0; j<betaHistograms_total(i).dim1(); j++)
          totalWeight += (betaHistograms_total(i))(j);
        
        *beta_strm << "#\t" << totalWeight << endl;
        
        for (int j=0; j<betaHistograms_total(i).dim1(); j++)
        {
          rValue = (j+0.5)*dr;
          normalHist = (betaHistograms_total(i))(j)/totalWeight;
          dividedHist = normalHist/(4*PI*rValue*rValue*dr);
          orbital = sqrt(dividedHist);
          *beta_strm << rValue << "\t" << (betaHistograms_total(i))(j);
          *beta_strm << "\t" << normalHist << "\t" << dividedHist;
          *beta_strm << "\t" << orbital << endl;
        }
        delete beta_strm;
        beta_strm = 0;
	betaHistograms_total(i).deallocate();
      }
    }
    alphaHistograms_total.deallocate();
    betaHistograms_total.deallocate();
  }

  if (nalpha > 1 || nbeta > 1)
    {
      if (Input.flags.writePllxCorrelationDiagram == 1)
	{
	  string pllxFile = baseFileName + ".pllx.CorrelationDiagram";
	  ofstream * pllxStrm = new ofstream(pllxFile.c_str());
	  pllxStrm->precision(15);

	  // We just print out the raw 2D histograms

	  double min = Input.flags.pllxCorrelationDiagramMin;
	  double max = Input.flags.pllxCorrelationDiagramMax;
	  double dr = (max-min)/pllxCorrelationDiagram_total.dim1();
	  double rvalue = min;

	  for (int i=0; i<pllxCorrelationDiagram_total.dim1(); i++)
	    {
	      rvalue = min + (i+0.5)*dr;
	      *pllxStrm << "\t" << rvalue;
	    }
	  *pllxStrm << endl;

	  for (int i=0; i<pllxCorrelationDiagram_total.dim1(); i++)
	    {
	      rvalue = min + (i+0.5)*dr;
	      *pllxStrm << rvalue;
	      for (int j=0; j<(pllxCorrelationDiagram_total(i)).dim1(); j++)
		*pllxStrm << "\t" << (pllxCorrelationDiagram_total(i))(j);
	      *pllxStrm << endl;
	      pllxCorrelationDiagram_total(i).deallocate();
	    }
	  delete pllxStrm;
	  pllxStrm = 0;
	  pllxCorrelationDiagram_total.deallocate();
	}
      if (Input.flags.writePllyCorrelationDiagram == 1)
	{
	  string pllyFile = baseFileName + ".plly.CorrelationDiagram";
	  ofstream * pllyStrm = new ofstream(pllyFile.c_str());
	  pllyStrm->precision(15);

	  // We just print out the raw 2D histograms

	  double min = Input.flags.pllyCorrelationDiagramMin;
	  double max = Input.flags.pllyCorrelationDiagramMax;
	  double dr = (max-min)/pllyCorrelationDiagram_total.dim1();
	  double rvalue = min;

	  for (int i=0; i<pllyCorrelationDiagram_total.dim1(); i++)
	    {
	      rvalue = min + (i+0.5)*dr;
	      *pllyStrm << "\t" << rvalue;
	    }
	  *pllyStrm << endl;

	  for (int i=0; i<pllyCorrelationDiagram_total.dim1(); i++)
	    {
	      rvalue = min + (i+0.5)*dr;
	      *pllyStrm << rvalue;
	      for (int j=0; j<(pllyCorrelationDiagram_total(i)).dim1(); j++)
		*pllyStrm << "\t" << (pllyCorrelationDiagram_total(i))(j);
	      *pllyStrm << endl;
	      pllyCorrelationDiagram_total(i).deallocate();
	    }
	  delete pllyStrm;
	  pllyStrm = 0;
	  pllyCorrelationDiagram_total.deallocate();
	}
      if (Input.flags.writePllzCorrelationDiagram == 1)
	{
	  string pllzFile = baseFileName + ".pllz.CorrelationDiagram";
	  ofstream * pllzStrm = new ofstream(pllzFile.c_str());
	  pllzStrm->precision(15);

	  // We just print out the raw 2D histograms

	  double min = Input.flags.pllzCorrelationDiagramMin;
	  double max = Input.flags.pllzCorrelationDiagramMax;
	  double dr = (max-min)/pllzCorrelationDiagram_total.dim1();
	  double rvalue = min;

	  for (int i=0; i<pllzCorrelationDiagram_total.dim1(); i++)
	    {
	      rvalue = min + (i+0.5)*dr;
	      *pllzStrm << "\t" << rvalue;
	    }
	  *pllzStrm << endl;

	  for (int i=0; i<pllzCorrelationDiagram_total.dim1(); i++)
	    {
	      rvalue = min + (i+0.5)*dr;
	      *pllzStrm << rvalue;
	      for (int j=0; j<(pllzCorrelationDiagram_total(i)).dim1(); j++)
		*pllzStrm << "\t" << (pllzCorrelationDiagram_total(i))(j);
	      *pllzStrm << endl;
	      pllzCorrelationDiagram_total(i).deallocate();
	    }
	  delete pllzStrm;
	  pllzStrm = 0;
	  pllzCorrelationDiagram_total.deallocate();
	}
    }
  if (nalpha > 0 && nbeta > 0)
    {
      if (Input.flags.writeOppxCorrelationDiagram == 1)
	{
	  string oppxFile = baseFileName + ".oppx.CorrelationDiagram";
	  ofstream * oppxStrm = new ofstream(oppxFile.c_str());
	  oppxStrm->precision(15);

	  // We just print out the raw 2D histograms

	  double min = Input.flags.oppxCorrelationDiagramMin;
	  double max = Input.flags.oppxCorrelationDiagramMax;
	  double dr = (max-min)/oppxCorrelationDiagram_total.dim1();
	  double rvalue = min;

	  for (int i=0; i<oppxCorrelationDiagram_total.dim1(); i++)
	    {
	      rvalue = min + (i+0.5)*dr;
	      *oppxStrm << "\t" << rvalue;
	    }
	  *oppxStrm << endl;

	  for (int i=0; i<oppxCorrelationDiagram_total.dim1(); i++)
	    {
	      rvalue = min + (i+0.5)*dr;
	      *oppxStrm << rvalue;
	      for (int j=0; j<(oppxCorrelationDiagram_total(i)).dim1(); j++)
		*oppxStrm << "\t" << (oppxCorrelationDiagram_total(i))(j);
	      *oppxStrm << endl;
	      oppxCorrelationDiagram_total(i).deallocate();
	    }
	  delete oppxStrm;
	  oppxStrm = 0;
	  oppxCorrelationDiagram_total.deallocate();
	}
      if (Input.flags.writeOppyCorrelationDiagram == 1)
	{
	  string oppyFile = baseFileName + ".oppy.CorrelationDiagram";
	  ofstream * oppyStrm = new ofstream(oppyFile.c_str());
	  oppyStrm->precision(15);

	  // We just print out the raw 2D histograms

	  double min = Input.flags.oppyCorrelationDiagramMin;
	  double max = Input.flags.oppyCorrelationDiagramMax;
	  double dr = (max-min)/oppyCorrelationDiagram_total.dim1();
	  double rvalue = min;

	  for (int i=0; i<oppyCorrelationDiagram_total.dim1(); i++)
	    {
	      rvalue = min + (i+0.5)*dr;
	      *oppyStrm << "\t" << rvalue;
	    }
	  *oppyStrm << endl;

	  for (int i=0; i<oppyCorrelationDiagram_total.dim1(); i++)
	    {
	      rvalue = min + (i+0.5)*dr;
	      *oppyStrm << rvalue;
	      for (int j=0; j<(oppyCorrelationDiagram_total(i)).dim1(); j++)
		*oppyStrm << "\t" << (oppyCorrelationDiagram_total(i))(j);
	      *oppyStrm << endl;
	      oppyCorrelationDiagram_total(i).deallocate();
	    }
	  delete oppyStrm;
	  oppyStrm = 0;
	  oppyCorrelationDiagram_total.deallocate();
	}
      if (Input.flags.writeOppzCorrelationDiagram == 1)
	{
	  string oppzFile = baseFileName + ".oppz.CorrelationDiagram";
	  ofstream * oppzStrm = new ofstream(oppzFile.c_str());
	  oppzStrm->precision(15);

	  // We just print out the raw 2D histograms

	  double min = Input.flags.oppzCorrelationDiagramMin;
	  double max = Input.flags.oppzCorrelationDiagramMax;
	  double dr = (max-min)/oppzCorrelationDiagram_total.dim1();
	  double rvalue = min;

	  for (int i=0; i<oppzCorrelationDiagram_total.dim1(); i++)
	    {
	      rvalue = min + (i+0.5)*dr;
	      *oppzStrm << "\t" << rvalue;
	    }
	  *oppzStrm << endl;

	  for (int i=0; i<oppzCorrelationDiagram_total.dim1(); i++)
	    {
	      rvalue = min + (i+0.5)*dr;
	      *oppzStrm << rvalue;
	      for (int j=0; j<(oppzCorrelationDiagram_total(i)).dim1(); j++)
		*oppzStrm << "\t" << (oppzCorrelationDiagram_total(i))(j);
	      *oppzStrm << endl;
	      oppzCorrelationDiagram_total(i).deallocate();
	    }
	  delete oppzStrm;
	  oppzStrm = 0;
	  oppzCorrelationDiagram_total.deallocate();
	}
    }

#undef PI
}

void QMCManager::writeEnergyResultsSummary( ostream & strm )
{
  // Print one iteration out
  int width = 19;
  double Eave = Properties_total.energy.getAverage();
  double Estd = Properties_total.energy.getStandardDeviation();
  
  if(  Estd > 1.0e6 )
    Estd = 0.0;
  
  strm << setw( 10 )  << iteration;
  
  strm << setprecision( 10 );
  
  strm << setw(width)  << Eave << setw(width)  << Estd;
  
  strm << setw(width) << QMCnode.getNumberOfWalkers();
  
  strm << setw(width) << Input.flags.energy_trial << setw(width)  << Input.flags.dt_effective;
  
  strm << setw(width) << QMCnode.getWeights();
  
  strm << setw(width) << Properties_total.energy.getNumberSamples();
  
  strm << endl << setprecision( 15 );
}

void QMCManager::writeTransientProperties( int label )
{
  // Create the properties file name
  string filename = Input.flags.base_file_name + ".properties." +
  StringManipulation::intToString( label );
  
  // Initilize and set the output formatting
  ofstream QMCprops( filename.c_str() );
  
  QMCprops.setf( ios::fixed,ios::floatfield );
  
  QMCprops.precision( 10 );
  
  QMCprops << *this;
  
  QMCprops.close();
}

void QMCManager::writeTimingData( ostream & strm )
{
  strm << globalTimers << endl;
}

void QMCManager::writeRestart()
{
  ofstream restart( Input.flags.restart_file_name.c_str() );
  restart.setf( ios::scientific,ios::floatfield );
  restart.precision( 15 );
  
  restart << Input << endl;
  restart.close();
}

void QMCManager::writeBFDensity()
{
  ofstream density( Input.flags.density_file_name.c_str() );
  density.setf( ios::scientific,ios::floatfield );
  density.precision( 15 );
  
  for ( int i=0; i<Input.WF.getNumberBasisFunctions(); i++ )
  {
    density << Properties_total.chiDensity( i )  << endl;
  }
  
  density.close();
}

void QMCManager::writeForces()
{
  ofstream forces( Input.flags.force_file_name.c_str() );
  forces.setf( ios::scientific,ios::floatfield );
  forces.precision( 15 );
  
  for (int i=0; i<fwProperties_total.nuclearForces.dim1(); i++)
    {
      forces << fwProperties_total.nuclearForces( i )  << endl;
    }
  forces.close();
}

void QMCManager::writeXML( ostream & strm )
{
  // Write out the random seed
  ran.writeXML(strm);
  
  // Write out the number of walkers
  strm << "<NumberOfWalkers>\n" << QMCnode.getNumberOfWalkers()
  << "\n</NumberOfWalkers>" << endl;
  
  // Write out if the node is equilibrating
  strm << "<Equilibrating>\n" << equilibrating
    << "\n</Equilibrating>" << endl;
  
  // Write out the QMCrun state
  QMCnode.toXML( strm );
}

void QMCManager::writeCheckpoint()
{
  // Create the checkpoint file name
  string filename = Input.flags.checkout_file_name + ".checkpoint." +
  StringManipulation::intToString( Input.flags.my_rank );

  /*
    This is the OS dependent part of storing the checkpoint files
    in a subdirectory. Seems pretty safe... It can be commented out
    if it doesn't work.
   */
  string command = "mkdir -p " + Input.flags.checkout_file_name;
  if(!system(command.c_str()))
     filename = Input.flags.checkout_file_name + "/" + filename;

  // Initilize and set the output formatting
  ofstream QMCcheckpoint( filename.c_str() );
  //  QMCcheckpoint.setf(ios::fixed,ios::floatfield);
  QMCcheckpoint.precision( 15 );
  
  writeXML( QMCcheckpoint );
  
  QMCcheckpoint.close();
}

void QMCManager::readXML( istream & strm )
{
  string temp;
  
  // Read the random seed
  // This is the first read from the checkpoint
  strm >> temp;
  ran.readXML(strm);  
  strm >> temp;

  // Read in the number of walkers
  strm >> temp >> temp;
  Input.flags.number_of_walkers = atoi( temp.c_str() );
  strm >> temp;
  
  // Read in if the node is equilibrating
  strm >> temp >> temp;  
  equilibrating = atoi( temp.c_str() );
  strm >> temp;

  if(  !equilibrating )
  {
    Input.flags.dt = Input.flags.dt_run;
  }
  
  QMCnode.readXML( strm );
}

void QMCManager::initializeCalculationState()
{
  // Create the checkpoint file name
  string filename = Input.flags.checkin_file_name + ".checkpoint." +
  StringManipulation::intToString( Input.flags.my_rank );
  
  // open the input stream
  ifstream qmcCheckpoint( filename.c_str() );

  /*
    If no checkpoint file exists in the current directory,
    then let's look in one particular subdirectory. During
    parallel runs, each processor creates a checkpoint file,
    cluttering the directory. This might be OS dependent a bit,
    but it preserved sanity.
   */
  if(!qmcCheckpoint)
    {
      filename = Input.flags.checkin_file_name + "/" + filename;
      qmcCheckpoint.clear();
      qmcCheckpoint.open(filename.c_str());
    }

  localTimers.getInitializationStopwatch()->start();
  
  if( qmcCheckpoint && Input.flags.use_available_checkpoints == 1 )
    {
      // There is a checkpoint file
      clog << "Reading in checkpointed file " << filename << "...";
      readXML( qmcCheckpoint );
      clog << " successful." << endl;
      writeEnergyResultsSummary( clog );
      if(Input.flags.zero_out_checkpoint_statistics == 1)
	clog << "Will zero the checkpoint." << endl;
    }
  else
    {
      // There is not a checkpoint file
      QMCnode.randomlyInitializeWalkers();
    }

  localTimers.getInitializationStopwatch() ->stop();
  
  qmcCheckpoint.close();
}

void QMCManager::receiveSignal(signalType signal)
{
  /*
    We will only print a receipt for the signal if
    we are the root node since this can result in a
    lot of text for jobs with many processors.
   */
  switch(signal){
  case SIG_REDUCE:
    if(QMCManager::PRINT_SIG_INFO)
      clog << "QMCManager procedure for SIG_REDUCE: reduce all properties and dump results, then resume normally." << endl;
    QMCManager::REDUCE_ALL_NOW = true;
    break;
  case SIG_INCREASE:
    if(QMCManager::PRINT_SIG_INFO)
      clog << "QMCManager procedure for SIG_INCREASE: increase Input.flags.max_time_steps by 10%." << endl;
    QMCManager::INCREASE_TIME = true;
    break;
  case SIG_QUIT:
    if(QMCManager::PRINT_SIG_INFO)
      clog << "QMCManager procedure for SIG_QUIT: reduce and gracefully end." << endl;
    QMCManager::SIGNAL_SAYS_QUIT = true;
    break;
  case SIG_NOTHING:
    break;
  default:
    clog << "Warning: unknown signal sent." << endl;
    break;
  }
}

void QMCManager::checkTerminationCriteria()
{
  checkMaxStepsTerminationCriteria();
  checkMaxTimeTerminationCriteria();
  checkConvergenceBasedTerminationCriteria();
}

void QMCManager::checkMaxStepsTerminationCriteria()
{
  if(  Input.flags.max_time_steps == 0 )
    return;
  if(  iteration >= Input.flags.max_time_steps )
  {
    done = true;
  }
}

void QMCManager::checkMaxTimeTerminationCriteria()
{
  if(  Input.flags.max_time == 0 )
    return;
  static unsigned long total_time =
    (unsigned long)(Input.flags.max_time / Input.flags.dt_run);
  cout << "total time " << total_time << endl;
  if(  Properties_total.energy.getNumberSamples()  >=
       total_time )
  {
    done = true;
  }
}

void QMCManager::checkConvergenceBasedTerminationCriteria()
{
  if( Properties_total.energy.getNumberSamples()  > 1 &&
      Properties_total.energy.getStandardDeviation()  <=
      Input.flags.desired_convergence )
  {
    done = true;
  }
}

void QMCManager::updateEstimatedEnergy( QMCProperties* Properties )
{
  // Update the estimated energy

  if(  !equilibrating && Properties->energy.getNumberSamples()  > 1 )
  {
    Input.flags.energy_estimated = Properties->energy.getAverage();
  }
  else
  {
    // if equilibrating the estimated energy value is taken to be
    // the input energy value
  }
}

void QMCManager::updateTrialEnergy( double weights, int nwalkers_init )
{
  // Update the trial energy
  double ratio = weights / nwalkers_init;
  double logRatio = 0;
  if(ratio > 1e-250)
    {
      logRatio = log( ratio );
    } else {
      //It's pointless to keep lowering the trial energy
      logRatio = 0;
    }

  if(  equilibrating )
  {
    /*
     This section uses an estimate of the energy based on the change
     in the walkers weights.  This eliminates the upward bias introduced
     into the average energy calculated during the equilibration.  The
     estimator was derived by David R. "Chip" Kent IV and is listed
     in his notebook and possibly thesis.
     */

    //See formula 11 from UNR93
    Input.flags.energy_trial = Input.flags.energy_estimated -
      1.0 / Input.flags.dt_effective * logRatio;
  }
  else
  {
    if ( Input.flags.lock_trial_energy == 0 )
    Input.flags.energy_trial = Input.flags.energy_estimated -
      Input.flags.population_control_parameter * logRatio;
  }
}

void QMCManager::updateEffectiveTimeStep( QMCProperties* Properties )
{
  // Update the dt_effective
  QMCDerivativeProperties derivativeProperties( Properties,
  						&fwProperties_total, Input.flags.dt );
  Input.flags.dt_effective = derivativeProperties.getEffectiveTimeStep();
}

void QMCManager::synchronizationBarrier()
{
#ifdef PARALLEL
  localTimers.getCommunicationSynchronizationStopwatch() ->start();
  MPI_Barrier( MPI_COMM_WORLD );
  localTimers.getCommunicationSynchronizationStopwatch() ->stop();
#endif
}

string QMCManager::sendAllProcessorsInputFileName( char **argv )
{
#ifdef PARALLEL
  // Send the file name to be run to all the processors
  const int maximum_filename_characters = 1024;
  char c_runfilename[ maximum_filename_characters ];
  
  if(  Input.flags.my_rank == 0 )
  {
    // Copy the runfile name from the command line of the root processor
    // to the char array that will be broadcast by MPI
      strncpy( c_runfilename, argv[ 1 ], maximum_filename_characters );
  }
  
  if(  MPI_Bcast( c_runfilename, maximum_filename_characters, MPI_CHAR, 0,
                  MPI_COMM_WORLD )  )
  {
    cerr << "ERROR: Error broadcasting the runfile name to all processors"
    << endl;
    exit( 1 );
  }
  
  string runfile = c_runfilename;
  
#else
  string runfile = argv[ 1 ];
#endif
  
  return runfile;
}

QMCInput * QMCManager::getInputData()
{
  return &Input;
}

ostream * QMCManager::getResultsOutputStream()
{
  return qmcRslts;
}

void QMCManager::zeroOut()
{
  // Zero out the properties so that previous runs data aren't included
  equilibrationProperties.zeroOut();
  QMCnode.zeroOut();
  Properties_total.zeroOut();
  fwProperties_total.zeroOut();
}

ostream&  operator<<( ostream & strm, QMCManager & rhs )
{
  strm << "**************** Start of Results *****************" << endl;
  strm << rhs.Properties_total;
  strm << rhs.fwProperties_total;
  QMCDerivativeProperties derivativeProperties( &rhs.Properties_total,
						&rhs.fwProperties_total, rhs.Input.flags.dt );
  strm << derivativeProperties << endl;
  strm << "**************** End of Results *****************" << endl;
  return strm;
}
