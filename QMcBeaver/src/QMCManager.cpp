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

#include "QMCManager.h"

QMCManager::QMCManager()
{
  localTimers.getTotalTimeStopwatch()->start();
}

QMCManager::~QMCManager()
{
  finalizeOutputs();
}

void QMCManager::initialize(int argc, char **argv)
{
  initializeMPI();
  string runfile = sendAllProcessorsInputFileName(argv);

  // Load the input data
  Input.read(runfile);

  initializeOutputs();

  // Crappy random seed for parallel jobs

  Input.flags.iseed = Input.flags.iseed - Input.flags.my_rank;

  QMCnode.initialize(&Input);

  // Set equilibrating = true so that the run equilibrates before taking
  // data

  if( Input.flags.equilibrate_first_opt_step <= 0 )
    {
      equilibrating = false;
    }
  else
    {
      equilibrating = true;
    }

  done = false;


  iteration = 0;

  // Initialize the calculation state
  initializeCalculationState();
}

void QMCManager::initializeMPI()
{
#ifdef PARALLEL
  // Create MPI Community
  if( MPI_Comm_size(MPI_COMM_WORLD, &Input.flags.nprocs) )
    {
      cerr << "ERROR: MPI_Comm_size Error" << endl;
      exit(1);
    }

  // Set Processor Rank
  if( MPI_Comm_rank(MPI_COMM_WORLD, &Input.flags.my_rank) )
    {
      cerr << "ERROR: MPI_Comm_rank Error" << endl;
      exit(1);
    }
#else
  Input.flags.nprocs = 1;
  Input.flags.my_rank = 0;
#endif
}


void QMCManager::initializeOutputs()
{
  if( Input.flags.my_rank == 0 )
    {
      QMCCopyright copyright;

      // Allocate result stream
      qmcRslts = new ofstream(Input.flags.results_file_name.c_str());
      (*qmcRslts).setf(ios::fixed,ios::floatfield);
      (*qmcRslts).precision(10);
      
      if( !(*qmcRslts) )
	{
	  cerr << "ERROR opening " << Input.flags.results_file_name << endl;
	  exit(0);
	}
      
      *qmcRslts << copyright;
      
      
      // Allocate qmc output stream
      qmcOut = new ofstream(Input.flags.output_file_name.c_str());
      (*qmcOut).setf(ios::fixed,ios::floatfield);
      (*qmcOut).precision(10);
      
      if( !(*qmcOut) )
	{
	  cerr << "ERROR opening " << Input.flags.output_file_name << endl;
	  exit(0);
	}

      *qmcOut << copyright;
      
      cout.setf(ios::fixed,ios::floatfield);
      cout.precision(10);
      
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
  if( Input.flags.my_rank == 0 )
    {
      // Close and deallocate the result stream
      (*qmcRslts).close();
      delete qmcRslts;
      qmcRslts = 0;

      // Close and deallocate the qmc result stream
      (*qmcOut).close();
      delete qmcOut;
      qmcOut = 0;
    }
}

void QMCManager::finalize()
{
  //stop timing and package up the timer data
  localTimers.stop();

#ifdef PARALLEL
  MPI_Reduce(&localTimers,&globalTimers,1,QMCStopwatches::MPI_TYPE,
	     QMCStopwatches::MPI_REDUCE,0,MPI_COMM_WORLD);
#else
  globalTimers = localTimers;
#endif
}

void QMCManager::sendAllProcessorsACommand(int itag)
{
#ifdef PARALLEL
  localTimers.getSendCommandStopwatch()->start();

  MPI_Request *request = new MPI_Request[Input.flags.nprocs-1];

  // Send the command integer to each processor
  for(int i=1;i<Input.flags.nprocs;i++)
    {
      MPI_Isend(&itag,1,MPI_INT,i,itag,MPI_COMM_WORLD,&request[i-1]);
    }

  // Wait for all of the sends to be completed
  MPI_Status *status = new MPI_Status[Input.flags.nprocs-1];

  if( MPI_Waitall(Input.flags.nprocs-1, request, status) )
    {
      cerr << "ERROR: MPI_Waitall Error (QMCManager::MPI_command_send)" 
	   << endl;
    }

  delete [] request;
  delete [] status;
 
  localTimers.getSendCommandStopwatch()->stop();
#endif
}

void QMCManager::gatherProperties()
{
#ifdef PARALLEL
  localTimers.getGatherPropertiesStopwatch()->start();

  Properties_total.zeroOut();

  MPI_Reduce(QMCnode.getProperties(),&Properties_total,1,
	     QMCProperties::MPI_TYPE, QMCProperties::MPI_REDUCE,0,
	     MPI_COMM_WORLD);

  localTimers.getGatherPropertiesStopwatch()->stop();
#else
  Properties_total = *QMCnode.getProperties();
#endif
}

int QMCManager::pollForACommand()
{
#ifdef PARALLEL
  // returns the tag if there is one otherwise returns -1.
  localTimers.getCommandPollingStopwatch()->start();

  int poll_result = 0;
  int itag;
  MPI_Status status;
  MPI_Request request;

  MPI_Iprobe( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, 
	      &poll_result,&status);
  if(poll_result==1)
    {
      MPI_Irecv(&itag,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&request);

      // Wait to get the message!
      if( MPI_Waitall(1, &request, &status) )
	{
	  cerr << "ERROR: MPI_Waitall ERROR (QMCManager::MPI_poll)" << endl;
	}
    }
  else
    {
      itag = -1;
    }

  localTimers.getCommandPollingStopwatch()->stop();
  return itag;
#else
  return -1;
#endif
}


void QMCManager::run()
{
  equilibrationProperties.zeroOut();

  done           = false;
  iteration      = 0;

  
  if( Input.flags.equilibrate_every_opt_step > 0 )
    {
      equilibrating = true;
    }

  ofstream *config_strm_out = 0;
  if( Input.flags.optimize_Psi || Input.flags.print_configs == 1 )
    {
      // Initialize the output configuration stream if it is needed

      config_strm_out = new ofstream(Input.flags.config_file_name.c_str());
      (*config_strm_out).setf(ios::fixed,ios::floatfield);
      (*config_strm_out).precision(10);

      if( !(*config_strm_out) )
	{
	  cerr << "ERROR opening " << Input.flags.config_file_name << endl;
	  exit(0);
	}
    }
  else
    {
      config_strm_out = new ofstream();
    }

  // Create the file to write all energies out to
  ofstream *energy_strm_out = 0;
  if( Input.flags.write_all_energies_out == 1 )
    {
      // Create the file name
      char my_rank_str[32];
      snprintf( my_rank_str, 32, "%d", Input.flags.my_rank);
      string filename = Input.flags.base_file_name + ".energy." + 
	my_rank_str;
      
      // Initilize and set the output formatting
      energy_strm_out = new ofstream(filename.c_str());
      energy_strm_out->precision(15);
    }

  
  // Check to make sure a valid parallelization algorithm is chosen
   if( Input.flags.parallelization_method == "pure_iterative" )
    {
      cerr << "FIX Add back in PI parallel" << endl;
    }
   else if( Input.flags.parallelization_method == "manager_worker" )
     {
     }
   else
     {
       cerr << "ERROR: Not a valid form of parallelization (";
       cerr << Input.flags.parallelization_method << ") in " << endl;
       cerr << "QMCManager::run()" << endl;
       cerr << "exiting now" << endl;
      exit(1);
     }

   while(!done)
     {
       //each process keeps track of an iteration counter
       iteration++;
	  
       if(equilibrating) localTimers.getInitializationStopwatch()->start();
       else              localTimers.getPropagationStopwatch()->start();
	  
       QMCnode.step();

       updateEstimatedEnergy();
       updateTrialEnergy();
       updateEffectiveTimeStep();

       if( !equilibrating && 
	   ( Input.flags.optimize_Psi || Input.flags.print_configs == 1 ) &&
	   iteration%Input.flags.print_config_frequency == 0 )
	 {
	   QMCnode.writeCorrelatedSamplingConfigurations(*config_strm_out); 
	 }	 

       if( Input.flags.checkpoint == 1 && 
	   iteration%Input.flags.checkpoint_interval == 0 )
	 {
	   writeCheckpoint();
	 }
       
       if( Input.flags.write_all_energies_out == 1 )
	 {
	   QMCnode.writeEnergies(*energy_strm_out);
	 }

       if(equilibrating) 
	 {
	   localTimers.getInitializationStopwatch()->stop();
	   equilibrationProperties = equilibrationProperties + 
	     *QMCnode.getProperties();
	   QMCnode.zeroOut();
	 }
       else
	 {
	   localTimers.getPropagationStopwatch()->stop();
	 }

       //--------------------------------------------------

       if( Input.flags.my_rank == 0 )
	 {
	   if(iteration%Input.flags.mpireduce_interval == 0 )
	     {
	       checkTerminationCriteria();
	     }

	   if( done )
	     {
	       sendAllProcessorsACommand(TERMINATE);
	       gatherProperties();
	     }
	   else if(iteration%Input.flags.mpireduce_interval == 0 )
	     {
	       sendAllProcessorsACommand(REDUCE);
	       gatherProperties();
	     }
	      
	   if( Input.flags.my_rank == 0 &&
	       Input.flags.print_transient_properties &&
	       iteration%Input.flags.print_transient_properties_interval == 0)
	     {
	       writeTransientProperties(iteration);
	     }

	   if(iteration%Input.flags.output_interval == 0) 
	     {
	       writeEnergyResultsSummary(cout);
	       writeEnergyResultsSummary(*qmcOut);
	     }
	 }
       else
	 {
	   int poll_result = WORK_STEP;
	   if(iteration%Input.flags.mpipoll_interval == 0 )
	     {
	       poll_result = pollForACommand();
	     }

	   if(poll_result == REDUCE || poll_result == TERMINATE) 
	     {
	       gatherProperties();
	     }
	      
	   if(poll_result == TERMINATE) done = true;
	 }
	  
       if(equilibrating) 
	 {
	   equilibration_step();
	 }    
	  
     }




  if( Input.flags.optimize_Psi || Input.flags.print_configs == 1 )
    {
      (*config_strm_out).close();
    }
  delete config_strm_out;
  config_strm_out = 0;

  if( Input.flags.write_all_energies_out == 1 )
    {
      (*energy_strm_out).close();
      delete energy_strm_out;
      energy_strm_out = 0;
    }
}


void QMCManager::optimize()
{
  localTimers.getOptimizationStopwatch()->start();

  if( Input.flags.optimize_Psi )
    {
      QMCCorrelatedSamplingVMCOptimization::optimize(&Input);

      // Print out the optimized parameters
      if( Input.flags.my_rank == 0 )
	{
	  *qmcRslts << Input.JP << endl;
	}
    }

  localTimers.getOptimizationStopwatch()->stop();
}


void QMCManager::equilibration_step()
{
  // adjust the time steps and zeros everything out during equilibration
  if( Input.flags.equilibration_function == "step")
    {
      // step function
      if( iteration >= Input.flags.equilibration_steps )
	{
	  QMCnode.zeroOut();
	  equilibrating = false;
	  Input.flags.dt = Input.flags.dt_run;
	}
    }
  else if( Input.flags.equilibration_function == "ramp" )
    {
      double ddt = (Input.flags.dt_equilibration-Input.flags.dt_run)/
	(Input.flags.equilibration_steps-1);
      // ramp function
      Input.flags.dt -= ddt;
      if( Input.flags.dt <= Input.flags.dt_run  )
	{
	  Input.flags.dt = Input.flags.dt_run;
	  
	  QMCnode.zeroOut();
	  equilibrating = false;
	}
    }
  else if( Input.flags.equilibration_function == "CKAnnealingEquilibration1" )
    {
      if( Input.flags.CKAnnealingEquilibration1_parameter >= 
	  Input.flags.equilibration_steps )
	{
	  cerr << "ERROR: For CKAnnealingEquilibration1, "
	       << "CKAnnealingEquilbration1_parameter must be less than "
	       << "equilibration_steps!" << endl;
	  exit(0);
	}

      // step function
      if( iteration >= Input.flags.equilibration_steps )
	{
	  QMCnode.zeroOut();
	  equilibrating = false;
	  Input.flags.dt = Input.flags.dt_run;
	}
      else
	{
	  if( iteration == 1 )
	    {
	      Input.flags.old_walker_acceptance_parameter += -1000;
	    }
	  else if( iteration == 
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
      exit(1);
    }
}

void QMCManager::writeEnergyResultsSummary(ostream & strm)
{
  // Print one iteration out

  double Eave = Properties_total.energy.getAverage();
  double Estd = Properties_total.energy.getStandardDeviation();

  if( Estd > 1.0e6 ) Estd = 0.0;

  strm << iteration << "\t" << Eave << "\t" << Estd << "\t" 
       << Eave-Estd << "\t" << Eave+Estd << "\t" 
       << QMCnode.getNumberOfWalkers() << "\t" << Input.flags.energy_trial
       << "\t" << Input.flags.dt_effective << "\t" << QMCnode.getWeights()
       << endl;
}

void QMCManager::writeTransientProperties(int label)
{
  // Create the properties file name
  string filename = Input.flags.base_file_name + ".properties." + 
    StringManipulation::intToString(label);
  
  // Initilize and set the output formatting
  ofstream QMCprops(filename.c_str());
  QMCprops.setf(ios::fixed,ios::floatfield);
  QMCprops.precision(10);
  
  QMCprops << *this;
  
  QMCprops.close();
}

void QMCManager::writeTimingData(ostream & strm)
{
  strm << globalTimers << endl;
}

void QMCManager::writeRestart()
{
  ofstream restart(Input.flags.restart_file_name.c_str());
  restart.setf(ios::scientific,ios::floatfield);
  restart.precision(15);

  restart << Input << endl;
  restart.close();
}

void QMCManager::writeXML(ostream & strm)
{
  // Write out the random seed
  if (Input.flags.iseed > 0)
    {
      strm << "<iseed>\n" << -1*Input.flags.iseed << "\n</iseed>" << endl;
    }
  else if (Input.flags.iseed <= 0)
    {  
      strm << "<iseed>\n" << Input.flags.iseed << "\n</iseed>" << endl;
    }

  // Write out the number of walkers
  strm << "<NumberOfWalkers>\n" << QMCnode.getNumberOfWalkers()
       << "\n</NumberOfWalkers>" << endl;

  // Write out if the node is equilibrating
  strm << "<Equilibrating>\n" << equilibrating 
       << "\n</Equilibrating>" << endl;

  // Write out the QMCrun state
  QMCnode.toXML(strm);
}

void QMCManager::writeCheckpoint()
{
  // Create the checkpoint file name
  string filename = Input.flags.base_file_name + ".checkpoint." + 
    StringManipulation::intToString(Input.flags.my_rank);
  
  // Initilize and set the output formatting
  ofstream QMCcheckpoint(filename.c_str());
  //  QMCcheckpoint.setf(ios::fixed,ios::floatfield);
  QMCcheckpoint.precision(15);

  writeXML(QMCcheckpoint);

  QMCcheckpoint.close();
}

void QMCManager::readXML(istream & strm)
{  
  string temp;

  // Read the random seed
  strm >> temp;
  strm >> temp;
  Input.flags.iseed = atoi(temp.c_str());
  if (Input.flags.iseed > 0) Input.flags.iseed *= -1;
  strm >> temp;

  // Read in the number of walkers
  strm >> temp >> temp;
  Input.flags.number_of_walkers = atoi(temp.c_str());
  strm >> temp;

  // Read in if the node is equilibrating
  strm >> temp >> temp;
  equilibrating = atoi(temp.c_str());
  strm >> temp;

  if( !equilibrating )
    {
      Input.flags.dt = Input.flags.dt_run;
    }

  QMCnode.readXML(strm);
}

void QMCManager::initializeCalculationState()
{
  // Create the checkpoint file name
  string filename = Input.flags.base_file_name + ".checkpoint." + 
    StringManipulation::intToString(Input.flags.my_rank);

   // open the input stream
  ifstream qmcCheckpoint(filename.c_str());

  if( qmcCheckpoint && Input.flags.use_available_checkpoints == 1 )
    {
      // There is a checkpoint file
      readXML(qmcCheckpoint);
    } 
  else
    {
      // There is not a checkpoint file
      QMCnode.randomlyInitializeWalkers();
    }

  qmcCheckpoint.close();
}

void QMCManager::checkTerminationCriteria()
{
  checkMaxStepsTerminationCriteria();
  checkConvergenceBasedTerminationCriteria();
}

void QMCManager::checkMaxStepsTerminationCriteria()
{  
  if( Properties_total.energy.getNumberSamples() >= 
      Input.flags.max_time_steps )
    {
      done = true;
    }
}

void QMCManager::checkConvergenceBasedTerminationCriteria()
{
  if(Properties_total.energy.getNumberSamples() > 1 && 
     Properties_total.energy.getStandardDeviation() <= 
     Input.flags.desired_convergence )
    {
      done = true;
    }
}

void QMCManager::updateEstimatedEnergy()
{
  // Update the estimated energy

  if( !equilibrating && 
      QMCnode.getProperties()->energy.getNumberSamples() > 1 )
    {
      Input.flags.energy_estimated = 
	QMCnode.getProperties()->energy.getAverage();
    }
  else
    {
      // if equilibrating the estimated energy value is taken to be
      // the input energy value
    }
}

void QMCManager::updateTrialEnergy()
{
  // Update the trial energy

  if( !equilibrating && 
      QMCnode.getProperties()->energy.getNumberSamples() > 1 )
    {
      // determine the total of all the walker weights
      double total_weights = QMCnode.getWeights();

      Input.flags.energy_trial = Input.flags.energy_estimated - 
	Input.flags.population_control_parameter * 
	log( total_weights / Input.flags.number_of_walkers_initial );
    }
  else
    {
      /*
	This section uses an estimate of the energy based on the change
	in the walkers weights.  This eleminates the upward bias introduced
	into the average energy calculated during the equilibration.  The
	estimator was derived by David R. "Chip" Kent IV and is listed
	in his notebook and possibly thesis.
       */

      static const double originalEest = Input.flags.energy_estimated;

      // determine the total of all the walker weights
      double total_weights = QMCnode.getWeights();

      Input.flags.energy_trial = originalEest - 
	1.0 / Input.flags.dt_effective * 
	log( total_weights / Input.flags.number_of_walkers_initial );
    }
}

void QMCManager::updateEffectiveTimeStep()
{
  // Update the dt_effective

  if( !equilibrating && 
      QMCnode.getProperties()->energy.getNumberSamples() > 1 )
    {
      QMCDerivativeProperties derivativeProperties( QMCnode.getProperties(), 
						    Input.flags.dt);

      Input.flags.dt_effective = derivativeProperties.getEffectiveTimeStep();
    }
  else
    {
      if( equilibrationProperties.energy.getNumberSamples() > 1 )
	{
	  QMCDerivativeProperties derivativeProperties( 
			     &equilibrationProperties, Input.flags.dt);

	  Input.flags.dt_effective = 
	    derivativeProperties.getEffectiveTimeStep();
	}
    }
}

void QMCManager::synchronizationBarrier()
{
#ifdef PARALLEL
  localTimers.getCommunicationSynchronizationStopwatch()->start();

  MPI_Barrier(MPI_COMM_WORLD);

  localTimers.getCommunicationSynchronizationStopwatch()->stop();
#endif
}


string QMCManager::sendAllProcessorsInputFileName(char **argv)
{
#ifdef PARALLEL
  // Send the file name to be run to all the processors
  const int maximum_filename_characters = 100;
  char c_runfilename[maximum_filename_characters];
  
  if( Input.flags.my_rank == 0 )
    {
      // Copy the runfile name from the command line of the root processor
      // to the char array that will be broadcast by MPI
      strcpy(c_runfilename,argv[1]);
    }

  if( MPI_Bcast(c_runfilename, maximum_filename_characters, MPI_CHAR, 0,
		MPI_COMM_WORLD) )
    {
      cerr << "ERROR: Error broadcasting the runfile name to all processors" 
	   << endl;
      exit(1);
    }  
  
  string runfile = c_runfilename;
  
#else
  string runfile = argv[1];
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
}

ostream&  operator<<(ostream & strm, QMCManager & rhs)
{
  strm << "**************** Start of Results *****************" << endl;
  
  strm << rhs.Properties_total;
  
  QMCDerivativeProperties derivativeProperties(&rhs.Properties_total,
					       rhs.Input.flags.dt);

  strm << derivativeProperties << endl;

  strm << "**************** End of Results *****************" << endl;

  return strm;
}
