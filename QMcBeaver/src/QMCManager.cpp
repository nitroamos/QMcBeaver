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
  
  if (Input.flags.calculate_bf_density == 1)
    {
      bool calcDensity = true;
      Properties_total.setCalcDensity(calcDensity,
                                      Input.WF.getNumberBasisFunctions());
    }
    
  QMCnode.initialize(&Input);
  
  // Set equilibrating = true so that the run equilibrates before taking
  // data
  
  if( Input.flags.equilibrate_first_opt_step == 0 )
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
      
      if(true)
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
  
  if (Input.flags.use_equilibration_array == 1)
    {
      QMCnode.stopTimers();
      *localTimers.getPropagationStopwatch() =
        *QMCnode.getPropagationStopwatch();
      *localTimers.getEquilibrationStopwatch() =
        *localTimers.getEquilibrationStopwatch() +
        *QMCnode.getEquilibrationStopwatch();
    }
    
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
  Properties_total.zeroOut();
#ifdef PARALLEL
  localTimers.getGatherPropertiesStopwatch()->start();
  
  MPI_Reduce(QMCnode.getProperties(), &Properties_total, 1,
             QMCProperties::MPI_TYPE, QMCProperties::MPI_REDUCE, 0, MPI_COMM_WORLD);
             
  localTimers.getGatherPropertiesStopwatch()->stop();
#else
  Properties_total = *QMCnode.getProperties();
#endif
}

void QMCManager::synchronizeDMCEnsemble()
{
  Properties_total.zeroOut();
#ifdef PARALLEL
  localTimers.getCommunicationSynchronizationStopwatch()->start();
  
  // Collect the global properties and send to all nodes
  
  MPI_Allreduce(QMCnode.getProperties(), &Properties_total, 1,
                QMCProperties::MPI_TYPE, QMCProperties::MPI_REDUCE, MPI_COMM_WORLD);
                
  localTimers.getCommunicationSynchronizationStopwatch()->stop();
#else
  Properties_total = *QMCnode.getProperties();
#endif
}

void QMCManager::gatherDensities()
{
#ifdef PARALLEL
  localTimers.getGatherPropertiesStopwatch()->start();
  
  Array1D<QMCProperty> localChiDensity = QMCnode.getProperties()->chiDensity;
  
  MPI_Reduce(QMCnode.getProperties()->chiDensity.array(),
             Properties_total.chiDensity.array(),Input.WF.getNumberBasisFunctions(),
             QMCProperty::MPI_TYPE,QMCProperty::MPI_REDUCE,0,MPI_COMM_WORLD);
             
  localTimers.getGatherPropertiesStopwatch()->stop();
  
#else
  for (int i=0; i<Input.WF.getNumberBasisFunctions(); i++)
    {
      Properties_total.chiDensity(i) = QMCnode.getProperties()->chiDensity(i);
    }
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
  
  done      = false;
  iteration = 0;
  
  
  if( Input.flags.equilibrate_every_opt_step > 0 )
    {
      equilibrating = true;
      Input.flags.dt = Input.flags.dt_equilibration;
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
#if defined(_WIN32) && !defined(__CYGWIN__)
      _snprintf( my_rank_str, 32, "%d", Input.flags.my_rank);
#else
snprintf( my_rank_str, 32, "%d", Input.flags.my_rank);
#endif
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
  {}
  else
    {
      cerr << "ERROR: Not a valid form of parallelization (";
      cerr << Input.flags.parallelization_method << ") in " << endl;
      cerr << "QMCManager::run()" << endl << "exiting now" << endl;
      exit(1);
    }
    
  while(!done)
    {
      //each process keeps track of an iteration counter
      iteration++;
      
      if(equilibrating) localTimers.getEquilibrationStopwatch()->start();
      else
        {
          if (Input.flags.use_equilibration_array == 1) QMCnode.startTimers();
          else localTimers.getPropagationStopwatch()->start();
        }
        
      bool writeConfigs = !equilibrating &&
              ( Input.flags.optimize_Psi || Input.flags.print_configs == 1 ) &&
                          iteration%Input.flags.print_config_frequency == 0;
                          
      QMCnode.step(writeConfigs);
      
      if( equilibrating )
        {
          equilibrationProperties = equilibrationProperties +
                                    *QMCnode.getProperties();
          updateEffectiveTimeStep(&equilibrationProperties);
          if( Input.flags.run_type == "diffusion")
            updateTrialEnergy(QMCnode.getWeights(),
                              Input.flags.number_of_walkers_initial);
        }
      else if( !equilibrating && Input.flags.run_type == "diffusion")
        {
          if( Input.flags.synchronize_dmc_ensemble == 1 )
            {
              double weight = QMCnode.getWeights();
              weight *= QMCnode.getPopulationSizeBiasCorrectionFactor();
              Properties_total.newSample(QMCnode.getTimeStepProperties(),
                                         weight,QMCnode.getNumberOfWalkers());
              updateEstimatedEnergy(&Properties_total);
              updateEffectiveTimeStep(&Properties_total);
            }
          else
            {
              updateEstimatedEnergy(QMCnode.getProperties());
              updateEffectiveTimeStep(QMCnode.getProperties());
            }
          updateTrialEnergy(QMCnode.getWeights(),
                            Input.flags.number_of_walkers_initial);
        }
      else
        // If this is a variational calculation, the estimated energy and the
        // trial energy aren't used.
        updateEffectiveTimeStep(QMCnode.getProperties());
        
      if( writeConfigs )
        {
          QMCnode.writeCorrelatedSamplingConfigurations(*config_strm_out);
        }
        
      if( Input.flags.checkpoint == 1 &&
          iteration%Input.flags.checkpoint_interval == 0 )
        {
          writeCheckpoint();
        }
        
      if( !equilibrating && Input.flags.write_all_energies_out == 1 )
        {
          QMCnode.writeEnergies(*energy_strm_out);
        }
        
      if( !equilibrating && Input.flags.write_electron_densities == 1 )
        {
          QMCnode.calculateElectronDensities();
          if (iteration%Input.flags.write_electron_densities_interval == 0)
            QMCnode.writeElectronDensityHistograms();
        }
        
      if( Input.flags.use_hf_potential == 1 )
	{
	  QMCnode.updateHFPotential();
	}

      if( equilibrating )
        {
          localTimers.getEquilibrationStopwatch()->stop();
          QMCnode.zeroOut();
        }
      else
        {
          if (Input.flags.use_equilibration_array == 1) QMCnode.stopTimers();
          else localTimers.getPropagationStopwatch()->stop();
        }
        
      //--------------------------------------------------
      
      if( Input.flags.my_rank == 0 )
        {
          if(iteration % Input.flags.mpireduce_interval == 0)
            {
              sendAllProcessorsACommand(QMC_REDUCE);
              gatherProperties();
              checkTerminationCriteria();
            }
            
          if(!equilibrating && Input.flags.run_type == "diffusion" &&
             Input.flags.synchronize_dmc_ensemble == 1 &&
             iteration % Input.flags.synchronize_dmc_ensemble_interval == 0)
            {
              sendAllProcessorsACommand(QMC_SYNCHRONIZE);
              synchronizeDMCEnsemble();
              checkTerminationCriteria();
            }
            
          if( done )
            {
              sendAllProcessorsACommand(QMC_TERMINATE);
              gatherProperties();
              if (Input.flags.calculate_bf_density == 1)
                gatherDensities();
              if (Input.flags.write_electron_densities == 1)
                QMCnode.writeElectronDensityHistograms();
            }
            
          if( Input.flags.print_transient_properties &&
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
          int poll_result = QMC_WORK_STEP;
          if(iteration%Input.flags.mpipoll_interval == 0 )
            {
              poll_result = pollForACommand();
            }
            
          if(poll_result == QMC_REDUCE)
            {
              gatherProperties();
            }
            
          else if(poll_result == QMC_TERMINATE)
            {
              done = true;
              gatherProperties();
              if (Input.flags.calculate_bf_density == 1)
                gatherDensities();
              if (Input.flags.write_electron_densities == 1)
                QMCnode.writeElectronDensityHistograms();
            }
          else if(poll_result == QMC_SYNCHRONIZE)
            {
              synchronizeDMCEnsemble();
            }
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
  
  int configsToSkip = 0;
  
  if (Input.flags.use_equilibration_array == 1)
    {
      configsToSkip = 1 + ( iteration - Input.flags.equilibration_steps -
                            QMCnode.getProperties()->energy.getNumberSamples() ) /
                      Input.flags.print_config_frequency;
    }
    
  if( Input.flags.optimize_Psi )
    {
      QMCCorrelatedSamplingVMCOptimization::optimize(&Input,configsToSkip);
      
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
  int width = 19;
  double Eave = Properties_total.energy.getAverage();
  double Estd = Properties_total.energy.getStandardDeviation();
  if( Estd > 1.0e6 )
    Estd = 0.0;
    
  strm << setw(10) << iteration;
  strm << setprecision(10);
  strm << setw(width) << Eave << setw(width) << Estd << setw(width);
  strm << Eave-Estd << setw(width) << Eave+Estd << setw(width);
  strm << QMCnode.getNumberOfWalkers() << setw(width);
  strm << Input.flags.energy_trial << setw(width) << Input.flags.dt_effective;
  strm << setw(width) << QMCnode.getWeights() << setw(width);
  strm << Properties_total.energy.getNumberSamples() << endl;
  strm << setprecision(15);
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

void QMCManager::writeBFDensity()
{
  ofstream density(Input.flags.density_file_name.c_str());
  density.setf(ios::scientific,ios::floatfield);
  density.precision(15);
  
  for (int i=0; i<Input.WF.getNumberBasisFunctions(); i++)
    {
      density << Properties_total.chiDensity(i) << endl;
    }
  density.close();
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
  
  localTimers.getInitializationStopwatch()->start();
  
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
    
  localTimers.getInitializationStopwatch()->stop();
  
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

void QMCManager::updateEstimatedEnergy(QMCProperties* Properties)
{
  // Update the estimated energy
  
  if( !equilibrating && Properties->energy.getNumberSamples() > 1 )
    {
      Input.flags.energy_estimated = Properties->energy.getAverage();
    }
  else
    {
      // if equilibrating the estimated energy value is taken to be
      // the input energy value
    }
}

void QMCManager::updateTrialEnergy(double weights, int nwalkers_init)
{
  // Update the trial energy
  if( equilibrating )
    {
      /*
      This section uses an estimate of the energy based on the change
      in the walkers weights.  This eliminates the upward bias introduced
      into the average energy calculated during the equilibration.  The
      estimator was derived by David R. "Chip" Kent IV and is listed
      in his notebook and possibly thesis.
       */
      Input.flags.energy_trial = Input.flags.energy_estimated -
               1.0 / Input.flags.dt_effective * log( weights / nwalkers_init );
    }
  else
    {
      if (Input.flags.lock_trial_energy == 0)
	Input.flags.energy_trial = Input.flags.energy_estimated -
     Input.flags.population_control_parameter * log( weights / nwalkers_init );
    }
}

void QMCManager::updateEffectiveTimeStep(QMCProperties* Properties)
{
  // Update the dt_effective
  
  QMCDerivativeProperties derivativeProperties(Properties, Input.flags.dt);
  Input.flags.dt_effective = derivativeProperties.getEffectiveTimeStep();
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
