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

#ifndef QMCFLAGS_H
#define QMCFLAGS_H

#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h>

using namespace std;

/**
   Information specifying how to perform this QMC calculation.
*/

class QMCFlags
{
 public:
  /**
     Creates an instance of this class.
  */
  QMCFlags();

  /**
     Base file name for the calculation.  This is the file name without
     any extension.
  */
  string base_file_name;

  /**
     Input file name for the calculation.  This is the base file 
     name with the ".ckmf" extension.
  */
  string input_file_name;

  /**
     Output file name for the calculation.  This is the base file 
     name with the ".qmc" extension.
  */
  string output_file_name;

  /**
     Restart input file name for the calculation.  This is the base file 
     name with the ".01.ckmf" extension.
  */
  string restart_file_name;

  /**
     Output results file name for the calculation.  This is the base file 
     name with the ".rslts" extension.
  */
  string results_file_name;

  /**
     Saved walker file name for the calculation.  This is the base file 
     name with the ".<processor>.cfgs" extension.
  */
  string config_file_name;

  /**
     Method for initializing the walkers for the calculation.  
     @see QMCInitializeWalkerFactory
  */
  string walker_initialization_method;

  
  /**
     Number of walkers to generate when initializing the walkers.  The
     desired number of walkers will be chosen from the generated walkers.
  */
  int walker_initialization_combinations;

  /**
     Type of QMC calculation.  For exampe: variational, diffusion, etc.
  */
  string run_type;

  /**
     Scratch directory for the calculation.  Temporary files which will
     not be needed later will be written here.
  */
  string temp_dir;

  /**
     Parallelization algorithm to use.
  */
  string parallelization_method;

  /**
     Seed used in generating pseudo-random numbers for the calculation.
  */
  long iseed;

  /**
     Diffusion Green's function to be used for the calculation.  
     Only importance sampled Green's functions should be used with DMC.
  */
  string sampling_method;

  /**
     Should the "quantum force" for the calculation be modified.
     Modifications exist which, for example, can reduce the time-step
     bias in DMC.
  */
  string QF_modification_type;

  /**
     Should the "local energy" for the calculation be modified.
     Modifications exist which, for example, can reduce the time-step
     bias in DMC.
  */
  string energy_modification_type;

  /**
     If Umrigar's 1993 QMC algorithm (J. Chem. Phys. 99, 2865(1993))
     is used with the equal-electron modifications to the "local energy"
     and "quantum force," this is the parameter describing how the
     modifications will occur.
  */
  double umrigar93_equalelectrons_parameter;

  /**
     Algorithm used to reweight the walkers after a time step.
  */
  string walker_reweighting_method;

  /**
     Method for branching walkers.  Possible methods include not
     branching, branching with all walkers getting unit weight, 
     branching with walkers getting fractional weights, and similar 
     methods.
  */
  string branching_method;

  /**
     If a walker's weight exceeds this threshold, it will branch
     when the branching method allows the walkers to have fractional
     weights.
  */
  double branching_threshold;

  /**
     If a walker's weight is less than this threshold, it will fuse
     with another low weight walker when the branching method allows 
     the walkers to have fractional weights.
  */
  double fusion_threshold;

  /**
     Parameter which can be used to prevent persistent configurations
     from interfering with a calculation.  A description of this is
     found in J. Chem. Phys. 99, 2865(1993).
  */
  int old_walker_acceptance_parameter;

  /**
     Time step being used for the calculation.
  */
  double dt;

  /**
     Time step to be used for the non-equilibration portion of the
     calculation.
  */
  double dt_run;

  /**
     Last updated value of the effective time step.
  */
  double dt_effective;

  /**
     Calculation terminates if the standard deviation of the energy
     falls below this threshold.
  */
  double desired_convergence;

  /**
     Maximum number of time steps the calculation will run for.
  */
  long max_time_steps;

  /**
     Current number of walkers.
  */
  long number_of_walkers;

  /**
     Initial number of walkers.  Branching calculations will try to 
     maintain this many walkers.
  */
  long number_of_walkers_initial;

  /**
     Use interpolation of the radial basis functions in an attempt
     to speed up the calculation.
  */
  int use_basis_function_interpolation;

  /**
     Number of grid points to use in interpolating the radial basis
     functions.
  */
  int number_basis_function_interpolation_grid_points;

  /**
     Smallest point to use in interpolating the radial basis functions.
  */
  double basis_function_interpolation_first_point;

  /**
     Number of time steps taken between producing output.
  */
  int output_interval;

  /**
     1 if the accumulated statistics in the input checkpoint file
     are to be discarded, and 0 otherwise.
  */
  int zero_out_checkpoint_statistics;

  /**
     Number of time steps to take for the calculation to equilibrate.
  */
  unsigned int equilibration_steps;

  /**
     Initial time step used during the equilibration of the calculation.
  */
  double dt_equilibration;

  /**
     Procedure used to modify the time step while the calculation is
     equilibrating.
  */
  string equilibration_function;

  /**
     Number of time steps to use the CKAnnealingEquilibration1 before
     just equilibrating with the time step of the non-equilibration
     portion of the calculation.  The CKAnnealingEquilibration1 
     algorithm accepts essentially any proposed move.  If importance
     sampling is used, this can rapidly reach a "good" section of
     configuration space.
  */
  unsigned int CKAnnealingEquilibration1_parameter;

  /**
     Number of time steps taken on the root processor before the
     results are collected from all other processors.
  */
  int mpireduce_interval;

  /**
     Number of time steps taken on the worker processors before
     checking for commands from the root processor.
  */
  int mpipoll_interval;

  /**
     Number of time steps taken before this processor writes a 
     checkpoint file.
  */
  int checkpoint_interval;

  /**
     1 if the calculation is to checkpoint, and 0 otherwise.
  */
  int checkpoint;

  /**
     Use existing checkpoints when performing this calculation.
  */
  int use_available_checkpoints;

  /**
     Number of time steps between printing the global properties as
     the calculation progresses.
  */
  int print_transient_properties_interval;

  /**
     1 if the global properties of the calculation are output as the
     calculation progresses, and 0 otherwise.
  */
  int print_transient_properties;

  /**
     Number of atoms in the system.
  */
  int Natoms;

  /**
     Charge of the system in atomic units.  (e.g. +1 if a neutral system
     looses an electron)
  */
  int charge;

  /**
     Number of orbitals in the trial wavefunction.
  */
  int Norbitals;

  /**
     Number of basis functions used to create the orbitals for the
     trial wavefunction.
  */
  int Nbasisfunc;

  /**
     Trial energy used for branching QMC calculations.
  */
  double energy_trial;

  /**
     Last estimate of the energy of the system.
  */
  double energy_estimated;

  /**
     Original estimate of the energy of the system read from the 
     input file.
  */
  double energy_estimated_original;

  /**
     1 if the correction for the population size bias is used, and 0
     otherwise.  See J. Chem. Phys. 99, 2865(1993).
  */
  int correct_population_size_bias;

  /**
     1 if walkers from the calculation are written out to file, and 0
     otherwise.
  */
  int print_configs;

  /**
     Number of time steps between writing walkers out to file.
  */
  int print_config_frequency;

  /**
     1 if the wavefunction is to be optimized during a VMC wavefunction,
     and 0 otherwise.
  */
  int optimize_Psi;
  
  /**
     Number of vmc-optimize iterations used when optimizing the 
     wavefunction.
  */
  int max_optimize_Psi_steps;

  /**
     Objective function parameter used in forming a barrier around 
     a valid region of parameter space.  Not all objective functions
     use such a barrier.  A valid region of parameter space is one where
     the wavefunction is not majorly different from the one used
     to generate the walkers when correlated sampling is used.
  */
  double optimize_Psi_barrier_parameter;

  /**
     Objective function to optimize when optimizing a VMC wavefunction.
  */
  string optimize_Psi_criteria;

  /**
     Numerical optimization algorithm to use when optimizing a VMC
     wavefunction.
     @see QMCOptimizationFactory
  */
  string optimize_Psi_method;

  /**
     Objective function to use when calculating derivatives for
     optimizing a VMC wavefunction.
  */
  string numerical_derivative_surface;

  /**
     Maximum number of iterations to use when optimizing a VMC 
     wavefunction with one set of correlated-sampling configurations.
  */
  int optimization_max_iterations;

  /**
     Tolerance used to determine when a numerical optimization 
     with one set of correlated-sampling configurations has converged.
  */
  double optimization_error_tolerance;

  /**
     Population size used when the CKGeneticAlgorithm1 is used to
     optimize a VMC wavefunction.
     @see CKGeneticAlgorithm1
  */
  int ck_genetic_algorithm_1_population_size;

  /**
     Mutation rate used when the CKGeneticAlgorithm1 is used to
     optimize a VMC wavefunction.
     @see CKGeneticAlgorithm1
  */
  double ck_genetic_algorithm_1_mutation_rate;

  /**
     Standard deviation used when initializing the CKGeneticAlgorithm1.
     @see CKGeneticAlgorithm1
  */
  double ck_genetic_algorithm_1_initial_distribution_deviation;

  /**
     Algorithm used to determine the step length used when a line
     search algorithm is used to optimize a VMC wavefunction.
  */
  string line_search_step_length;
 
  /**
     Parameter used for the logarithmic barrier when D.R. Kent IV's 
     algorithm for optimizing possibly singular Jastrow functions is 
     used.  See David Randall Kent IV's thesis, New Quantum Monte Carlo 
     Algorithms to Efficiently Utilize Massively Parallel Computers, 
     from the California Institute of Technology for more details.  
     This parameter should be a small positive number.
  */
  double singularity_penalty_function_parameter;

  /**
     1 if up and down electrons use the same Jastrow parameters,
     and 0 otherwise.
  */
  int link_Jastrow_parameters;

  /**
     1 if every VMC calculation used to generated correlated-sampling
     configurations for an optimization equilibrates for the specified 
     equilibration time, and 0 otherwise.
  */
  int equilibrate_every_opt_step;

  /**
     1 if the first VMC calculation used to generated correlated-sampling
     configurations for an optimization equilibrates for the specified 
     equilibration time, and 0 otherwise.  If a non-VMC calculation is
     performed, 1 if the calculation equilibrates, and 0 otherwise.
  */
  int equilibrate_first_opt_step;

  /**
     Parameter used in controlling the number of walkers or total 
     weight of all walkers during a branching QMC calculation.
  */
  double population_control_parameter;

  /**
     1 if the local energy for every walker for every time step is
     to be written out, and 0 otherwise.
  */
  int write_all_energies_out;

  /**
     Are Chip and Mike cool?  Answer: Yea Baby!
  */
  string chip_and_mike_are_cool;
  
  /**
     MPI rank of this processor.
  */
  int my_rank;

  /**
     Number of processors used in this calculation.
  */
  int nprocs;

  /**
     Load this object's state from a QMC input file and initialize 
     the object.

     @param InFileName QMC input file to load this object's state from.
  */
  void read_flags(string InFileName);

  /**
     From the input file name, generate other file names that will be
     used during the calculation. (e.g. base_file_name)

     @param runfile QMC input file name
  */
  void set_filenames(string runfile);

  /**
     Write this object's state out to a stream. The same format is 
     used as in the QMC input file.
  */  
  friend ostream& operator <<(ostream& strm, QMCFlags& flags);
};

#endif


















