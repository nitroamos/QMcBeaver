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
#include <vector>
#include <climits>
#include <algorithm>

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
    This facilitates easy testing of parameters through the input
    file. Enter as many as desired, separated by a space.
  */
  vector<long> programmersLongs;
  
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
     Name of the checkpoint file to read in.
     The checkpoint files will automatically have a suffix added to them
     according to parallel node rank. It will search both the current
     directory for the checkpoint file, as well as in a directory also
     called checkin_file_name.

     There is some use in permitting a different input than output since
     several runs might want to use the same input, and you don't want to
     overwrite the original checkpoints. For example, the checkpoint represents
     a (probably) equilibrated distribution.

     If this flag is not specified, it will assume the name of the
     ckmf file.
  */
  string checkin_file_name;

  /**
     Name of the checkpoint file to write out. It will try to put
     the checkpoint files into a directory with the name
     checkout_file_name.

     If this flag is not specified, it will assume the name of the
     ckmf file.
  */
  string checkout_file_name;

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
     Basis function density file name for the calculation.  This is the base 
     file name with the ".density" extension.
  */
  string density_file_name;

  /**
     Storage for the calculated forces.  This is the base 
     file name with the ".force" extension.
  */
  string force_file_name;

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
     Whether to update one electron at a time (1) or all at once (0).

     The discussions that I've found on the issue are 
     pg 97, pgs 273-278 of the Hammond, Lester, Reynolds QMC book
     They say acceptance probability is "reduced considerably" if all
     electrons are moved at once.

     pg 144 of the Umrigar "VMC basics and applications to atoms and molecules" 
     from the NATO series, 1998...
     He also recommends moving one electron at a time.

     On the other hand, if we move one electron per iteration,
     then the amount of system memory required scales with the
     number of walkers.

     You will want to link LAPACK if you move all at once. If you
     move one at a time, we use a different method to update the inverse.
  */
  int one_e_per_iter;

  /**
     Using QMCSurfer instead of QMCWalker. This will let you "surf"
     the wavefunction.

  */
  int use_surfer;

  /**
     Minimize the amount of output.
  */
  int set_debug;

  /**
    These are debugging flags. They might be used in QMCSurfer.
  */
  int use_jastrow;
  int detailed_energies;

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
     Implemented here is the method described in:
     S. Chiesa, D.M. Ceperley, and S. Zhang, Phys. Rev. Lett. 94, 036404 (2005)
     This selects the method used in calculating the forces in a molecule.
     Choices are none, bare_hellmann_feynman, poly_hellmann_feynman,
     and slaterpoly_hellmann_feynman.
  */
  string nuclear_derivatives;
  
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
     Trail's eps^2 value for his weighting method found in
     Alternative sampling for variational quantum Monte Carlo
     Phys Rev E 77, 016704 (2008)

     He recommends that eps2 = variance
  */
  double trail_eps2;

  /**
     The e_0 value used in:
     Alternative sampling for variational quantum Monte Carlo
     Phys Rev E 77, 016704 (2008)

     It is supposed to hold the best known VMC energy.
  */
  double e_0;

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
     1 if the global results are to be collected and broadcast to all nodes and
     used to calculate the estimated energy, trial energy, and effective time
     step.
  */
  int synchronize_dmc_ensemble;

  /**
     The interval at which the DMC ensemble is to be synchronized.
  */
  int synchronize_dmc_ensemble_interval;

  /**
     Parameter which can be used to prevent persistent configurations
     from interfering with a calculation.  A description of this is
     found in J. Chem. Phys. 99, 2865(1993).
  */
  int old_walker_acceptance_parameter;

  /**
     Warning verbosity level. A higher integer should correspond to more
     warnings printed.
   */
  int warn_verbosity;

  /**
     Controls
   */
  double rel_cutoff;

  /**
     Branching recommendations
   */
  int limit_branching;
  
  /**
     Time step being used for the calculation.

     Notes: There is sort of a balance here. If dt
     is large, then your error will be large. If dt
     is small, it will take longer to explore all the
     space.

     Officially, you're supposed to estimate an extrapolation
     to dt = 0. That is, run it for (for example) 1e-3, 
     then 1e-4, then 5e-5, etc until you can guess what
     the value at zero would be.

     (AGA note, maybe allow this to be input as an array, and run
     all the calculations simultaneously or in succession...)
  */
  double dt;

  /**
     These two parameters are used by Umrigar's metropolis acceleration
     scheme, labeled here umrigar93_accelerated_sampling.
     See QMCWalker for details.
  */
  double accel_delta;
  double accel_tm;

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

     Notes: It's difficult to determine a priori how many time
     steps it will take to reach a given convergence, so it's often
     easier to control the length of the run by setting this to 0.0
     and choosing an appropriate max_time_steps.
  */
  double desired_convergence;

  /**
     Maximum number of time steps the calculation will run for.
     Set max_time_steps = 0 if you do not want it to be used as termination
     criteria.
  */
  unsigned long max_time_steps;
  unsigned long original_max_time_steps;

  /**
     Maximum amount of time the calculation will run for. The number
     of steps that this will limit to will be max_time / dt. The main
     point of this parameter is to have one that will change automatically
     depending on dt. This makes it easier to use the same input files
     for different dt.

     Set max_time = 0 if you do not want it to be used as termination
     criteria.

     There is an ambiguitiy as far as what this parameter means when doing
     a parallel calculation. Should this limit the amount of time that
     each processor individually covers or the "collective" time that all
     processors cover? I've chosen the later since the former would involve
     a loss of parallelization efficiency.
  */
  double max_time;

  /**
     This parameter is used in QMCRun to process several walkers 
     simultaneously. The value of this macro is how many walkers to treat at 
     once. If this parameter is set to 1, then the code will perform the way it
     did previously. Also, the amount of memory the program requires is 
     dependent on this because all the walkers have to be stored.
    
     This used instead of the WALKERS_PER_PASS macro.
     DON"T CHANGE THIS IN THE MIDDLE OF A CALCULATION!
  */
  long walkers_per_pass;

  /**
    Implemented here is the future walking method described in:
    J. Casulleras and J. Boronat, Phys. Rev. B 52, 3654 (1995)

    This vector represents future walking projection time.
    Inspired by the decorrelation algorithm,
    it might be best to set this to something exponential like
    1.0 2.0 4.0 8.0 16.0 32.0

    A 0 is automatically included in the set since 0 represents
    no future walking. Preliminary results seem to indicate
    that future walking converges somewhere between
    1 and 10 Hartrees^-1. For smaller dt, this requires more
    iterations.

    In memory, we actually store it as block iteration length,
    since this is what QMCWalker uses.

    Block Length = Time / dt    
  */
  vector<int> future_walking;

  long gpu_walkers_per_pass;

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
    This is the index of the Lebedev-Laikov grid to use for the run
    on each nucleus.
    6  pts for pseudo_gridLevel = 0
    14 pts for pseudo_gridLevel = 1
    26 pts for pseudo_gridLevel = 2
    38 pts for pseudo_gridLevel = 3
    50 pts for pseudo_gridLevel = 4
    etc
    see QMCMolecule::readPseudoPotential for the other sizes

    There are 32 levels to choose from but 14 or 26 points should
    be sufficient.
  */
  int pseudo_gridLevel;
  
  /**
    If the local part of the pseudopotential is sufficiently small,
    then we can forgo the calculation of the nonlocal part.
  */
  double pseudo_cutoff;

  /**
    This is an internal flag, not an input flag.
    It is turned on if we read any ECPs.
    It helps control memory allocation.
  */
  int use_pseudopotential;

  /**
     Number of time steps taken between producing output.
  */
  int output_interval;

  /**
     1 if the accumulated statistics in the input checkpoint file
     are to be discarded, and 0 otherwise.

     This is useful if you want to keep the (equilibrated) walker
     positions and weights, but don't want the accumulated statistics.
  */
  int zero_out_checkpoint_statistics;

  /**
     Number of time steps to take for the calculation to equilibrate.
  */
  unsigned long equilibration_steps;

  /**
     Initial time step used during the equilibration of the calculation.

     Notes: You might want to set this higher than you would set your
     regular dt since the point is to get the electrons in approximately
     the right place, not to make measurements.
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
     1 if QMCEquilibrationArray is to be used.  This defines an array of 
     QMCProperties objects, where the ith element starts collecting statistics
     on the (2^i)th iteration.  The element with the lowest standard deviation
     for the total energy is chosen, automatically determining the correct
     number of equilibrtion steps.
     0 if all the statistics are to be collected into one QMCProperties object.
  */

  int use_equilibration_array;

  /**
     Number of time steps taken on the root processor before the
     results are collected from all other processors.
     
     Notes: the longer this is, the longer it will take for the
     calculation to converge. The short this is, the more time
     the calculation will spend shuttling data back and forth.

     This balance might be molecule specific since the communication
     time depends on how much data to send.
  */
  int mpireduce_interval;

  /**
     Number of time steps taken on the worker processors before
     checking for commands from the root processor.

     Notes: if this is too large, then the root processor might
     end of spending a long time waiting. There isn't
     a significant penalty to making this small. Perhaps this
     is the sort of parameter we want to remain opaque to the user.

     It should probably be pretty small relative to mpireduce_interval.
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
     To make them smaller, write only the energy data to the checkpoint files.
  */
  int checkpoint_energy_only;

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
     Number of determinants in the SCF part of the wavefunction.
  */
  int Ndeterminants;

  /**
    Restricted or unrestricted trial functions can be used.
  */
  string trial_function_type;

  /**
     1 if the basis function density is to be calculated, 0 if not.
  */
  int calculate_bf_density;

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

     Notes: turning this on is required for wavefunction optimization.
     However, it can *easily* take up all of your compute time since it
     obviously requires harddrive access. If you select optimization,
     this will automatically be set to 1.

     Also, depending on print_config_frequency and molecule size, these
     files can easily get up to be several gigabytes in size. The files
     can be output in ASCII or binary -- that option is set in QMCConfigIO.
     
     Eventually, we'll probably add HDF5 functionality, but last I (AGA) checked
     they were still working on their C++ API.
  */
  int print_configs;

  /**
     Number of time steps between writing walkers out to file.
  */
  int print_config_frequency;

  /**
     1 if the wavefunction is to be optimized during a VMC wavefunction,
     and 0 otherwise.

     Currently, we can optimize:
     1) EE Jastrow parameters
     2) NE Jastrow parameters
     3) CI coefficients

     Notes: This will rerun the entire VMC calculation each time it
     optimizes. That is, your calculation time might go up by a factor
     of max_optimize_Psi_steps, depending on convergence.
  */
  int optimize_Psi;
  
  /**
     Number of vmc-optimize iterations used when optimizing the 
     wavefunction.

     Notes: Wavefunction convergence can take a long time, depending
     on your Jastrow parameters.
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
     This is an internal flag. Some flags will require us to
     calculate the analytic wavefunction parameter derivatives,
     others won't. QMCFlags will figure it out and set this parameter
     accordingly.
  */
  int calculate_Derivatives;

  /**
    These next flags control which part of the wavefunction
    we're optimizing. If the optimization method is set to
    automatic, then these parameters will automatically be changed.
  */
  int optimize_EE_Jastrows;
  int optimize_EN_Jastrows;
  int optimize_NEE_Jastrows;
  int optimize_CI;
  int optimize_Orbitals;

  /**
    Some Jastrows (e.g. the cambridge ones) use cutoffs which can be difficult
    to optimize. this one parameter will let you decide whether to optimize them.
  */
  int optimize_L;

  int use_three_body_jastrow;

  /**
    If this flag is true, we allow the terms in the three body Jastrow that
    overlap the NE Jastrow to be nonzero.
  */
  int reproduce_NE_with_NEE_jastrow;

  /**
    If this flag is true, we allow the terms in the three body Jastrow that
    overlap the EE Jastrow to be nonzero.
  */
  int reproduce_EE_with_NEE_jastrow;

  /**
    This controls the initial value used by QMCLineSearch for
    modifying the hessian. If this is large, then the optimization
    will be more like steepest descent, which is slower but more
    stable.

    If it's smaller, then optimization should happen faster.
  */
  double a_diag;

  /**
    A parameter used by the Linearize step length method, must
    be between 0 and 1.

    The method is most stable when ksi = 0.
  */
  double ksi;

  /**
     Objective function to optimize when optimizing a VMC wavefunction.
  */
  string optimize_Psi_criteria;

  /**
     Numerical optimization algorithm to use when optimizing a VMC
     wavefunction.

     Notes: for QMC, if your wavefunction is an eigenfunction, your
     variance will be exactly zero. That's why you might want to minimize
     your variance.

     If you set this to "automatic", then the code will attempt to
     select the parameter choices that probably work best. However,
     you might want to tweek some stuff anyway. For example, the
     a_diag parameter in QMCLineSearch.cpp
     the run lengths in QMcBeaver.cpp
     
     @see QMCOptimizationFactory
     @see QMCLineSearch
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

     Notes: This is how many iterations the optimization might take to
     converge a set of parameters over your wavefunction during a
     given optimization cycle.
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
     1 if up and down electrons use the same Jastrow parameters for
     their electron-nucleus terms, and 0 otherwise.
  */
  int link_Jastrow_parameters;
  int link_NEE_Jastrows;

  /**
     1 if up and down electrons use the same orbital parameters.
     
     Just because the input file used the same orbitals from the
     ckmf file doesn't mean they can't be optimized independently.

     The restart file will have more orbitals than the input file.
  */
  int link_Orbital_parameters;
  int constrain_Orbital_zeros;
  int constrain_Orbital_same;
  /**
     Just because the determinants shared orbitals when they
     were defined using the input file, that doesn't mean we
     can't optimize them independently.

     The restart file will have more orbitals than the input file.
  */
  int link_Determinant_parameters;

  /**
     1 if Gaussian orbitals are to be replaced by exponentials that 
     satisfy the electron nucleus cusp consitions near nuclei and 0 otherwise.
  */
  int replace_electron_nucleus_cusps;

  /**
     1 if we'll use the energy cutoff recommended by
     DePasquale, Rothstein, Vrbik JCP, 89, 3629, 1988
     0 otherwise
   */
  string energy_cutoff_type;

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
     1 if the distance between each pair of electrons for every time step is to
     be collected and written out in a histogram, and 0 otherwise.
  */
  int write_electron_densities;

  /**
     The maximum distance, in atomic units, of the histograms for the electron
     densities.
  */
  double max_pair_distance;

  /**
     1 if the x coordinates of each parallel spin pair of electrons are to be
     collected in a 2D histogram to make a correlation diagram, and 0 
     otherwise.
  */
  int writePllxCorrelationDiagram;

  /**
     The minimum value for the x coordinate of the 2D parallel spin correlation
     diagram.
  */
  double pllxCorrelationDiagramMin;

  /**
     The maximum value for the x coordinate of the 2D parallel spin correlation
     diagram.
  */
  double pllxCorrelationDiagramMax;

  /**
     1 if the y coordinates of each parallel spin pair of electrons are to be
     collected in a 2D histogram to make a correlation diagram, and 0 
     otherwise.
  */
  int writePllyCorrelationDiagram;

  /**
     The minimum value for the y coordinate of the 2D parallel spin correlation
     diagram.
  */
  double pllyCorrelationDiagramMin;

  /**
     The maximum value for the y coordinate of the 2D parallel spin correlation
     diagram.
  */
  double pllyCorrelationDiagramMax;

  /**
     1 if the z coordinates of each parallel spin pair of electrons are to be
     collected in a 2D histogram to make a correlation diagram, and 0 
     otherwise.
  */
  int writePllzCorrelationDiagram;

  /**
     The minimum value for the z coordinate of the 2D parallel spin correlation
     diagram.
  */
  double pllzCorrelationDiagramMin;

  /**
     The maximum value for the z coordinate of the 2D parallel spin correlation
     diagram.
  */
  double pllzCorrelationDiagramMax;

  /**
     1 if the x coordinates of each opposite spin pair of electrons are to be
     collected in a 2D histogram to make a correlation diagram, and 0 
     otherwise.
  */
  int writeOppxCorrelationDiagram;

  /**
     The minimum value for the x coordinate of the 2D opposite spin correlation
     diagram.
  */
  double oppxCorrelationDiagramMin;

  /**
     The maximum value for the x coordinate of the 2D opposite spin correlation
     diagram.
  */
  double oppxCorrelationDiagramMax;

  /**
     1 if the y coordinates of each opposite spin pair of electrons are to be
     collected in a 2D histogram to make a correlation diagram, and 0 
     otherwise.
  */
  int writeOppyCorrelationDiagram;

  /**
     The minimum value for the y coordinate of the 2D opposite spin correlation
     diagram.
  */
  double oppyCorrelationDiagramMin;

  /**
     The maximum value for the y coordinate of the 2D opposite spin correlation
     diagram
  */
  double oppyCorrelationDiagramMax;

  /**
     1 if the z coordinates of each opposite spin pair of electrons are to be
     collected in a 2D histogram to make a correlation diagram, and 0 
     otherwise.
  */
  int writeOppzCorrelationDiagram;

  /**
     The minimum value for the z coordinate of the 2D opposite spin correlation
     diagram.
  */
  double oppzCorrelationDiagramMin;

  /**
     The maximum value for the z coordinate of the 2D opposite spin correlation
     diagram.
  */
  double oppzCorrelationDiagramMax;

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
    Hartree-Fock calculations using QMC, aimed at producing basis independent 
    reference wavefunctions.  Computes V_eff for each electron by averaging 
    over multiple walkers and configurations.  
  */
  int use_hf_potential;

  /** 
    Number of electron samples to average v_eff over.
  */
  int hf_num_average;

  /**
    Option to keep trial energy fixed to its original value.
  */
  int lock_trial_energy;

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
     Check for some obvious flaws with parameter
     combinations. If something is irreparably wrong, then
     it will return false.

     Otherwise, it will try to correct any problems automatically
     and return true;
  */
  bool checkFlags();

  /**
     If we are going to allow the GPU to
     calculate a subset of walkers out of the
     walkers_per_pass batch, then we need to
     know how many it is going to process.

     If we never set this parameter in the
     input file, then set it to be the entire
     walkers_per_pass batch.

     If QMC_GPU is not set (we aren't using
     the GPU, then return ZERO, since the
     GPU will be calculating ZERO of the
     results).

     @return number of walkers the GPU will process per pass
  */
  int getNumGPUWalkers()
  {
#if defined QMC_GPU
    if(gpu_walkers_per_pass < 0)
    {
      return walkers_per_pass;
    } else {
      return gpu_walkers_per_pass;
    }
#else
    return 0;
#endif
  }
  /**
     Write this object's state out to a stream. The same format is 
     used as in the QMC input file.
  */  
  friend ostream& operator <<(ostream& strm, QMCFlags& flags);
};

#endif


















