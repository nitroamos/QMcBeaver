
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

#include "QMCFlags.h"
#include "QMCLineSearchStepLengthSelectionAlgorithm.h"
#include "QMCLineSearchStepLengthSelectionFactory.h"

QMCFlags::QMCFlags()
{
  my_rank = 0;
  nprocs  = 1;
}

void QMCFlags::read_flags(string InFileName)
{ 
  /*************************************************************************
    Set parameter defaults
  *************************************************************************/
  //QMC Parameters
  run_type                       = "variational";
  dt                             = 0.01;
  accel_delta                    = 5;
  accel_tm                       = 0.5*3.14159265359;

  desired_convergence            = 0.0;
  max_time_steps                 = 2000000;
  max_time                       = -1.0;
  number_of_walkers              = 100;
  one_e_per_iter                 = 0;
  use_surfer                     = 0;
  set_debug                      = 0;

  //Initialization parameters
  walker_initialization_method   = "dans_walker_initialization";
  walker_initialization_combinations = 3;
  dt_equilibration               = 0.01;
  equilibration_steps            = 2000;
  equilibration_function         = "ramp";
  CKAnnealingEquilibration1_parameter = 500;
  use_equilibration_array        = 0;

  //Weights, branching, fusion
  branching_method               = "nonunit_weight_branching";
  walker_reweighting_method      = "umrigar93_probability_weighted";
  branching_threshold            = 2.0;
  fusion_threshold               = 0.45;
  correct_population_size_bias   = 1;
  population_control_parameter   = 1.0;
  old_walker_acceptance_parameter = 50;
  warn_verbosity                 = 1;
  rel_cutoff                     = 100.0;
  limit_branching                = 1;

  branch_age_toolazy             = 4;

  //*
  branch_dWgrowth_toofast        = 5.0;
  branch_dR_badE                 = 100.0;
  branch_W_tooheavy              = 100.0;
  /*/
  branch_dWgrowth_toofast        = 5.0;
  branch_dR_badE                 = 100.0;
  branch_W_tooheavy              = 4.0;
  //*/

  //Green's function parameters
  sampling_method                = "umrigar93_importance_sampling";
  QF_modification_type           = "umrigar93_unequalelectrons";
  energy_modification_type       = "umrigar93";
  energy_cutoff_type             = "umrigar93";
  umrigar93_equalelectrons_parameter = 0.5;  // in (0,1]
  synchronize_dmc_ensemble       = 0;
  synchronize_dmc_ensemble_interval = 1000;

  //Wavefunction parameters
  trial_function_type            = "restricted";
  Ndeterminants                  = 1;

  //Other QMcBeaver improvements or added functionality
  replace_electron_nucleus_cusps = 1;
  print_replacement_orbitals     = 0;
  calculate_bf_density           = 0;
  use_hf_potential               = 0;
  hf_num_average                 = 100;
  lock_trial_energy              = 0;
  nuclear_derivatives            = "none";
  vector<double> fwInvHartrees;
  fwInvHartrees.clear();
  fwInvHartrees.push_back(0);

  //Computation parameters
  parallelization_method         = "manager_worker";
  iseed                          = 0;
  use_basis_function_interpolation = 0;
  number_basis_function_interpolation_grid_points = 1000;
  basis_function_interpolation_first_point        = 1e-10;
  pseudo_gridLevel               = 1;
  pseudo_cutoff                  = 1e-4;
  use_pseudopotential            = 0;
  walkers_per_pass               = 1;
  mpireduce_interval             = 100;
  mpipoll_interval               = 5;

  //Output parameters
  output_interval                = 1000;
  checkpoint_interval            = 100000;
  checkpoint                     = 0;
  use_available_checkpoints      = 0;
  checkin_file_name              = "";
  checkout_file_name             = "";
  zero_out_checkpoint_statistics = 0;
  checkpoint_energy_only         = 0;
  print_transient_properties     = 0;
  print_transient_properties_interval = 10000;
  print_configs                  = 0;
  print_config_frequency         = 50;
  temp_dir                       = "./";
  write_all_energies_out         = 0;
  write_electron_densities       = 0;
  max_pair_distance              = -1;

  writePllxCorrelationDiagram    = 0;
  writePllyCorrelationDiagram    = 0;
  writePllzCorrelationDiagram    = 0;
  writeOppxCorrelationDiagram    = 0;
  writeOppyCorrelationDiagram    = 0;
  writeOppzCorrelationDiagram    = 0;

  //Wavefunction optmization parameters
  optimize_Psi                   = 0;
  max_optimize_Psi_steps         = 30;
  optimize_Psi_barrier_parameter = 1.0;
  calculate_Derivatives          = 0;
  optimize_Psi_criteria          = "generalized_eigenvector";
  optimize_Psi_method            = "automatic";
  a_diag                         = -1e-5;
  ksi                            = 0.5;
  equilibrate_every_opt_step     = 1;
  equilibrate_first_opt_step     = 1;
  numerical_derivative_surface   = "umrigar88";
  line_search_step_length        = "None";
  optimization_max_iterations    = 100;
  optimization_error_tolerance   = 0.001;
  ck_genetic_algorithm_1_population_size = 1000;
  ck_genetic_algorithm_1_mutation_rate = 0.2;
  ck_genetic_algorithm_1_initial_distribution_deviation = 1.0;
  singularity_penalty_function_parameter = 1.0e-6;

  optimize_UD_Jastrows           = 1;
  optimize_UU_Jastrows           = 1;
  optimize_DD_Jastrows           = 1;
  optimize_EN_Jastrows           = 1;
  optimize_NEE_Jastrows          = 1;
  optimize_L                     = 1;
  optimize_CI                    = 1;
  optimize_Orbitals              = 0;
  link_Jastrow_parameters        = 1;
  link_NEE_Jastrows              = 2;
  link_Orbital_parameters        = 1;
  constrain_Orbital_zeros        = 1;
  constrain_Orbital_same         = 0;
  link_Determinant_parameters    = 1;

  use_three_body_jastrow         = 0;
  use_jastrow                    = 1;
  detailed_energies              = -1;

  reproduce_NE_with_NEE_jastrow  = 0;
  reproduce_EE_with_NEE_jastrow  = 0;

  //Difficult to categorize
  programmersLongs.clear();
  chip_and_mike_are_cool         = "false";

  /*************************************************************************
    Read in all the user specified parameters
  *************************************************************************/
  ifstream input_file(InFileName.c_str());
  string temp_string;

  if(!input_file)
    {
      cerr << "ERROR: Can't open input " << InFileName.c_str() << endl;
      exit(1);
    }

  input_file >> temp_string;
  while(temp_string != "&flags")
    {
      input_file >> temp_string;
    }

  //skipping the line that says "&flags"
  input_file >> temp_string;

  while((temp_string != "&") && (input_file.eof() != 1))
    {
      if(temp_string == "run_type")
        {
          input_file >> run_type;
        }
      else if(temp_string == "use_surfer")
        {
          input_file >> temp_string;
	  use_surfer = atoi(temp_string.c_str());
        }
      else if(temp_string == "set_debug")
        {
          input_file >> temp_string;
	  set_debug = atoi(temp_string.c_str());
        }
      else if(temp_string == "one_e_per_iter")
        {
          input_file >> temp_string;
	  one_e_per_iter = atoi(temp_string.c_str());
        }
      else if(temp_string == "temp_dir")
        {
          input_file >> temp_dir;
        }
      else if(temp_string == "parallelization_method")
        {
          input_file >> parallelization_method;
        }
      else if(temp_string == "walker_initialization_method")
        {
          input_file >> walker_initialization_method;
        }
      else if(temp_string == "walker_initialization_combinations")
        {
          input_file >> temp_string;
          walker_initialization_combinations = atoi(temp_string.c_str());
        }
      else if(temp_string == "iseed")
        {
	  // At this point, my_rank and nprocs have already been set.  We check
	  // the iseed to be sure that it does not overflow the int data type.

          input_file >> temp_string;
          iseed = atol(temp_string.c_str());
        }
      else if(temp_string == "sampling_method")
        {
          input_file >> sampling_method;
        }
      else if(temp_string == "nuclear_derivatives")
        {
          input_file >> nuclear_derivatives;
        }
      else if(temp_string == "QF_modification_type")
        {
          input_file >> QF_modification_type;
        }
      else if(temp_string == "energy_modification_type")
        {
          input_file >> energy_modification_type;
        }
      else if(temp_string == "energy_cutoff_type")
        {
          input_file >> energy_cutoff_type;
        }
      else if(temp_string == "umrigar93_equalelectrons_parameter")
        {
          input_file >> temp_string;
          umrigar93_equalelectrons_parameter = atof(temp_string.c_str());
        }
      else if(temp_string == "walker_reweighting_method")
        {
          input_file >> walker_reweighting_method;
        }
      else if(temp_string == "branching_method")
        {
          input_file >> branching_method;
        }
      else if(temp_string == "branching_threshold")
        {
          input_file >> temp_string;
          branching_threshold = atof(temp_string.c_str());
        }
      else if(temp_string == "fusion_threshold")
        {
          input_file >> temp_string;
          fusion_threshold = atof(temp_string.c_str());
        }
      else if(temp_string == "synchronize_dmc_ensemble")
        {
          input_file >> temp_string;
          synchronize_dmc_ensemble = atoi(temp_string.c_str());
        }
      else if(temp_string == "synchronize_dmc_ensemble_interval")
        {
          input_file >> temp_string;
          synchronize_dmc_ensemble_interval = atoi(temp_string.c_str());
        }
      else if(temp_string == "old_walker_acceptance_parameter")
        {
          input_file >> temp_string;
          old_walker_acceptance_parameter = atoi(temp_string.c_str());
        }
      else if(temp_string == "warn_verbosity")
        {
          input_file >> temp_string;
          warn_verbosity = atoi(temp_string.c_str());
        }
      else if(temp_string == "rel_cutoff")
        {
          input_file >> temp_string;
          rel_cutoff = atof(temp_string.c_str());
        }
      else if(temp_string == "limit_branching")
        {
          input_file >> temp_string;
          limit_branching = atoi(temp_string.c_str());
        }
      else if(temp_string == "use_basis_function_interpolation")
        {
          input_file >> temp_string;
          use_basis_function_interpolation = atoi(temp_string.c_str());
        }
      else if(temp_string == "number_basis_function_interpolation_grid_points")
        {
          input_file >> temp_string;
          number_basis_function_interpolation_grid_points =
            atoi(temp_string.c_str());
        }
      else if(temp_string == "basis_function_interpolation_first_point")
        {
          input_file >> temp_string;
          basis_function_interpolation_first_point = atof(temp_string.c_str());
        }
      else if(temp_string == "dt")
        {
          input_file >> temp_string;
          dt = atof(temp_string.c_str());
        }
      else if(temp_string == "accel_delta")
        {
          input_file >> temp_string;
          accel_delta = atof(temp_string.c_str());
        }
      else if(temp_string == "accel_tm")
        {
          input_file >> temp_string;
          accel_tm = atof(temp_string.c_str());
        }
      else if(temp_string == "desired_convergence")
        {
          input_file >> temp_string;
          desired_convergence = atof(temp_string.c_str());
        }
      else if(temp_string == "max_time_steps")
        {
          input_file >> temp_string;
          max_time_steps = atol(temp_string.c_str());
	  original_max_time_steps = max_time_steps;
        }
      else if(temp_string == "max_time")
        {
          input_file >> temp_string;
          max_time = atof(temp_string.c_str());
        }
      else if(temp_string == "future_walking")
        {
          input_file >> ws;//this removes any leading whitespace
          while(input_file.peek() >= '0' && input_file.peek() <= '9')
            {
              input_file >> temp_string;
              fwInvHartrees.push_back(atof(temp_string.c_str()));
              input_file >> ws;
            }
        }
      else if(temp_string == "walkers_per_pass")
        {
          input_file >> temp_string;
          walkers_per_pass = atol(temp_string.c_str());
        }
      else if(temp_string == "gpu_walkers_per_pass")
        {
          input_file >> temp_string;
          gpu_walkers_per_pass = atol(temp_string.c_str());
        }
      else if(temp_string == "number_of_walkers")
        {
          input_file >> temp_string;
          number_of_walkers = atol(temp_string.c_str());
        }
      else if(temp_string == "output_interval")
        {
          input_file >> temp_string;
          output_interval = atoi(temp_string.c_str());
        }
      else if(temp_string == "checkpoint_interval")
        {
          input_file >> temp_string;
          checkpoint_interval = atoi(temp_string.c_str());
        }
      else if(temp_string == "checkpoint")
        {
          input_file >> temp_string;
          checkpoint = atoi(temp_string.c_str());
        }
      else if(temp_string == "use_available_checkpoints")
        {
          input_file >> temp_string;
          use_available_checkpoints = atoi(temp_string.c_str());
        }
      else if(temp_string == "checkpoint_input_name")
        {
          input_file >> checkin_file_name;
        }
      else if(temp_string == "checkpoint_output_name")
        {
          input_file >> checkout_file_name;
        }
      else if(temp_string == "checkpoint_energy_only")
	{
	  input_file >> temp_string;
	  checkpoint_energy_only = atoi(temp_string.c_str());
	}
      else if(temp_string == "equilibration_steps")
        {
          input_file >> temp_string;
          equilibration_steps = atoi(temp_string.c_str());
        }
      else if(temp_string == "dt_equilibration")
        {
          input_file >> temp_string;
          dt_equilibration = atof(temp_string.c_str());
        }
      else if(temp_string == "equilibration_function")
        {
          input_file >> equilibration_function;
        }
      else if(temp_string == "CKAnnealingEquilibration1_parameter")
        {
          input_file >> temp_string;
          CKAnnealingEquilibration1_parameter = atoi(temp_string.c_str());
        }
      else if (temp_string == "use_equilibration_array")
        {
          input_file >> temp_string;
          use_equilibration_array = atoi(temp_string.c_str());
        }
      else if(temp_string == "mpireduce_interval")
        {
          input_file >> temp_string;
          mpireduce_interval = atoi(temp_string.c_str());
        }
      else if(temp_string == "mpipoll_interval")
        {
          input_file >> temp_string;
          mpipoll_interval = atoi(temp_string.c_str());
        }
      else if(temp_string == "atoms")
        {
          input_file >> temp_string;
          Natoms = atoi(temp_string.c_str());
        }
      else if(temp_string == "charge")
        {
          input_file >> temp_string;
          charge = atoi(temp_string.c_str());
        }
      else if(temp_string == "norbitals")
        {
          input_file >> temp_string;
          Norbitals = atoi(temp_string.c_str());
        }
      else if(temp_string == "nbasisfunc")
        {
          input_file >> temp_string;
          Nbasisfunc = atoi(temp_string.c_str());
        }
      else if(temp_string == "ndeterminants")
        {
          input_file >> temp_string;
          Ndeterminants = atoi(temp_string.c_str());
        }
      else if(temp_string == "trial_function_type")
        {
          input_file >> trial_function_type;
        }
      else if(temp_string == "pseudo_gridLevel")
        {
          input_file >> temp_string;
	  pseudo_gridLevel = atoi(temp_string.c_str());
        }
      else if(temp_string == "pseudo_cutoff")
        {
          input_file >> temp_string;
	  pseudo_cutoff = atof(temp_string.c_str());
        }
      else if(temp_string == "calculate_bf_density")
        {
          input_file >> temp_string;
          calculate_bf_density = atoi(temp_string.c_str());
        }
      else if(temp_string == "energy")
        {
          input_file >> temp_string;
          energy_trial = atof(temp_string.c_str());
        }
      else if(temp_string == "correct_population_size_bias")
        {
          input_file >> temp_string;
          correct_population_size_bias = atoi(temp_string.c_str());
        }
      else if(temp_string == "print_transient_properties")
        {
          input_file >> temp_string;
          print_transient_properties = atoi(temp_string.c_str());
        }
      else if(temp_string == "print_transient_properties_interval")
        {
          input_file >> temp_string;
          print_transient_properties_interval = atoi(temp_string.c_str());
        }
      else if(temp_string == "print_configs")
        {
          input_file >> temp_string;
          print_configs = atoi(temp_string.c_str());
        }
      else if(temp_string == "print_config_frequency")
        {
          input_file >> temp_string;
          print_config_frequency = atoi(temp_string.c_str());
        }
      else if(temp_string == "optimize_Psi")
        {
          input_file >> temp_string;
          optimize_Psi = atoi(temp_string.c_str());
        }
      else if(temp_string == "optimize_EE_Jastrows")
        {
          input_file >> temp_string;
          optimize_UD_Jastrows = atoi(temp_string.c_str());
          optimize_UU_Jastrows = atoi(temp_string.c_str());
          optimize_DD_Jastrows = atoi(temp_string.c_str());
        }
      else if(temp_string == "optimize_UD_Jastrows")
        {
          input_file >> temp_string;
          optimize_UD_Jastrows = atoi(temp_string.c_str());
        }
      else if(temp_string == "optimize_UU_Jastrows")
        {
          input_file >> temp_string;
          optimize_UU_Jastrows = atoi(temp_string.c_str());
        }
      else if(temp_string == "optimize_DD_Jastrows")
        {
          input_file >> temp_string;
          optimize_DD_Jastrows = atoi(temp_string.c_str());
        }
      else if(temp_string == "optimize_EN_Jastrows")
        {
          input_file >> temp_string;
          optimize_EN_Jastrows = atoi(temp_string.c_str());
        }
      else if(temp_string == "optimize_NEE_Jastrows")
        {
          input_file >> temp_string;
          optimize_NEE_Jastrows = atoi(temp_string.c_str());
        }
      else if(temp_string == "use_three_body_jastrow")
        {
          input_file >> temp_string;
          use_three_body_jastrow = atoi(temp_string.c_str());
        }
      else if(temp_string == "reproduce_NE_with_NEE_jastrow")
	{
	  input_file >> temp_string;
	  reproduce_NE_with_NEE_jastrow = atoi(temp_string.c_str());
	}
      else if(temp_string == "reproduce_EE_with_NEE_jastrow")
	{
	  input_file >> temp_string;
	  reproduce_EE_with_NEE_jastrow = atoi(temp_string.c_str());
	}
      else if(temp_string == "optimize_L")
	{
	  input_file >> temp_string;
	  optimize_L = atoi(temp_string.c_str());
	}
      else if(temp_string == "optimize_CI")
        {
          input_file >> temp_string;
          optimize_CI = atoi(temp_string.c_str());
        }
      else if(temp_string == "optimize_Orbitals")
        {
          input_file >> temp_string;
          optimize_Orbitals = atoi(temp_string.c_str());
        }
      else if(temp_string == "optimize_Psi_barrier_parameter")
        {
          input_file >> temp_string;
          optimize_Psi_barrier_parameter = atof(temp_string.c_str());
        }
      else if(temp_string == "max_optimize_Psi_steps")
        {
          input_file >> temp_string;
          max_optimize_Psi_steps = atoi(temp_string.c_str());
        }
      else if(temp_string == "optimize_Psi_criteria")
        {
          input_file >> optimize_Psi_criteria;
        }
      else if(temp_string == "optimize_Psi_method")
        {
          input_file >> optimize_Psi_method;
        }
      else if(temp_string == "a_diag")
        {
          input_file >> temp_string;
	  a_diag = atof(temp_string.c_str());
        }
      else if(temp_string == "ksi")
        {
          input_file >> temp_string;
	  ksi = atof(temp_string.c_str());
        }
      else if(temp_string == "singularity_penalty_function_parameter")
        {
          input_file >> temp_string;
          singularity_penalty_function_parameter = atof(temp_string.c_str());
        }
      else if(temp_string == "line_search_step_length")
        {
          input_file >> line_search_step_length;
        }
      else if(temp_string == "optimization_max_iterations")
        {
          input_file >> temp_string;
          optimization_max_iterations = atoi(temp_string.c_str());
        }
      else if(temp_string == "optimization_error_tolerance")
        {
          input_file >> temp_string;
          optimization_error_tolerance = atof(temp_string.c_str());
        }
      else if(temp_string == "ck_genetic_algorithm_1_population_size")
        {
          input_file >> temp_string;
          ck_genetic_algorithm_1_population_size = atoi(temp_string.c_str());
        }
      else if(temp_string == "ck_genetic_algorithm_1_mutation_rate")
        {
          input_file >> temp_string;
          ck_genetic_algorithm_1_mutation_rate = atof(temp_string.c_str());
        }
      else if(temp_string ==
              "ck_genetic_algorithm_1_initial_distribution_deviation")
        {
          input_file >> temp_string;
          ck_genetic_algorithm_1_initial_distribution_deviation =
            atof(temp_string.c_str());
        }
      else if(temp_string == "numerical_derivative_surface")
        {
          input_file >> numerical_derivative_surface;
        }
      else if(temp_string == "link_Jastrow_parameters")
        {
          input_file >> temp_string;
          link_Jastrow_parameters = atoi(temp_string.c_str());
        }
      else if(temp_string == "link_NEE_Jastrows")
        {
          input_file >> temp_string;
          link_NEE_Jastrows = atoi(temp_string.c_str());
        }
      else if(temp_string == "link_Orbital_parameters")
        {
          input_file >> temp_string;
          link_Orbital_parameters = atoi(temp_string.c_str());
        }
      else if(temp_string == "constrain_Orbital_zeros")
        {
          input_file >> temp_string;
          constrain_Orbital_zeros = atoi(temp_string.c_str());
        }
      else if(temp_string == "constrain_Orbital_same")
        {
          input_file >> temp_string;
          constrain_Orbital_same = atoi(temp_string.c_str());
        }
      else if(temp_string == "link_Determinant_parameters")
        {
          input_file >> temp_string;
          link_Determinant_parameters = atoi(temp_string.c_str());
        }
      else if(temp_string == "replace_electron_nucleus_cusps")
	{
	  input_file >> temp_string;
	  replace_electron_nucleus_cusps = atoi(temp_string.c_str());
	}
      else if(temp_string == "print_replacement_orbitals")
	{
	  input_file >> temp_string;
	  print_replacement_orbitals = atoi(temp_string.c_str());
	}
      else if(temp_string == "equilibrate_every_opt_step")
        {
          input_file >> temp_string;
          equilibrate_every_opt_step = atoi(temp_string.c_str());
        }
      else if(temp_string == "equilibrate_first_opt_step")
        {
          input_file >> temp_string;
          equilibrate_first_opt_step = atoi(temp_string.c_str());
        }
      else if(temp_string == "population_control_parameter")
        {
          input_file >> temp_string;
          population_control_parameter = atof(temp_string.c_str());
        }
      else if(temp_string == "write_electron_densities")
        {
          input_file >> temp_string;
          write_electron_densities = atoi(temp_string.c_str());
        }
      else if(temp_string == "max_pair_distance")
	{
	  input_file >> temp_string;
	  max_pair_distance = atof(temp_string.c_str());
	}
      else if(temp_string == "writePllxCorrelationDiagram")
	{
	  input_file >> temp_string;
	  writePllxCorrelationDiagram = atoi(temp_string.c_str());
	}
      else if(temp_string == "pllxCorrelationDiagramMin")
	{
	  input_file >> temp_string;
	  pllxCorrelationDiagramMin = atof(temp_string.c_str());
	}
      else if(temp_string == "pllxCorrelationDiagramMax")
	{
	  input_file >> temp_string;
	  pllxCorrelationDiagramMax = atof(temp_string.c_str());
	}
      else if(temp_string == "writePllyCorrelationDiagram")
	{
	  input_file >> temp_string;
	  writePllyCorrelationDiagram = atoi(temp_string.c_str());
	}
      else if(temp_string == "pllyCorrelationDiagramMin")
	{
	  input_file >> temp_string;
	  pllyCorrelationDiagramMin = atof(temp_string.c_str());
	}
      else if(temp_string == "pllyCorrelationDiagramMax")
	{
	  input_file >> temp_string;
	  pllyCorrelationDiagramMax = atof(temp_string.c_str());
	}
      else if(temp_string == "writePllzCorrelationDiagram")
	{
	  input_file >> temp_string;
	  writePllzCorrelationDiagram = atoi(temp_string.c_str());
	}
      else if(temp_string == "pllzCorrelationDiagramMin")
	{
	  input_file >> temp_string;
	  pllzCorrelationDiagramMin = atof(temp_string.c_str());
	}
      else if(temp_string == "pllzCorrelationDiagramMax")
	{
	  input_file >> temp_string;
	  pllzCorrelationDiagramMax = atof(temp_string.c_str());
	}
      else if(temp_string == "writeOppxCorrelationDiagram")
	{
	  input_file >> temp_string;
	  writeOppxCorrelationDiagram = atoi(temp_string.c_str());
	}
      else if(temp_string == "oppxCorrelationDiagramMin")
	{
	  input_file >> temp_string;
	  oppxCorrelationDiagramMin = atof(temp_string.c_str());
	}
      else if(temp_string == "oppxCorrelationDiagramMax")
	{
	  input_file >> temp_string;
	  oppxCorrelationDiagramMax = atof(temp_string.c_str());
	}
      else if(temp_string == "writeOppyCorrelationDiagram")
	{
	  input_file >> temp_string;
	  writeOppyCorrelationDiagram = atoi(temp_string.c_str());
	}
      else if(temp_string == "oppyCorrelationDiagramMin")
	{
	  input_file >> temp_string;
	  oppyCorrelationDiagramMin = atof(temp_string.c_str());
	}
      else if(temp_string == "oppyCorrelationDiagramMax")
	{
	  input_file >> temp_string;
	  oppyCorrelationDiagramMax = atof(temp_string.c_str());
	}
      else if(temp_string == "writeOppzCorrelationDiagram")
	{
	  input_file >> temp_string;
	  writeOppzCorrelationDiagram = atoi(temp_string.c_str());
	}
      else if(temp_string == "oppzCorrelationDiagramMin")
	{
	  input_file >> temp_string;
	  oppzCorrelationDiagramMin = atof(temp_string.c_str());
	}
      else if(temp_string == "oppzCorrelationDiagramMax")
	{
	  input_file >> temp_string;
	  oppzCorrelationDiagramMax = atof(temp_string.c_str());
	}
      else if(temp_string == "write_all_energies_out")
        {
          input_file >> temp_string;
          write_all_energies_out = atoi(temp_string.c_str());
        }
      else if(temp_string == "zero_out_checkpoint_statistics" )
        {
          input_file >> temp_string;
          zero_out_checkpoint_statistics = atoi(temp_string.c_str());
        }
      else if(temp_string == "programmers_longs" )
        {
          input_file >> ws;//this removes any leading whitespace
          while(input_file.peek() >= '0' && input_file.peek() <= '9')
            {
              input_file >> temp_string;
              programmersLongs.push_back(atol(temp_string.c_str()));
              input_file >> ws;
            }
        }
      else if(temp_string == "chip_and_mike_are_cool")
        {
          input_file >> chip_and_mike_are_cool;
        }
      else if(temp_string == "use_hf_potential")
        {
          input_file >> temp_string;
          use_hf_potential = atoi(temp_string.c_str());
        }
      else if(temp_string == "hf_num_average")
        {
          input_file >> temp_string;
          hf_num_average = atoi(temp_string.c_str());
        }
      else if(temp_string == "lock_trial_energy")
        {
          input_file >> temp_string;
          lock_trial_energy = atoi(temp_string.c_str());
        }
      else if(temp_string.find("#",0) != string::npos)
        {
	  //it's a comment if one of the first 3 characters is a '#'
	  //IMPORTANT NOTE: this only checks for comments where parameter
	  //names are expected, not where parameter values are read in!
          while(input_file.peek() != '\n')
            {
	      //just peal one character off at a time so that we
	      //don't have to worry about whitespace being skipped
	      input_file.get();
            }
        }
      else
        {
          clog << "Warning: Unknown input flag: " << temp_string << endl;
	  exit(0);
        }
      input_file >> temp_string;
    }
  input_file.close();

  if(!checkFlags())
    {
      clog << "Error: serious flaws with parameter choices; no fix available." << endl;
      exit(1);
    }

  /*************************************************************************
    Finally, fill in the parameters to get us ready for the calculation.
  *************************************************************************/
  /*
    This should make all the elements in the future walking vector unique
    and sorted.
  */
  sort(fwInvHartrees.begin(), fwInvHartrees.end());
  vector<double>::iterator last = unique(fwInvHartrees.begin(), fwInvHartrees.end());
  fwInvHartrees.erase(last,fwInvHartrees.end());
  
  //since it's better to input future walking lengths in terms of
  // Hatrees^{-1} rather than block length, we need to convert it
  for(unsigned int i=0; i<fwInvHartrees.size(); i++)
    {
      future_walking.push_back( (int)(fwInvHartrees[i]/dt) );
    }

  set_filenames(InFileName);

  dt_effective = dt;
  dt_run = dt;
  dt = dt_equilibration;
  number_of_walkers_initial = number_of_walkers;
  energy_estimated          = energy_trial;
  energy_estimated_original = energy_estimated;

  if(Ndeterminants == 1)
    optimize_CI = 0;
}

void QMCFlags::set_filenames(string runfile)
{
  string tempstring;
  string basename;
  int restart_number = 0;

  input_file_name   = runfile;
  string file_name  = input_file_name.substr(0,input_file_name.size()-5);
  base_file_name    = file_name;
  output_file_name  = file_name + ".qmc";
  results_file_name = file_name + ".rslts";
  density_file_name = file_name + ".density";
  force_file_name   = file_name + ".force";

  if(checkin_file_name == "")
    {
      checkin_file_name = base_file_name;
    }
  if(checkout_file_name == "")
    {
      checkout_file_name = base_file_name;
    }

  char my_rank_string[32];
#ifdef _WIN32
  _snprintf( my_rank_string, 32, "%d", my_rank );
#else
  snprintf( my_rank_string, 32, "%d", my_rank );
#endif
  config_file_name = temp_dir + base_file_name + "." +my_rank_string + ".cfgs";

  basename = file_name;
  if(file_name[file_name.size()-3] == '.')
    {
      tempstring = file_name.substr(file_name.size()-2,file_name.size()-1);
      restart_number = atoi(tempstring.c_str());
      basename = input_file_name.substr(0,file_name.size()-3);
    }

  restart_number += 1;

  char Tens, Ones;
  int tens = restart_number/10;
  int ones = restart_number%10;

  switch(tens)
    {
        case 0:     Tens = '0'; break;
        case 1:     Tens = '1'; break;
        case 2:     Tens = '2'; break;
        case 3:     Tens = '3'; break;
        case 4:     Tens = '4'; break;
        case 5:     Tens = '5'; break;
        case 6:     Tens = '6'; break;
        case 7:     Tens = '7'; break;
        case 8:     Tens = '8'; break;
        case 9:     Tens = '9'; break;
        default:
        {
          cerr << "Error in naming restart!" << endl;
          exit(1);
        }
    }
  switch(ones)
    {
        case 0:     Ones = '0'; break;
        case 1:     Ones = '1'; break;
        case 2:     Ones = '2'; break;
        case 3:     Ones = '3'; break;
        case 4:     Ones = '4'; break;
        case 5:     Ones = '5'; break;
        case 6:     Ones = '6'; break;
        case 7:     Ones = '7'; break;
        case 8:     Ones = '8'; break;
        case 9:     Ones = '9'; break;
        default:
        {
          cerr << "Error in naming restart!" << endl;
          exit(1);
        }
    }

  restart_file_name = basename + '.' + Tens + Ones + ".ckmf";
}

ostream& operator <<(ostream& strm, QMCFlags& flags)
{
  strm.setf( ios::fixed );
  strm << "&flags" << endl;

  /*
    I tried to categorize all the parameters so that similar parameters are
    close to each other when output. This categorization is somewhat arbitrary...

    After categorization, I attempted to order the categories according to how
    likely someone is to want to change something in that category. Obviously,
    there is no perfect way of doing this.
   */
  
  strm << "# Parameters for QMC\n";
  strm << "run_type\n " << flags.run_type << endl;
  strm << "dt\n " << flags.dt_run << endl;
  strm << "dt_equilibration\n " << flags.dt_equilibration << endl;
  strm << "number_of_walkers\n " << flags.number_of_walkers_initial<< endl;
  //strm << "max_time_steps\n " << flags.original_max_time_steps << endl;
  strm << "max_time_steps\n " << flags.max_time_steps << endl;
  strm << "equilibration_steps\n " << flags.equilibration_steps << endl;
  strm << "desired_convergence\n " << flags.desired_convergence << endl;
  strm << "iseed\n " << flags.iseed << endl;
  strm << "optimize_Psi\n " << flags.optimize_Psi << endl;
  strm << "max_time\n " << flags.max_time << endl;
  strm << "one_e_per_iter\n " << flags.one_e_per_iter << endl;
  strm << "output_interval\n " << flags.output_interval << endl;

  strm << "\n# Parameters for wavefunction optimization\n";
  strm << "optimize_UD_Jastrows\n " << flags.optimize_UD_Jastrows << endl;
  strm << "optimize_UU_Jastrows\n " << flags.optimize_UU_Jastrows << endl;
  strm << "optimize_DD_Jastrows\n " << flags.optimize_DD_Jastrows << endl;
  strm << "optimize_EN_Jastrows\n " << flags.optimize_EN_Jastrows << endl;
  strm << "optimize_NEE_Jastrows\n " << flags.optimize_NEE_Jastrows << endl;
  strm << "optimize_L\n " << flags.optimize_L << endl;
  strm << "optimize_CI\n " << flags.optimize_CI << endl;
  strm << "optimize_Orbitals\n " << flags.optimize_Orbitals << endl;
  strm << "optimize_Psi_method\n " << flags.optimize_Psi_method << endl;
  strm << "optimize_Psi_criteria\n " << flags.optimize_Psi_criteria << endl;
  strm << "a_diag\n " << flags.a_diag << endl;
  strm << "ksi\n " << flags.ksi << endl;
  strm << "max_optimize_Psi_steps\n " << flags.max_optimize_Psi_steps << endl;
  strm << "equilibrate_first_opt_step\n "
       << flags.equilibrate_first_opt_step << endl;
  strm << "equilibrate_every_opt_step\n "
       << flags.equilibrate_every_opt_step << endl;
  strm << "optimization_max_iterations\n "
       << flags.optimization_max_iterations << endl;
  strm << "optimization_error_tolerance\n "
       << flags.optimization_error_tolerance << endl;
  strm << "singularity_penalty_function_parameter\n "
       << flags.singularity_penalty_function_parameter << endl;
  strm << "optimize_Psi_barrier_parameter\n "
       << flags.optimize_Psi_barrier_parameter << endl;
  strm << "numerical_derivative_surface\n "
       << flags.numerical_derivative_surface << endl;
  strm << "line_search_step_length\n "
       << flags.line_search_step_length << endl;
  strm << "ck_genetic_algorithm_1_population_size\n "
       << flags.ck_genetic_algorithm_1_population_size << endl;
  strm << "ck_genetic_algorithm_1_mutation_rate\n "
       << flags.ck_genetic_algorithm_1_mutation_rate << endl;
  strm << "ck_genetic_algorithm_1_initial_distribution_deviation\n "
       << flags.ck_genetic_algorithm_1_initial_distribution_deviation << endl;

  strm << "\n# Parameters specific to the Green's function\n";
  strm << "sampling_method\n " << flags.sampling_method << endl;
  strm << "QF_modification_type\n " << flags.QF_modification_type << endl;
  strm << "umrigar93_equalelectrons_parameter\n "
       << flags.umrigar93_equalelectrons_parameter << endl;
  strm << "warn_verbosity\n " << flags.warn_verbosity << endl;
  strm << "rel_cutoff\n " << flags.rel_cutoff << endl;
  strm << "limit_branching\n " << flags.limit_branching << endl;
  strm << "energy_modification_type\n "
       << flags.energy_modification_type << endl;
  strm << "energy_cutoff_type\n "
       << flags.energy_cutoff_type << endl;
  strm << "lock_trial_energy\n " << flags.lock_trial_energy << endl;
  strm << "synchronize_dmc_ensemble\n "
       << flags.synchronize_dmc_ensemble << endl;
  strm << "synchronize_dmc_ensemble_interval\n "
       << flags.synchronize_dmc_ensemble_interval << endl;
  
  strm << "\n# Parameters specific to weights, branching, and fusion\n";
  strm << "walker_reweighting_method\n "
       << flags.walker_reweighting_method << endl;
  strm << "branching_method\n " << flags.branching_method << endl;
  strm << "branching_threshold\n " << flags.branching_threshold << endl;
  strm << "fusion_threshold\n " << flags.fusion_threshold << endl;
  strm << "population_control_parameter\n "
       << flags.population_control_parameter << endl;
  strm << "correct_population_size_bias\n "
       << flags.correct_population_size_bias << endl;
  strm << "old_walker_acceptance_parameter\n "
       << flags.old_walker_acceptance_parameter << endl;

  strm << "\n# Parameters for initialization\n";
  strm << "use_equilibration_array\n " << flags.use_equilibration_array << endl;
  strm << "equilibration_function\n " << flags.equilibration_function << endl;
  strm << "CKAnnealingEquilibration1_parameter\n "
       << flags.CKAnnealingEquilibration1_parameter << endl;
  strm << "walker_initialization_method\n "
       << flags.walker_initialization_method << endl;
  strm << "walker_initialization_combinations\n "
       << flags.walker_initialization_combinations << endl;
    
  strm << "\n# Parameters for added functionality/improvements\n";
  strm << "calculate_bf_density\n " << flags.calculate_bf_density << endl;
  strm << "use_hf_potential\n " << flags.use_hf_potential << endl;
  strm << "hf_num_average\n " << flags.hf_num_average << endl;
  strm << "replace_electron_nucleus_cusps\n " 
       << flags.replace_electron_nucleus_cusps << endl;
  strm << "print_replacement_orbitals\n " 
       << flags.print_replacement_orbitals << endl;
  strm << "nuclear_derivatives\n " << flags.nuclear_derivatives << endl;
  if(flags.future_walking.size() > 0)
    {
      strm << "future_walking\n";
      for(unsigned int i=0; i<flags.future_walking.size(); i++)
	{
	  strm << " " << (double)(flags.future_walking[i]*flags.dt);
	}
      strm << endl;
    }

  strm << "\n# Parameters relating to output\n";
  strm << "checkpoint\n " << flags.checkpoint << endl;
  strm << "checkpoint_interval\n " << flags.checkpoint_interval << endl;

  if(flags.checkpoint)
      strm << "use_available_checkpoints\n "
	   << 1 << endl;
  else
      strm << "use_available_checkpoints\n "
	   << flags.use_available_checkpoints << endl;

  strm << "checkpoint_input_name\n "
       << flags.checkout_file_name << endl;
  strm << "zero_out_checkpoint_statistics\n "
       << flags.zero_out_checkpoint_statistics << endl;
  strm << "checkpoint_energy_only\n " << flags.checkpoint_energy_only << endl;
  strm << "print_configs\n " << flags.print_configs << endl;
  strm << "print_config_frequency\n " << flags.print_config_frequency << endl;
  strm << "temp_dir\n " << flags.temp_dir << endl;
  strm << "write_all_energies_out\n " << flags.write_all_energies_out << endl;
  strm << "write_electron_densities\n "
       << flags.write_electron_densities << endl;
  strm << "max_pair_distance\n " << flags.max_pair_distance << endl;
  strm << "print_transient_properties\n "
       << flags.print_transient_properties << endl;
  strm << "print_transient_properties_interval\n "
       << flags.print_transient_properties_interval << endl;

  if(false){
    strm << "writePllxCorrelationDiagram\n " << flags.writePllxCorrelationDiagram
	 << endl;
    strm << "pllxCorrelationDiagramMin\n " << flags.pllxCorrelationDiagramMin
	 << endl;
    strm << "pllxCorrelationDiagramMax\n " << flags.pllxCorrelationDiagramMax
	 << endl;
    strm << "writePllyCorrelationDiagram\n " << flags.writePllyCorrelationDiagram
	 << endl;
    strm << "pllyCorrelationDiagramMin\n " << flags.pllyCorrelationDiagramMin
	 << endl;
    strm << "pllyCorrelationDiagramMax\n " << flags.pllyCorrelationDiagramMax
	 << endl;
    strm << "writePllzCorrelationDiagram\n " << flags.writePllzCorrelationDiagram
	 << endl;
    strm << "pllzCorrelationDiagramMin\n " << flags.pllzCorrelationDiagramMin
	 << endl;
    strm << "pllzCorrelationDiagramMax\n " << flags.pllzCorrelationDiagramMax
	 << endl;
    
    strm << "writeOppxCorrelationDiagram\n " << flags.writeOppxCorrelationDiagram
	 << endl;
    strm << "oppxCorrelationDiagramMin\n " << flags.oppxCorrelationDiagramMin
	 << endl;
    strm << "oppxCorrelationDiagramMax\n " << flags.oppxCorrelationDiagramMax
	 << endl;
    strm << "writeOppyCorrelationDiagram\n " << flags.writeOppyCorrelationDiagram
	 << endl;
    strm << "oppyCorrelationDiagramMin\n " << flags.oppyCorrelationDiagramMin
	 << endl;
    strm << "oppyCorrelationDiagramMax\n " << flags.oppyCorrelationDiagramMax
	 << endl;
    strm << "writeOppzCorrelationDiagram\n " << flags.writeOppzCorrelationDiagram
	 << endl;
    strm << "oppzCorrelationDiagramMin\n " << flags.oppzCorrelationDiagramMin
	 << endl;
    strm << "oppzCorrelationDiagramMax\n " << flags.oppzCorrelationDiagramMax
	 << endl;
  }
  
  strm << "\n# Parameters for computation/MPI\n";
  strm << "parallelization_method\n " << flags.parallelization_method << endl;
  strm << "mpireduce_interval\n " << flags.mpireduce_interval << endl;
  strm << "mpipoll_interval\n " << flags.mpipoll_interval << endl;
  strm << "walkers_per_pass\n " << flags.walkers_per_pass << endl;
  //strm << "gpu_walkers_per_pass\n " << flags.gpu_walkers_per_pass << endl;
  strm << "use_basis_function_interpolation\n "
       << flags.use_basis_function_interpolation << endl;
  strm << "number_basis_function_interpolation_grid_points\n "
       << flags.number_basis_function_interpolation_grid_points << endl;
  strm << "basis_function_interpolation_first_point\n "
       << flags.basis_function_interpolation_first_point << endl;

  strm << "\n# Parameters for the molecule and wavefunction\n";
  strm << "atoms\n " << flags.Natoms << endl;
  strm << "charge\n " << flags.charge << endl;

  if(flags.checkpoint)
    //This is if you want the checkpointed file to be able to duplicate exactly
    //a run that wasn't checkpointed.
    strm << "energy\n " << flags.energy_estimated_original << endl;
  else
    strm << "energy\n " << flags.energy_trial << endl;

  strm << "trial_function_type\n " << flags.trial_function_type << endl;
  strm << "pseudo_gridLevel\n " << flags.pseudo_gridLevel << endl;
  strm << "pseudo_cutoff\n " << flags.pseudo_cutoff << endl;
  strm << "norbitals\n " << flags.Norbitals << endl;
  strm << "nbasisfunc\n " << flags.Nbasisfunc << endl;
  strm << "ndeterminants\n " << flags.Ndeterminants << endl;
  strm << "link_Jastrow_parameters\n "
       << flags.link_Jastrow_parameters << endl;
  strm << "link_NEE_Jastrows\n "
       << flags.link_NEE_Jastrows << endl;
  strm << "link_Orbital_parameters\n "
       << flags.link_Orbital_parameters << endl;
  strm << "link_Determinant_parameters\n "
       << flags.link_Determinant_parameters << endl;
  strm << "reproduce_NE_with_NEE_jastrow\n " 
       << flags.reproduce_NE_with_NEE_jastrow << endl;
  strm << "reproduce_EE_with_NEE_jastrow\n "
       << flags.reproduce_EE_with_NEE_jastrow << endl;
  
  strm << "\n# Other parameters\n";
  strm << "chip_and_mike_are_cool\n " << flags.chip_and_mike_are_cool << endl;
  
  strm << "& " << endl;
  strm.unsetf( ios::fixed );
  return strm;
}

bool QMCFlags::checkFlags()
{
  /*************************************************************************
    Our warning section.
    Make sure that all the input parameters are compatible with each other
  *************************************************************************/
  if(run_type != "variational" && run_type != "diffusion")
    {
      clog << "ERROR: Unknown run_type: " << run_type << endl;
      return false;
    }
  
  if(max_time_steps <= equilibration_steps)
    {
      clog << "ERROR: (max_time_steps = " << max_time_steps
	   << ") < (equilibration_steps = " << equilibration_steps << ")" << endl; 
      return false;
    }

  if(limit_branching == 1 && run_type == "diffusion")
    {
      clog << "Laziness filter: " << branch_age_toolazy << endl;
      clog << "Fast growth:     " << branch_dWgrowth_toofast << endl;
      clog << "Bad Energy:      " << branch_dR_badE << endl;
      clog << "Too heavy:       " << branch_W_tooheavy << endl << endl;
    }

  if(dt_equilibration < 0)
    {
      dt_equilibration = dt;
    }

  if(sampling_method == "umrigar93_accelerated_sampling")
    {
      if(run_type != "variational"){
	clog << "Error: sampling method " << sampling_method <<
	  " can only be used with run_type = variational.";
	return false;
      }
      clog << "Using sampling_method = " << sampling_method << endl;
      clog << " accel_delta = " << accel_delta << endl;
      clog << " accel_tm    = " << accel_tm << endl;
    }

  if(dt > dt_equilibration)
    {
      clog << "Warning: dt (" << dt << ") > dt_equilibration (" << dt_equilibration << "), setting dt_equilibration = dt." << endl;
      dt_equilibration = dt;
    }

  if(rel_cutoff < 10.0)
    {
      clog << "Warning: the filter rel_cutoff = " << rel_cutoff << " might be too low..." << endl;
    }

  if(checkpoint == 1)
    {
      string filename =
	globalInput.flags.temp_dir
	+ "/"
	+ globalInput.flags.checkout_file_name;

      clog << "Notice: checkpoint files will be written into " << filename << endl;
    }

  if(one_e_per_iter == 1 && optimize_Psi == 1)
    {
      /*
	Currently:
	1) correlated sampling doesn't work with one_e_per_iter
	2) I haven't added updates for the parameter derivatives

	Maybe someone could use a different optimization scheme, but at
	the moment, it doesn't seem worth it.
      */
      clog << "Error: the code is unable to optimize psi using one_e_per_iter.\n";
      return false;
    }

  /**
     In allowing a harmonic oscillator model, we use some of the regular
     molecular parameters to help define our HO system.
  */
  if(trial_function_type == "harmonicoscillator")
    {
      if(walker_initialization_method != "amos_boring_initialization")
	{
	  walker_initialization_method = "amos_boring_initialization";
	  clog << "Warning: Setting walker_initialization_method to \"" << 
	    walker_initialization_method << "\" since trial_function_type = " <<
	    trial_function_type << endl;
	}

      if(Natoms != 0)
	{
	  clog << "Warning: Setting Natoms to 0 since trial_function_type = "
	       << trial_function_type << endl;
	  Natoms = 0;
	}

      if(Nbasisfunc <= 0)
	{
	  clog << "Error: set nbasisfunc to the number of dimensions for trial_function_type = "
	       << trial_function_type << endl;
	  return false;
	}

      if(charge == 0)
	{
	  clog << "Error: set charge to the number of electrons for trial_function_type = "
	       << trial_function_type << endl;
	  return false;
	}
    }

  bool needPrintedConfigs = false;
  if(optimize_Psi == 1)
    {
      if(run_type != "variational")
	{
	  clog << "Warning: attempting to optimize the wavefunction for"
	       << " run_type = " << run_type << "!" << endl;
	  //why not optimize using diffusion walkers?
	  //clog << "Setting run_type = \"variational\"" << endl;
	  //run_type = "variational";
	}

      /*
	Just check to be sure the parameter is valid so that the user
	doesn't have to wait until the end of the run to find out.
      */
      QMCLineSearchStepLengthSelectionAlgorithm *stepAlg =
	QMCLineSearchStepLengthSelectionFactory::factory(line_search_step_length);
      
      delete stepAlg;
      stepAlg = 0;

      if(optimize_Psi_method == "analytical_energy_variance" ||
	 optimize_Psi_method == "automatic")
	{
	  needPrintedConfigs = false;

	  if(optimize_Psi_criteria != "analytical_energy_variance" &&
	     optimize_Psi_criteria != "generalized_eigenvector")
	    {
	      /*
		Using analytical_energy_variance for method currently means we have
		"analytical" means of determining the hessian. This requires the analytical
		measurements of the gradients, which are currently only available if
		optimize_Psi_criteria = "analytical_energy_variance"

		As for the converse, there's no problem with using analytical gradients with
		other hessians (i.e. a different optimize_Psi_method). For example, Steepest_Descent.
	      */
	      clog << "Warning: if optimize_Psi_method = " << optimize_Psi_method
		   << " then optimize_Psi_criteria must = \"analytical_energy_variance\"" << endl;
	      clog << "Warning: setting optimize_Psi_criteria = \"analytical_energy_variance\"" << endl;
	      optimize_Psi_criteria = "analytical_energy_variance";
	    }

	  if(line_search_step_length == "Wolfe" ||
	     line_search_step_length == "MikesBracketing")
	    {
	      /*
		I don't think this parameter combination would be interesting, so warn.
	      */
	      clog << "Warning: you probably want \"None\" as your line_search_step_length!" << endl;
	      needPrintedConfigs = true;
	    }	  

	  if(ksi < 0.0 || ksi > 1)
	    {
	      clog << "Error: bad value for ksi = " << ksi << endl;
	      return false;
	    }
	}
      else
	{
	  needPrintedConfigs = true;

	  if(optimize_Psi_method == "Steepest_Descent" &&
	     (optimize_Psi_criteria == "analytical_energy_variance" ||
	      optimize_Psi_criteria == "generalized_eigenvector"))
	    needPrintedConfigs = false;

	  if(optimize_Psi_method == "BFGSQuasiNewton" && 
	     optimization_max_iterations <= 1)
	    {
	      clog << "Warning: BFGSQuasiNewton needs optimization_max_iterations > " << optimization_max_iterations << endl;
	      clog << "Setting optimization_max_iterations = 2" << endl;
	      optimization_max_iterations = 2;
	    }
	}

      if(replace_electron_nucleus_cusps == 1 &&
	 optimize_Psi_method == "automatic")
	{
	  /*
	    jastrows don't necessarily conflict with cusp replacement. is it worth programming something
	    to verify no conflict?
	  */
	  //clog << "Notice: the we won't optimize the orbitals since we're using cusp replacement." << endl;
	}

      if(link_Orbital_parameters == 0)
	{
	  if(trial_function_type == "restricted")
	    {
	      clog << "Notice: optimization with link_Orbital_parameters == " << link_Orbital_parameters
		   << " will convert your\n"
		   << "   trial_function_type == \"" << trial_function_type << "\" into\n"
		   << "   trial_function_type == \"unrestricted\"" << endl;
	    }
	} else {
	  if(trial_function_type == "unrestricted")
	    {
	      clog << "Warning: link_Orbital_parameters == " << link_Orbital_parameters
		   << " has no effect on a\n"
		   << "   trial_function_type == \"" << trial_function_type << "\"\n";
	      //Obviously you can't link the orbitals if they're unlinked to begin with
	    }
	}

    } else {
      //we're not optimizing
      needPrintedConfigs = false;
    }
  
  if(needPrintedConfigs && print_configs == 0)
    {
      clog << "Warning: based on your parameter choices, you need to have print_configs == 1."
	   << " It will be set to 1 now." << endl;
      print_configs = 1;
    }
  else if(!needPrintedConfigs && print_configs == 1)
    {
      clog << "Warning: print_configs == 1 is very expensive, and none of your parameter choices required it." << endl;
    }

  if(print_configs == 1)
    {
      clog << "Using " << temp_dir << " as scratch space for writing config files." << endl;
    }
  else if(optimize_Psi == 1)
    {
      calculate_Derivatives = 1;
    }

  if(nuclear_derivatives != "none" && Natoms == 0)
    {
      clog << "Warning: Setting nuclear_derivatives to \"none\" since Natoms = "
	   << Natoms << endl;
      nuclear_derivatives = "none";
    }

  if(sampling_method == "umrigar93_importance_sampling" && Natoms == 0)
    {
      clog << "Warning: Setting sampling_method to \"no_importance_sampling\" since " <<
	"umrigar93_importance_sampling requires Natoms > 0\n" << endl;
      sampling_method = "no_importance_sampling";
    }

  if(walkers_per_pass > number_of_walkers)
    {
      clog << "Warning: you requested more walkers_per_pass (" << walkers_per_pass
	   << ") than walkers (" << number_of_walkers << "), resetting walkers_per_pass." << endl;
      walkers_per_pass = number_of_walkers;
    }

  if(pseudo_gridLevel > 31){
    clog << "Error: our Lebedev-Laikov grids are only defined up to pseudo_gridLevel = 31.\n";
    clog << "You requested pseudo_gridLevel = " << pseudo_gridLevel << endl;
    return false;
  } else {

  }

#ifdef QMC_GPU
  if( getNumGPUWalkers() == 0)
  {
    clog << "A GPU run with no GPU walkers? Please set the gpu_walkers_per_pass parameter.\n";
    return false;
  }

  if( getNumGPUWalkers() + walkers_per_pass != number_of_walkers)
  {
    //clog << "Setting all walkers to be GPU walkers.\n";
    gpu_walkers_per_pass = -1;
  }

  cout << "Running GPUQMC with " << walkers_per_pass << " walkers per pass, " <<
    getNumGPUWalkers() << " of them on the GPU\n";

#endif

  if(use_surfer && iseed == 0)
    {
      iseed = 1000;
      clog << "Warning: you want a fixed seed for surfing... setting to " << iseed << endl;

    }
  if(chip_and_mike_are_cool != "Yea_Baby!")
    {
      clog << "ERROR: Incorrect value for chip_and_mike_are_cool set" << endl;
      return false;
    }
  return true;
}



