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

QMCFlags::QMCFlags()
{
  my_rank = 0;
  nprocs  = 1;
}

void QMCFlags::read_flags(string InFileName)
{
  ifstream input_file(InFileName.c_str());
  string temp_string;

  if(!input_file)
   {
      cerr << "ERROR: Can't open input " << InFileName.c_str() << endl;
#ifdef _WIN32
		getchar();
#endif  
      exit(1);
    }

  //***** Default Flag Values ********
  walker_initialization_method = "mikes_jacked_initialization";
  walker_initialization_combinations = 3;
  temp_dir = "/temp1/";
  parallelization_method = "manager_worker";
  iseed  = -5135696;
  sampling_method = "importance_sampling";
  QF_modification_type = "none";
  energy_modification_type = "none";
  umrigar93_equalelectrons_parameter = 0.5;  // in (0,1] 
  walker_reweighting_method = "umrigar93_probability_weighted";
  branching_method = "nonunit_weight_branching";
  branching_threshold = 2.0;
  fusion_threshold = 0.5;
  old_walker_acceptance_parameter = 50;
  zero_out_checkpoint_statistics = 0;
  Ndeterminants = 1;
  dt     = 0.001;
  dt_equilibration               = 0.02;

  use_basis_function_interpolation = 0;
  number_basis_function_interpolation_grid_points = 1000;
  basis_function_interpolation_first_point = 1e-10;

  equilibration_steps            = 10000;
  equilibration_function         = "ramp";
  CKAnnealingEquilibration1_parameter = 500;
  use_equilibration_array        = 0;

  correct_population_size_bias   = 0;

  desired_convergence            = 0.0;
  max_time_steps                 = 1000000;
  number_of_walkers              = 1;
  population_control_parameter   = 1.0;

  output_interval                = 1000;
  mpireduce_interval             = 1000;
  mpipoll_interval               = 1;
  checkpoint_interval            = 1000;
  checkpoint                     = 1;
  use_available_checkpoints      = 0;

  print_transient_properties     = 1;
  print_transient_properties_interval = 10000;
  print_configs                  = 0;
  print_config_frequency         = 50;

  optimize_Psi                   = 0;
  max_optimize_Psi_steps         = 10;
  optimize_Psi_barrier_parameter = 1.0;
  optimize_Psi_criteria          = "umrigar88";
  optimize_Psi_method            = "BFGSQuasiNewton";  
  numerical_derivative_surface   = "umrigar88";
  line_search_step_length = "Wolfe";

  optimization_max_iterations = 100;
  optimization_error_tolerance = 0.001;

  ck_genetic_algorithm_1_population_size = 30;
  ck_genetic_algorithm_1_mutation_rate = 0.2;
  ck_genetic_algorithm_1_initial_distribution_deviation = 1.0;

  singularity_penalty_function_parameter = 1.0e-6;

  link_Jastrow_parameters        = 0;
  equilibrate_every_opt_step     = 0;
  equilibrate_first_opt_step     = 1;
  write_all_energies_out         = 0;
  chip_and_mike_are_cool         = "false";

  //**********************************

  input_file >> temp_string;
  while(temp_string != "&flags")
    {
      input_file >> temp_string;
    }

  input_file >> temp_string;

  while((temp_string != "&") && (input_file.eof() != 1))
    {
      if(temp_string == "run_type")
        {
          input_file >> run_type;
	  if(run_type != "variational" && run_type != "diffusion") 
	    {
	      cerr << "ERROR: Unknown run_type: " << run_type << endl;
	      exit(1);
	    }
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
          input_file >> temp_string;
          iseed = atoi(temp_string.c_str());
          if (iseed > 0) iseed *= -1;
        }
      else if(temp_string == "sampling_method")
        {
          input_file >> sampling_method;
        }
      else if(temp_string == "QF_modification_type")
        {
          input_file >> QF_modification_type;
        }
      else if(temp_string == "energy_modification_type")
        {
          input_file >> energy_modification_type;
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
      else if(temp_string == "old_walker_acceptance_parameter")
        {
          input_file >> temp_string;
          old_walker_acceptance_parameter = atoi(temp_string.c_str());
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
      else if(temp_string == "desired_convergence")
        {
          input_file >> temp_string;
          desired_convergence = atof(temp_string.c_str());
        }
      else if(temp_string == "max_time_steps")
        {
          input_file >> temp_string;
          max_time_steps = atoi(temp_string.c_str());
        }
      else if(temp_string == "number_of_walkers")
        {
          input_file >> temp_string;
          number_of_walkers = atoi(temp_string.c_str());
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
      else if(temp_string == "chip_and_mike_are_cool")
	{
	  input_file >> chip_and_mike_are_cool;
	}
      else
	{
	  cerr << "ERROR: Unknown input flag: " << temp_string << endl;
	  exit(1);
	}
      input_file >> temp_string;
    }
  input_file.close();

  if(run_type == "" )
    {
      cerr << "ERROR: run_type not set!" << endl;
      exit(1);
    }

  if(mpireduce_interval < output_interval)
    {
      cerr << "ERROR: mpireduce_interval < output_interval!" << endl;
      exit(1);
    }

  if(optimize_Psi == 1 && run_type != "variational")
    {
      cerr << "ERROR: attempting to optimize the wavefunction for"
	   << " run_type = " << run_type << "!" << endl;
    }

  set_filenames(InFileName);

  dt_effective = dt;
  dt_run = dt;
  dt = dt_equilibration;
  number_of_walkers_initial = number_of_walkers;
  energy_estimated = energy_trial;
  energy_estimated_original = energy_estimated;

  if(chip_and_mike_are_cool != "Yea_Baby!")
    {
      cerr << "ERROR: Incorrect value for chip_and_mike_are_cool set" << endl;
      exit(1);
    }
}

void QMCFlags::set_filenames(string runfile)
{
  string tempstring;
  string basename;
  int restart_number = 0;

  input_file_name   = runfile;
  string file_name  = input_file_name.substr(0,input_file_name.size()-5);
  output_file_name  = file_name + ".qmc";
  results_file_name = file_name + ".rslts";
  base_file_name    = file_name;

  char my_rank_string[32];
#ifdef _WIN32
	  _snprintf( my_rank_string, 32, "%d", my_rank );
#else    
	  snprintf( my_rank_string, 32, "%d", my_rank );
#endif
  config_file_name = temp_dir + base_file_name + "."+my_rank_string + ".cfgs";

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
 strm << "&flags" << endl;
 strm << "run_type\n " << flags.run_type << endl;
 strm << "temp_dir\n " << flags.temp_dir << endl;
 strm << "parallelization_method\n " << flags.parallelization_method << endl;
 strm << "walker_initialization_method\n " 
      << flags.walker_initialization_method << endl;
 strm << "walker_initialization_combinations\n " 
      << flags.walker_initialization_combinations << endl;
 if (flags.iseed > 0) strm << "iseed\n " << -1*flags.iseed << endl;
 else if (flags.iseed <= 0) strm << "iseed\n " << flags.iseed << endl;
 strm << "sampling_method\n " << flags.sampling_method << endl;
 strm << "QF_modification_type\n " << flags.QF_modification_type << endl;
 strm << "energy_modification_type\n " << flags.energy_modification_type 
      << endl;
 strm << "umrigar93_equalelectrons_parameter\n " 
      << flags.umrigar93_equalelectrons_parameter << endl;
 strm << "walker_reweighting_method\n " << flags.walker_reweighting_method 
      << endl;
 strm << "branching_method\n " << flags.branching_method << endl;
 strm << "branching_threshold\n " << flags.branching_threshold << endl;
 strm << "old_walker_acceptance_parameter\n " 
      << flags.old_walker_acceptance_parameter << endl;
 strm << "use_basis_function_interpolation\n " 
      << flags.use_basis_function_interpolation << endl;
 strm << "number_basis_function_interpolation_grid_points\n " 
      << flags.number_basis_function_interpolation_grid_points << endl;
 strm << "basis_function_interpolation_first_point\n " 
      << flags.basis_function_interpolation_first_point << endl;
 strm << "fusion_threshold\n " << flags.fusion_threshold << endl;
 strm << "dt\n " << flags.dt_run << endl;
 strm << "desired_convergence\n " << flags.desired_convergence << endl;
 strm << "max_time_steps\n " << flags.max_time_steps << endl;
 strm << "dt_equilibration\n " << flags.dt_equilibration << endl;
 strm << "equilibration_steps\n " << flags.equilibration_steps << endl; 
 strm << "equilibration_function\n " << flags.equilibration_function << endl;
 strm << "CKAnnealingEquilibration1_parameter\n " 
      << flags.CKAnnealingEquilibration1_parameter << endl;
 strm << "use_equilibration_array\n " << flags.use_equilibration_array << endl;
 strm << "number_of_walkers\n " << flags.number_of_walkers << endl;
 strm << "output_interval\n " << flags.output_interval << endl;
 strm << "mpireduce_interval\n " << flags.mpireduce_interval << endl;
 strm << "mpipoll_interval\n " << flags.mpipoll_interval << endl;
 strm << "checkpoint_interval\n " << flags.checkpoint_interval << endl;
 strm << "checkpoint\n " << flags.checkpoint << endl;
 strm << "use_available_checkpoints\n " << flags.use_available_checkpoints 
      << endl;
 strm << "atoms\n " << flags.Natoms << endl;
 strm << "charge\n " << flags.charge << endl;
 strm << "norbitals\n " << flags.Norbitals << endl;
 strm << "nbasisfunc\n " << flags.Nbasisfunc << endl;
 strm << "ndeterminants\n " << flags.Ndeterminants << endl;
 strm << "energy\n " << flags.energy_trial << endl;
 strm << "correct_population_size_bias\n " 
      << flags.correct_population_size_bias << endl;
 strm << "print_transient_properties\n " << flags.print_transient_properties 
      << endl;
 strm << "print_transient_properties_interval\n " 
      << flags.print_transient_properties_interval << endl;
 strm << "print_configs\n " << flags.print_configs << endl;
 strm << "print_config_frequency\n " << flags.print_config_frequency << endl;
 strm << "optimize_Psi\n " << flags.optimize_Psi << endl;
 strm << "max_optimize_Psi_steps\n " << flags.max_optimize_Psi_steps << endl;
 strm << "optimize_Psi_criteria\n " << flags.optimize_Psi_criteria << endl;
 strm << "optimize_Psi_method\n " << flags.optimize_Psi_method << endl;
 strm << "singularity_penalty_function_parameter\n " 
      << flags.singularity_penalty_function_parameter << endl;
 strm << "optimize_Psi_barrier_parameter\n " 
      << flags.optimize_Psi_barrier_parameter << endl;
 strm << "numerical_derivative_surface\n " 
      << flags.numerical_derivative_surface << endl;
 strm << "line_search_step_length\n " << flags.line_search_step_length << endl;
 strm << "optimization_max_iterations\n " 
      << flags.optimization_max_iterations << endl;
 strm << "optimization_error_tolerance\n " 
      << flags.optimization_error_tolerance << endl;
 strm << "ck_genetic_algorithm_1_population_size\n " 
      << flags.ck_genetic_algorithm_1_population_size << endl;
 strm << "ck_genetic_algorithm_1_mutation_rate\n " 
      << flags.ck_genetic_algorithm_1_mutation_rate << endl;
 strm << "ck_genetic_algorithm_1_initial_distribution_deviation\n "
      << flags.ck_genetic_algorithm_1_initial_distribution_deviation << endl;
 strm << "link_Jastrow_parameters\n " << flags.link_Jastrow_parameters << endl;
 strm << "equilibrate_first_opt_step\n " << flags.equilibrate_first_opt_step 
      << endl;
 strm << "equilibrate_every_opt_step\n " << flags.equilibrate_every_opt_step 
      << endl;
 strm << "population_control_parameter\n " 
      << flags.population_control_parameter << endl;
 strm << "write_all_energies_out\n " << flags.write_all_energies_out << endl;
 strm << "zero_out_checkpoint_statistics\n " 
      << flags.zero_out_checkpoint_statistics << endl;
 strm << "chip_and_mike_are_cool\n " << flags.chip_and_mike_are_cool << endl;
 strm << "& " << endl;
 return strm;
}





