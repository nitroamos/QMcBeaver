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

#include "QMCReadAndEvaluateConfigs.h"

QMCReadAndEvaluateConfigs::QMCReadAndEvaluateConfigs()
{}

QMCReadAndEvaluateConfigs::QMCReadAndEvaluateConfigs(QMCInput *In,
    int cfgsToSkip)
{
  initialize(In, cfgsToSkip);
}

void QMCReadAndEvaluateConfigs::initialize(QMCInput *In, int cfgsToSkip)
{
  Input = In;
  configsToSkip = cfgsToSkip;

  Jastrow.initialize(Input);

  Nelectrons = Input->WF.getNumberElectrons();
  Natoms     = Input->Molecule.getNumberAtoms();

  // allocate the empty R vector that contains electrons

  R.allocate(Nelectrons,3);
  R = 0.0;

  // allocate D2

  D2.allocate(Nelectrons,3);
}

// Calculate the properites from the configs for all the parameters in params
// using the root node
void QMCReadAndEvaluateConfigs::rootCalculateProperties(
  Array1D < Array1D<double> > & Params, Array1D<QMCProperties> & Properties)
{
  //Params holds all the vectors to be evaluated
  //Their scores be returned in the end

#ifdef PARALLEL
  // Root packs all of the different parameter sets into a large vector
  // This vector is sent to all processors where it is unpacked and
  // the function is evaluated locally

  // Broadcast 1 to signal the workers to execute workerCalculateProperties
  // 0 signals the workers to stop

  int WorkSignal = 1;
  MPI_Bcast(&WorkSignal,1,MPI_INT,0,MPI_COMM_WORLD);

  // Send the size of the Params to all of the nodes

  int ParamSize[2];
  ParamSize[0] = Params.dim1();
  ParamSize[1] = Params(0).dim1();
  MPI_Bcast(&ParamSize,2,MPI_INT,0,MPI_COMM_WORLD);

  int elements = Params.dim1()*Params(0).dim1();
  double *PackedParameters = new double[elements];

  // The root node packs all of the parameter sets into a vector

  for(int i=0;i<Params.dim1();i++)
    {
      for(int j=0;j<Params(i).dim1();j++)
        {
          PackedParameters[i*Params(0).dim1()+j] = Params(i)(j);
        }
    }

  // send PackedParameters to all cpu's

  MPI_Bcast(PackedParameters,elements,MPI_DOUBLE,0,MPI_COMM_WORLD);

  delete [] PackedParameters;

#endif

  Array1D < QMCProperties > local_properties(Params.dim1());
  Properties.allocate(Params.dim1());

  for(int i=0;i<Params.dim1();i++)
    {
      local_properties(i).zeroOut();
      Properties(i).zeroOut();
    }

  //this function does all the dirty work of reading in the configs
  //and analyzing them and puts the results in local_properties

  locally_CalculateProperties(Params,local_properties);

  //reduce on this properties

  MPI_reduce(local_properties,Properties);
}

// Calculate the properites from the configs for all the parameters in params
// using the root node
void QMCReadAndEvaluateConfigs::workerCalculateProperties()
{
  //Params holds all the vectors to be evaluated
  //Their scores be returned in the end

#ifdef PARALLEL
  // Root packs all of the different parameter sets into a large vector
  // This vector is sent to all processors where it is unpacked and
  // the function is evaluated locally

  // Receive the size of the Params from the root node and allocate structures

  int ParamSize[2];
  MPI_Bcast(ParamSize,2,MPI_INT,0,MPI_COMM_WORLD);

  Array1D< Array1D<double> > Params;
  Params.allocate(ParamSize[0]);

  for(int i=0; i<ParamSize[0]; i++)
    {
      Params(i).allocate(ParamSize[1]);
    }

  Array1D<QMCProperties> Properties(ParamSize[0]);

  int elements = Params.dim1()*Params(0).dim1();
  double *PackedParameters = new double[elements];

  // Everyone else MPI recv PackedParameters by participating
  //in the MPI_BCast

  MPI_Bcast(PackedParameters,elements,MPI_DOUBLE,0,MPI_COMM_WORLD);

  // Everyone unpacks PackedParameters

  for(int i=0;i<Params.dim1();i++)
    {
      for(int j=0;j<Params(i).dim1();j++)
        {
          Params(i)(j) = PackedParameters[i*Params(0).dim1()+j];
        }
    }

  delete [] PackedParameters;

  Array1D < QMCProperties > local_properties(Params.dim1());

  for(int i=0;i<Params.dim1();i++)
    {
      local_properties(i).zeroOut();
      Properties(i).zeroOut();
    }

  //this function does all the dirty work of reading in the configs
  //and analyzing them and puts the results in local_properties

  locally_CalculateProperties(Params,local_properties);

  //reduce on this properties

  MPI_reduce(local_properties,Properties);
#endif
}

// Calculate the properites from the configs for all the parameters in params
// on the current node
void QMCReadAndEvaluateConfigs::locally_CalculateProperties(
  Array1D < Array1D<double> > &Params, Array1D<QMCProperties> & Properties)
{  
  Input->outputer.open(Input->flags.config_file_name,false);

  //zero out the properties
  Properties.allocate(Params.dim1());
  for(int j=0;j<Properties.dim1();j++)
    {
      Properties(j).zeroOut();
    }

  // Skip the appropriate number of configs from the beginning of the file.

  for (int i=0; i<configsToSkip; i++)
    {
      Input->outputer.readCorrelatedSamplingConfiguration(R,D1,D2,lnJ,PE);
    }

  // read configurations until the file is empty

  while( !Input->outputer.eof() )
    {
      // read the next configuration

      Input->outputer.readCorrelatedSamplingConfiguration(R,D1,D2,lnJ,PE);
      
      // On some computers, it seems this class was reading one extra config
      // in, assiging zeros to all the parameters. Why was it doing this?
      // Probably some compilers had a different definition of eof() or
      // something like that.  Anyway, assuming that the check for lnJ==0
      // identifies these (it does on QSC) then this fix should work.

      if (lnJ != 0)
        {
          // loop over the input configurations and calculate the corresponding
          // properties

          for(int i=0; i<Params.dim1(); i++)
            {
              AddNewConfigToProperites(Params(i),Properties(i));
            }
        }
    }
  Input->outputer.close();
}

// given a set of parameters perform the necessary calculations and
// add the results to the properties.
void QMCReadAndEvaluateConfigs::AddNewConfigToProperites(
  Array1D<double> &Params,QMCProperties &Properties)
{
  // Calculate the jastrow values

  Input->JP.setParameterVector( Params );

  //there is a much better way to do this...
  Array1D< Array2D<double>* > temp;
  temp.allocate(1);
  temp(0) = &R;
  Jastrow.evaluate(temp,1,0);

  double E_Local;
  double logWeight;

  bool Am_I_Valid = true;

  if( Am_I_Valid == true )
    {
      // If the Jastrow is nonsingular perform the calculation as normal

      // calculate the local energy and weight for this configuration

      E_Local = calc_E_Local_current();
      logWeight  = calc_log_weight_current();

      // Limit the size of E_Local or the weights

      if(E_Local > MAXIMUM_ENERGY_VALUE)
        {
          E_Local = MAXIMUM_ENERGY_VALUE;
        }

      if(logWeight > MAXIMUM_LOG_WEIGHT_VALUE)
        {
          logWeight = MAXIMUM_LOG_WEIGHT_VALUE;
        }
    }
  else
    {
      // If the Jastrow is singular provide default values for the
      // local energy and weight;

      E_Local = MAXIMUM_ENERGY_VALUE;
      logWeight  = MAXIMUM_LOG_WEIGHT_VALUE;
    }

  // Place the results into the properties

  Properties.energy.newSample(E_Local,exp(logWeight));
  Properties.logWeights.newSample(logWeight,1.0);
}

// Calculate the local energy of the current configuration with the currently
// calculated jastrow
double QMCReadAndEvaluateConfigs::calc_E_Local_current()
{
  double E_Local_current = 0.0;

  //calc Grad_sum_u_current

  Array2D<double> * Grad_sum_u_current = Jastrow.getGradientLnJastrow(0);

  //calc Lap_sum_u_current

  double Lap_sum_u_current = Jastrow.getLaplacianLnJastrow(0);

  //this is only for clarity and can be optimized
  //calc Grad_PsiRatio_current_without_Jastrow

  Array2D<double> Grad_PsiRatio_current_without_Jastrow(Nelectrons,3);

  for(int i=0; i<Nelectrons; i++)
    {
      for(int j=0; j<3; j++)
        {
          Grad_PsiRatio_current_without_Jastrow(i,j)=D2(i,j);
        }
    }

  //look at QMCFunctions.calculate_Grad_PsiRatio_current()
  //calc Grad_PsiRatio_current

  Array2D<double> Grad_PsiRatio_current(Nelectrons,3);

  for(int i=0; i<Nelectrons; i++)
    {
      for(int j=0; j<3; j++)
        {
          Grad_PsiRatio_current(i,j)=
            Grad_PsiRatio_current_without_Jastrow(i,j) +(*Grad_sum_u_current)(i,j);
        }
    }

  //calc alpha_beta_sum_of_Lap_PsiRatio_current

  double alpha_beta_sum_of_Lap_PsiRatio_current=D1;

  //look at QMCFunctions.calculate_Laplacian_PsiRatio_current()
  //calc Laplacian_PsiRatio_current

  double Laplacian_PsiRatio_current = alpha_beta_sum_of_Lap_PsiRatio_current +
                                      Lap_sum_u_current;

  for(int i=0; i<Nelectrons; i++)
    {
      for(int j=0; j<3; j++)
        {
          Laplacian_PsiRatio_current += (*Grad_sum_u_current)(i,j) *
                                        (2*Grad_PsiRatio_current(i,j)-(*Grad_sum_u_current)(i,j));
        }
    }

  //look at QMCFunctions.calculate_E_Local_current()
  //calc E_Local_current

  E_Local_current = -0.5 * Laplacian_PsiRatio_current + PE;
  return E_Local_current;
}

// Calculate the weight of the current configuration with the currently
// calculated jastrow
double QMCReadAndEvaluateConfigs::calc_log_weight_current()
{
  double logweight = 2*(Jastrow.getLnJastrow(0) - lnJ);
  return logweight;
}

// perform an mpi reduce operation on the properties collected on each
// processor
void QMCReadAndEvaluateConfigs::MPI_reduce(
  Array1D <QMCProperties> &local_Properties,
  Array1D < QMCProperties> &global_Properties)
{
  global_Properties.allocate(local_Properties.dim1());

#ifdef PARALLEL

  MPI_Reduce(local_Properties.array(),global_Properties.array(),
	     local_Properties.dim1(),QMCProperties::MPI_TYPE,
	     QMCProperties::MPI_REDUCE,0,MPI_COMM_WORLD);

#else

  for(int i=0;i<local_Properties.dim1();i++)
    {
      global_Properties(i) = local_Properties(i);
    }
#endif
}


