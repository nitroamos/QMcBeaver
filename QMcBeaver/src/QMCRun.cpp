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

#include "QMCRun.h"
#include <map>

QMCRun::~QMCRun()
{
  delete QMF;
  QMF = 0;
}

QMCRun::QMCRun()
{
  populationSizeBiasCorrectionFactor = 1.0;
}

void QMCRun::propagateWalkers(bool writeConfigs, int iteration)
{
  int count = 0;
  int index = 0;
  int wpp = Input->flags.walkers_per_pass;
  
  /*The point here is to collect WALKERS_PER_PASS amount of walkers
    to process at once. Once we have that many (or the last remaining), we 
    finally run QMF->evaluate which fills up the dataPointers data structure 
    with values.  The actual data for dataPointers is stored in each walker, so
    initializePropagation asks for a pointer to it. The collection of pointers 
    is then passed around to everybody.
    The advantage to filling the array dynamically is that we don't have to 
    worry about branching -- walkers being deleted and created.
  */
  for(list<QMCWalker>::iterator wp=wlist.begin();wp!=wlist.end();++wp)
    {
      swTimers(0).start();
      wp->initializePropagation(dataPointers(index),rPointers(index), iteration);
      swTimers(0).stop();

      count++;
      index = count%wpp;
      if(index == 0 || count == (int)(wlist.size()))
	{
	  int num = index==0?wpp:index;
	  QMF->evaluate(dataPointers,rPointers,num);

	  if(globalInput.cs_Parameters.dim1() > 1)
	    {
	      //we don't need to do this if we're equilibrating, since we're just going
	      //to throw the data away
	      swTimers(1).start();
	      if(iteration >= 0)
		QMF->calculate_CorrelatedSampling(dataPointers,rPointers,num);
	      swTimers(1).lap(num);
	    }
	  else
	    {
	      for(int i=0; i<num; i++)
		{
		  dataPointers(i)->cs_LocalEnergy.deallocate();
		  dataPointers(i)->cs_Weights.deallocate();
		}
	    }

	}
    }

  /*At this point, all of the dataPointers have been filled, so we
    can tell each walker to finish processing each move. Note this will
    give different answers than before (when all processing was done before
    the next walker was evaluated) because random numbers are drawn in a
    different order.
   */
  swTimers(0).start();
  for(list<QMCWalker>::iterator wp=wlist.begin();wp!=wlist.end();++wp)
      wp->processPropagation(*QMF, writeConfigs);
  swTimers(0).lap(getNumberOfWalkers());
}

void QMCRun::branchWalkers()
{
  if( Input->flags.run_type == "variational" )
    {
      // No branching in VMC
    }
  else if( Input->flags.branching_method == "non_branching" )
    {
      // Non branching DMC -- The weights are varied instead of branching
    }
  else if( Input->flags.branching_method == "unit_weight_branching" )
    {
      // This is the method described in Lester and Hammonds book
      // It is Chip Kent's (my) experience, that this method
      // often has an exponentially growing or shrinking population
      // and has a large time step error
      
      unitWeightBranching();
    }
  else if( Input->flags.branching_method == "nonunit_weight_branching" )
    {
      // This method is described in Umrigar, Nightingale, and Runge
      // JCP 99(4) 2865; 1993 and elsewhere.  A walker is reweighted
      // and is divided if its weight exceeds a threshold and is fused
      // with another walker if its weight falls below a threshold.
      
      nonunitWeightBranching();
    }
  else if( Input->flags.branching_method == "ack_reconfiguration" )
    {
      
      ack_reconfiguration();
    }
  else
    {
      cerr << "ERROR: Unknown branching method!" << endl;
      exit(1);
    }
}

void QMCRun::zeroOut()
{
  for( list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end(); ++wp )
    wp->resetFutureWalking();
    
  if (Input->flags.use_equilibration_array == 1)
    {
      EquilibrationArray.zeroOut();
    } else {
      Properties.zeroOut();
      fwProperties.zeroOut();
    }
}

void QMCRun::initializeFunction()
{
  QMF = QMCFunctionsFactory::functionsFactory(Input, Input->flags.trial_function_type);

  if(Input->flags.trial_function_type ==   "restricted" ||
     Input->flags.trial_function_type == "unrestricted"){
    QMCSCFJastrow * qmfHFJ = static_cast<QMCSCFJastrow*>(QMF);
    qmfHFJ->initialize(Input,&HartreeFock);
  }

  //testPseudoPotential();
}

void QMCRun::testPseudoPotential()
{
  QMCInitializeWalker * IW =
    QMCInitializeWalkerFactory::initializeWalkerFactory(Input,
        Input->flags.walker_initialization_method);
  
  QMCWalker w;
  w.initialize(Input);
  Array2D<double> R = IW->initializeWalkerPosition();
  w.setR(R);

  //Pick the closest electron for our testing
  int elec = 0;
  for(int e=0; e<R.dim1(); e++){
    if(w.getWalkerData()->riA(elec,0) > w.getWalkerData()->riA(e,0)){
      elec = e;
    }
  }
  printf("\n\n\nStarting pseudopotential ratio check...\n");
  printf("       %2s:%1s %20s %20s/%20s %20s\n","gr","N","r","PsiP(gr)","Psi","ratio");
  QMCPotential_Energy::printElec = elec;
  QMF->evaluate(*w.getR(),*w.getWalkerData());
  QMCPotential_Energy::printElec = -2;
  printf("\n");

  double Psi = w.getWalkerData()->psi;
  Array1D<double> origR(globalInput.Molecule.getNumberAtoms());
  for(int nuc=0; nuc<globalInput.Molecule.getNumberAtoms(); nuc++)
    origR(nuc) = w.getWalkerData()->riA(elec,nuc);

  for(int nuc=0; nuc<globalInput.Molecule.getNumberAtoms(); nuc++){
    double r = origR(nuc);
    int whichRG = 39; //remember, QMCPotential_Energy chooses this randomly...
    Array2D<double> grid = globalInput.Molecule.getGrid(nuc,whichRG,r,true);

    Array1D<double> PsiP(grid.dim1());
    PsiP = 0.0;
    for(int gr=0; gr<grid.dim1(); gr++){

      for(int xyz=0; xyz<3; xyz++)
	R(elec,xyz) = grid(gr,xyz);

      //printf("x %20.10e y %20.10e z %20.10e\n",R(elec,0),R(elec,1),R(elec,2));
      
      w.setR(R); //Make sure the function accepts zeros as coordinates
      QMF->evaluate(*w.getR(),*w.getWalkerData());
      PsiP(gr) = w.getWalkerData()->psi;
      printf("QMCRUN %2i:%i %20.10f %20.10e/%20.10e %20.10e\n",gr,nuc,r,PsiP(gr),Psi,PsiP(gr)/Psi);
    }
  }
  exit(0);
}

void QMCRun::initialize(QMCInput *INPUT)
{
  Input = INPUT;

  initializeFunction();

  swTimers.allocate(4);
  swTimers(0).reset("Walker Propagation");
  swTimers(1).reset("Correlated Sampling");
  swTimers(2).reset("Calculate Observables");
  swTimers(3).reset("Entire Step");

  int wpp = Input->flags.walkers_per_pass;
  dataPointers.allocate(wpp);
  rPointers.allocate(wpp);

  if (Input->flags.calculate_bf_density == 1)
    {
      bool calcDensity = true;
      fwTimeStepProperties.setCalcDensity(calcDensity,
					  Input->WF.getNumberBasisFunctions());
      if (Input->flags.use_equilibration_array == 1)
        EquilibrationArray.setCalcDensity(calcDensity,
                                          Input->WF.getNumberBasisFunctions());
      else
        fwProperties.setCalcDensity(calcDensity,
				    Input->WF.getNumberBasisFunctions());
    }
 
  if(Input->flags.nuclear_derivatives != "none")
    {
      int dim1, dim2;
      if(Input->flags.nuclear_derivatives != "bin_force_density")
      {
        dim1 = Input->Molecule.getNumberAtoms();
        dim2 = 3;
      } else {
        dim1 = QMCNuclearForces::getNumBins();
        dim2 = 1;
      }
      
      bool calcForces = true;
      fwTimeStepProperties.setCalcForces(calcForces,dim1,dim2);
      if (Input->flags.use_equilibration_array == 1)
	{
	  EquilibrationArray.setCalcForces(calcForces,dim1,dim2);
	} else {
	  fwProperties.setCalcForces(calcForces,dim1,dim2);    
	}
    }
    
  if (Input->flags.use_equilibration_array == 1)
    {
      EquilibrationArray.zeroOut();
    } else {
      Properties.zeroOut();
      fwProperties.zeroOut();
    }

  int nalpha = Input->WF.getNumberElectrons(true);
  int nbeta = Input->WF.getNumberElectrons(false);

  if (Input->flags.write_electron_densities == 1)
    {
      // The pair density for each pair of particles is recorded in a 
      // histogram.  There are histograms for opposite and parallel spin
      // electrons, and for alpha and beta electrons and each unique nucleus in
      // the molecule.  The maximum distance is 20 divided by the largest
      // atomic charge in the molecule, and there are 5,000 bins.
      
      int max_Z = 0;
      
      if (Input->flags.max_pair_distance < 0)
	{
	  for (int i=0; i<Input->Molecule.getNumberAtoms(); i++)
	    if ( Input->Molecule.Z(i) > max_Z )
	      max_Z = Input->Molecule.Z(i);
          
	  max_pair_distance = 20.0/max_Z;
	}
      else
	max_pair_distance = Input->flags.max_pair_distance;

      dr = max_pair_distance/5000;
      
      // Histograms for parallel and opposite spin electron densities.

      if (nalpha > 1 || nbeta > 1)
	{
	  pllSpinHistogram.allocate(5000);
	  pllSpinHistogram = 0.0;
	}
	  
      if (nalpha > 0 && nbeta > 0)
	{
	  oppSpinHistogram.allocate(5000);
	  oppSpinHistogram = 0.0;
	}

      // Histograms for the one electron densities.
      int nucleiTypes = Input->Molecule.NucleiTypes.dim1();

      if (nalpha > 0)
	{
	  alphaHistograms.allocate(nucleiTypes);
	  for (int k=0; k<nucleiTypes; k++)
	    {
	      alphaHistograms(k).allocate(5000);
	      alphaHistograms(k) = 0.0;
	    }
	}

      if (nbeta > 0)
	{
	  betaHistograms.allocate(nucleiTypes);
	  for (int k=0; k<nucleiTypes; k++)
	    {
	      betaHistograms(k).allocate(5000);
	      betaHistograms(k) = 0.0;
	    }
	}
    }

  if (nalpha > 1 || nbeta > 1)
    {
      if (Input->flags.writePllxCorrelationDiagram == 1)
	{
	  if (Input->flags.pllxCorrelationDiagramMin >= 
	      Input->flags.pllxCorrelationDiagramMax)
	    {
	      cerr << "ERROR: pllxCorrelationDiagramMin = "
		   << Input->flags.pllxCorrelationDiagramMin << ", "
		   << "pllxCorrelationDiagramMax = "
		   << Input->flags.pllxCorrelationDiagramMax << endl;
	      exit(1);
	    }
	  pllxCorrelationDiagram.allocate(1000);
	  for (int i=0; i<pllxCorrelationDiagram.dim1(); i++)
	    {
	      pllxCorrelationDiagram(i).allocate(1000);
	      pllxCorrelationDiagram(i) = 0.0;
	    }
	}
      if (Input->flags.writePllyCorrelationDiagram == 1)
	{
	  if (Input->flags.pllyCorrelationDiagramMin >= 
	      Input->flags.pllyCorrelationDiagramMax)
	    {
	      cerr << "ERROR: pllyCorrelationDiagramMin = "
		   << Input->flags.pllyCorrelationDiagramMin << ", "
		   << "pllyCorrelationDiagramMax = "
		   << Input->flags.pllyCorrelationDiagramMax << endl;
	      exit(1);
	    }
	  pllyCorrelationDiagram.allocate(1000);
	  for (int i=0; i<pllyCorrelationDiagram.dim1(); i++)
	    {
	      pllyCorrelationDiagram(i).allocate(1000);
	      pllyCorrelationDiagram(i) = 0.0;
	    }
	}
      if (Input->flags.writePllzCorrelationDiagram == 1)
	{
	  if (Input->flags.pllzCorrelationDiagramMin >= 
	      Input->flags.pllzCorrelationDiagramMax)
	    {
	      cerr << "ERROR: pllzCorrelationDiagramMin = "
		   << Input->flags.pllzCorrelationDiagramMin << ", "
		   << "pllzCorrelationDiagramMax = "
		   << Input->flags.pllzCorrelationDiagramMax << endl;
	      exit(1);
	    }
	  pllzCorrelationDiagram.allocate(1000);
	  for (int i=0; i<pllzCorrelationDiagram.dim1(); i++)
	    {
	      pllzCorrelationDiagram(i).allocate(1000);
	      pllzCorrelationDiagram(i) = 0.0;
	    }
	}
    }

  if (nalpha > 0 && nbeta > 0)
    {
      if (Input->flags.writeOppxCorrelationDiagram == 1)
	{
	  if (Input->flags.oppxCorrelationDiagramMin >= 
	      Input->flags.oppxCorrelationDiagramMax)
	    {
	      cerr << "ERROR: oppxCorrelationDiagramMin = "
		   << Input->flags.oppxCorrelationDiagramMin << ", "
		   << "oppxCorrelationDiagramMax = "
		   << Input->flags.oppxCorrelationDiagramMax << endl;
	      exit(1);
	    }
	  oppxCorrelationDiagram.allocate(1000);
	  for (int i=0; i<oppxCorrelationDiagram.dim1(); i++)
	    {
	      oppxCorrelationDiagram(i).allocate(1000);
	      oppxCorrelationDiagram(i) = 0.0;
	    }
	}
      if (Input->flags.writeOppyCorrelationDiagram == 1)
	{
	  if (Input->flags.oppyCorrelationDiagramMin >= 
	      Input->flags.oppyCorrelationDiagramMax)
	    {
	      cerr << "ERROR: oppyCorrelationDiagramMin = "
		   << Input->flags.oppyCorrelationDiagramMin << ", "
		   << "oppyCorrelationDiagramMax = "
		   << Input->flags.oppyCorrelationDiagramMax << endl;
	      exit(1);
	    }
	  oppyCorrelationDiagram.allocate(1000);
	  for (int i=0; i<oppyCorrelationDiagram.dim1(); i++)
	    {
	      oppyCorrelationDiagram(i).allocate(1000);
	      oppyCorrelationDiagram(i) = 0.0;
	    }
	}
      if (Input->flags.writeOppzCorrelationDiagram == 1)
	{
	  if (Input->flags.oppzCorrelationDiagramMin >= 
	      Input->flags.oppzCorrelationDiagramMax)
	    {
	      cerr << "ERROR: oppzCorrelationDiagramMin = "
		   << Input->flags.oppzCorrelationDiagramMin << ", "
		   << "oppzCorrelationDiagramMax = "
		   << Input->flags.oppzCorrelationDiagramMax << endl;
	      exit(1);
	    }
	  oppzCorrelationDiagram.allocate(1000);
	  for (int i=0; i<oppzCorrelationDiagram.dim1(); i++)
	    {
	      oppzCorrelationDiagram(i).allocate(1000);
	      oppzCorrelationDiagram(i) = 0.0;
	    }
	}
    }

  if (Input->flags.use_hf_potential == 1)
    HartreeFock.Initialize(Input);
}

void QMCRun::randomlyInitializeWalkers()
{
  Array2D<double> temp_R;
  int initialization_try = 0;

  QMCInitializeWalker * IW =
    QMCInitializeWalkerFactory::initializeWalkerFactory(Input,
        Input->flags.walker_initialization_method);

  wlist.clear();
  Input->flags.number_of_walkers = 0;

  for (int i=0; i<Input->flags.number_of_walkers_initial; i++)
    {
      QMCWalker w;
      w.initialize(Input);

      temp_R = IW->initializeWalkerPosition();
      int numtries = 0;
      while(w.setR(temp_R) == false){
	if(numtries > 10){
	  cerr << "Error: initial walker " << i << " electronic position is bad, and we've tried " << numtries << " times." << endl;
	  exit(0);
	} else {
	  cerr << "Error: initial walker " << i << " electronic position is bad, retrying..." << endl;
	}
	cerr.flush();
	temp_R = IW->initializeWalkerPosition();
	numtries++;
      }

      QMF->evaluate(*w.getR(),*w.getWalkerData());
  
      /*
	A lot of walkers are going to start off with a poor
	rel_diff. We only want to delete the ones that are practically
	singular, because the other poor ones will be taken care of.

	On the other hand, if we delete them here, we get a fresh walker,
	as opposed to deleting it in QMCWalker, when it will be replaced
	by a different, high weight, walker.
       */
      double rel_diff = fabs( (w.getWalkerData()->localEnergy -
			       Input->flags.energy_estimated_original)/
			      Input->flags.energy_estimated_original);
      initialization_try = 1;
      while( (w.isSingular() || rel_diff > 1.0 || w.getWalkerData()->localEnergy > 0) && initialization_try < 1000)
	{	  
	  cerr << "Reinitializing Walker " << setw(4) << i
	       << " with energy " << setprecision(5) << setw(11) << w.getWalkerData()->localEnergy
	       << ", with psi " << setprecision(5) << setw(11) << (double)w.getWalkerData()->psi
	       << ", rel_diff = " << rel_diff <<  "..." << endl;
	  cerr.flush();
	  temp_R = IW->initializeWalkerPosition();
	  w.setR(temp_R);
	  QMF->evaluate(*w.getR(),*w.getWalkerData());
	  
	  rel_diff = fabs( (w.getWalkerData()->localEnergy -
			    Input->flags.energy_estimated_original)/
			   Input->flags.energy_estimated_original);

	  initialization_try++;
	}

      if(initialization_try > 1){
	cerr << "             Accepted     "
	     << "             " << setprecision(5) << setw(11) << w.getWalkerData()->localEnergy
	     << ",          " << setprecision(5) << setw(11) << (double)w.getWalkerData()->psi
	     << ",            " << rel_diff << endl;
	cerr.flush();
      }

      int cutoff = 100;
      if( initialization_try > cutoff )
	{
	  cerr << "ERROR: " << initialization_try << " consecutive singular configurations while "
	       << "trying to initialize walker!" << endl;
	}
      
      w.newID();
      wlist.push_back(w);
      Input->flags.number_of_walkers++;
    }
  delete IW;
  IW = 0;
}

void QMCRun::calculateObservables()
{
  // Get the observables calculated in each walker
  // Pre block all of the statistics from one time step together
  
  timeStepProperties.zeroOut();
  fwTimeStepProperties.zeroOut();

  for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end();++wp)
    {    
      wp->calculateObservables( timeStepProperties );
      wp->calculateObservables( fwTimeStepProperties );
      wp->calculateDerivatives( fwProperties );
    }

  // Add the pre blocked data from this time step to the accumulated
  // statistics
  timeStepProperties.growthRate.newSample(growthRate,1.0);

  double totalWeights = getWeights() * getPopulationSizeBiasCorrectionFactor();

  if (Input->flags.use_equilibration_array == 1)
    {
      EquilibrationArray.newSample(&timeStepProperties, totalWeights,
				   getNumberOfWalkers());
    } else {
      Properties.newSample(&timeStepProperties,
			   totalWeights,
			   getNumberOfWalkers());
      fwProperties.newSample(&fwTimeStepProperties,
			     totalWeights,
			     getNumberOfWalkers());
    }
}

void QMCRun::writeEnergies(ostream& strm)
{
  // Get the energy calculated in each walker and write out
  
  for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end(); ++wp)
    strm << wp->getLocalEnergyEstimator() << "\t" << wp->getWeight() << endl;
}

void QMCRun::calculateElectronDensities()
{
  for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end(); ++wp)
    wp->calculateElectronDensities(max_pair_distance, dr, pllSpinHistogram,
			    oppSpinHistogram, alphaHistograms, betaHistograms);

  int nalpha = Input->WF.getNumberElectrons(true);
  int nbeta = Input->WF.getNumberElectrons(false);

  if (nalpha > 1 || nbeta > 1)
    {
      if (Input->flags.writePllxCorrelationDiagram == 1)
	for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end(); ++wp)
	  wp->calculatePllCorrelationDiagram(0, 
Input->flags.pllxCorrelationDiagramMin, Input->flags.pllxCorrelationDiagramMax,
					               pllxCorrelationDiagram);

      if (Input->flags.writePllyCorrelationDiagram == 1)
	for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end(); ++wp)
	  wp->calculatePllCorrelationDiagram(1, 
Input->flags.pllyCorrelationDiagramMin, Input->flags.pllyCorrelationDiagramMax,
					               pllyCorrelationDiagram);

      if (Input->flags.writePllzCorrelationDiagram == 1)
	for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end(); ++wp)
	  wp->calculatePllCorrelationDiagram(2, 
Input->flags.pllzCorrelationDiagramMin, Input->flags.pllzCorrelationDiagramMax,
					               pllzCorrelationDiagram);
    }

  if (nalpha > 0 && nbeta > 0)
    {
      if (Input->flags.writeOppxCorrelationDiagram == 1)
	for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end(); ++wp)
	  wp->calculateOppCorrelationDiagram(0, 
Input->flags.oppxCorrelationDiagramMin, Input->flags.oppxCorrelationDiagramMax,
					               oppxCorrelationDiagram);

      if (Input->flags.writeOppyCorrelationDiagram == 1)
	for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end(); ++wp)
	  wp->calculateOppCorrelationDiagram(1, 
Input->flags.oppyCorrelationDiagramMin, Input->flags.oppyCorrelationDiagramMax,
					               oppyCorrelationDiagram);

      if (Input->flags.writeOppzCorrelationDiagram == 1)
	for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end(); ++wp)
	  wp->calculateOppCorrelationDiagram(2, 
Input->flags.oppzCorrelationDiagramMin, Input->flags.oppzCorrelationDiagramMax,
					               oppzCorrelationDiagram);
    }
}

Array1D<double>* QMCRun::getPllSpinHistogram()
{
  return &pllSpinHistogram;
}

Array1D<double>* QMCRun::getOppSpinHistogram()
{
  return &oppSpinHistogram;
}

Array1D< Array1D<double> >* QMCRun::getAlphaHistograms()
{
  return &alphaHistograms;
}

Array1D< Array1D<double> >* QMCRun::getBetaHistograms()
{
  return &betaHistograms;
}

Array1D< Array1D<double> >* QMCRun::getPllxCorrelationDiagram()
{
  return &pllxCorrelationDiagram;
}

Array1D< Array1D<double> >* QMCRun::getPllyCorrelationDiagram()
{
  return &pllyCorrelationDiagram;
}

Array1D< Array1D<double> >* QMCRun::getPllzCorrelationDiagram()
{
  return &pllzCorrelationDiagram;
}

Array1D< Array1D<double> >* QMCRun::getOppxCorrelationDiagram()
{
  return &oppxCorrelationDiagram;
}

Array1D< Array1D<double> >* QMCRun::getOppyCorrelationDiagram()
{
  return &oppyCorrelationDiagram;
}

Array1D< Array1D<double> >* QMCRun::getOppzCorrelationDiagram()
{
  return &oppzCorrelationDiagram;
}

double QMCRun::getdr()
{
  return dr;
}

void QMCRun::unitWeightBranching()
{
  // Make a list of walkers to add and delete
  list<QMCWalker> WalkersToAdd;
  list<list<QMCWalker>::iterator> WalkersToDelete;

  for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end();++wp)
    {
      int times_to_branch = 0;
      if(wp->branchRecommended(true))
	times_to_branch = int(wp->getWeight() + ran.unidev())-1;
	
      // Set the walkers weight back to unity
      wp->setWeight(1.0);
      
      if( times_to_branch == 0 )
        {
          // Don't need to copy the walker
        }
      else if( times_to_branch < 0 )
        {
          // Mark the walker to delete
          
          WalkersToDelete.push_back( wp );
        }
      else
        {
          // Mark the walker to add times_to_branch times
          
          for(int i=0; i<times_to_branch; i++)
            {
              WalkersToAdd.push_back( *wp );
	      wp->branchID();
            }
        }
    }
    
  for( list< list<QMCWalker>::iterator >::iterator
       wtd=WalkersToDelete.begin(); wtd!=WalkersToDelete.end(); ++wtd)
    {
      wlist.erase( *wtd );
    }
    
  wlist.splice( wlist.end(), WalkersToAdd );
}

void QMCRun::nonunitWeightBranching()
{
  // Make a list of walkers to add and delete
  list<QMCWalker> WalkersToAdd;
  list<list<QMCWalker>::iterator> WalkersToDelete;
  
  list<QMCWalker>::iterator TempWalkerToDelete;
  bool isTempWalkerToDelete = false;
  
  for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end();++wp)
    {
      if( wp->getWeight() > Input->flags.branching_threshold )
        {
	  if(wp->branchRecommended(true))
	    {
	      // Split this walker into two; each with half the weight
	      wp->setWeight( wp->getWeight()/2.0 );
	      
	      WalkersToAdd.push_back( *wp );
	      wp->branchID();
	    }
        }
      else if( wp->getWeight() < 0.01 )
	{
	  WalkersToDelete.push_back(wp);
	}
      else if( wp->getWeight() < Input->flags.fusion_threshold )
        {
          // This walker needs to be fused with another
          
          if( !isTempWalkerToDelete )
            {
              // No walker to immediately fuse this with so wait for one
              
              TempWalkerToDelete = wp;
              isTempWalkerToDelete = true;
            }
          else
            {
              double weight1 = TempWalkerToDelete->getWeight();
              double weight2 = wp->getWeight();
              double weight3 = weight1 + weight2;
              
              if( ran.unidev() < weight1/weight3 )
                {
                  // Keep TempWalkerToDelete and delete other
                  
                  TempWalkerToDelete->setWeight( weight3 );
                  WalkersToDelete.push_back( wp );
                }
              else
                {
                  // opposite of above
                  
                  wp->setWeight( weight3 );
                  WalkersToDelete.push_back( TempWalkerToDelete );
                }
                
              isTempWalkerToDelete = false;
            }
        }
      else
        {
          // walker is ok and don't do anything to it
        }
    }
    
  for(list< list<QMCWalker>::iterator >::iterator wtd=WalkersToDelete.begin();
      wtd!=WalkersToDelete.end(); ++wtd)
    {
      wlist.erase( *wtd );
    }
    
  wlist.splice( wlist.end(), WalkersToAdd );
}

void QMCRun::ack_reconfiguration()
{
  // This is the new target average
  double aveW = getWeights()/getNumberOfWalkers();

  map<double,list<QMCWalker>::iterator> aCandidates;
  map<double,list<QMCWalker>::iterator> dCandidates;
  double Nreconfp  = 0.0;
  double Nreconfn  = 0.0;
  for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end();++wp)
    {
      double w = wp->getWeight()/aveW;
      double sw = fabs(w-1.0);
      if(w >= 1.0){
	if(wp->branchRecommended(true)){
	  Nreconfp += sw;
	  aCandidates[Nreconfp] = wp;
	}
      } else {
	Nreconfn += sw;
	dCandidates[Nreconfn] = wp;
      }
    }
  
  //Remove Nreconf randomly selected low weight walkers, proportionally to their weight
  int Nreconf = min((int)(Nreconfn + ran.unidev()), (int)dCandidates.size());
  int Nreconf1 = Nreconf;
  //Let's ensure that these bad walkers get deleted
  list<list<QMCWalker>::iterator> del;
  for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end();++wp)
    if( wp->getWeight() < 1e-5 && Nreconf > 0){
      del.push_back(wp);  
      Nreconf--;
    }
  int Nreconf2 = Nreconf;
  int count = 0;
  while(Nreconf > 0)
    {
      count++;
      double ndx = Nreconfn*ran.unidev();
      for(map<double,list<QMCWalker>::iterator>::iterator it=dCandidates.begin();
	  it != dCandidates.end(); it++)
	if(ndx < (*it).first)
	  {
	    //Can't delete the same walker twice
	    if(find(del.begin(),del.end(),(*it).second) == del.end())
	      {
		del.push_back((*it).second);
		Nreconf--;
	      }
	    break;
	  }
 
      if(count > 50){
	cout << "Error in count: Nreconf = " << Nreconf << endl;
	break;
      }
    }      
  
  //Now add enough walkers to bring us back to the initial number.
  int numToAdd = min((int)(Input->flags.number_of_walkers_initial-getNumberOfWalkers()+del.size()),
		     (int)aCandidates.size());
  list<QMCWalker> add;
  list<list<QMCWalker>::iterator> addList;
  count = 0;
  while(numToAdd > 0)
    {
      count++;
      double ndx = Nreconfp*ran.unidev();
      for(map<double,list<QMCWalker>::iterator>::iterator it=aCandidates.begin();
	  it != aCandidates.end(); it++)
	if(ndx < (*it).first)
	  {
	    if(find(addList.begin(),addList.end(),(*it).second) == addList.end())
	      {
		addList.push_back((*it).second);
		add.push_back(*(*it).second);
		(*it).second->branchID();
		numToAdd--;
	      }
	    break;
	  }
 
      if(count > 50){
	cout << "Error in count: numToAdd = " << numToAdd << endl;
	break;
      }
    }      


  /*
  for(int i=0; i<numToAdd; i++)
    {
      double ndx = Nreconfp*ran.unidev();
      for(map<double,list<QMCWalker>::iterator>::iterator it=aCandidates.begin();
	  it != aCandidates.end();
	  ++it)
	if(ndx < (*it).first)
	  {
	    add.push_back( *(*it).second );
	    (*it).second->branchID();	  
	    break;
	  }
    }
  */
  if(Input->flags.number_of_walkers_initial != getNumberOfWalkers() +(int)add.size() - (int)del.size())
    {
      cout << "ERROR: walkers didn't rebalance (" <<
	getNumberOfWalkers() << ")" << endl;
      cout << "Add Size: " << add.size() << endl;
      cout << "Del Size: " << del.size() << endl;
      cout << "numToAdd: " << numToAdd << endl;
      cout << "Nreconfp: " << setprecision(20) << Nreconfp << endl;
      cout << "Nreconfn: " << setprecision(20) << Nreconfn << endl;
      cout << "Nreconf1: " << Nreconf1 << endl;
      cout << "Nreconf2: " << Nreconf2 << endl;
      cout << "aveW:     " << setprecision(20) << aveW << endl;

      cout << "added: " << endl;
      for(list<QMCWalker>::iterator wp=add.begin(); wp!=add.end();++wp)
	cout << "w= " << wp->getWeight() << " sw = " << fabs(wp->getWeight()/aveW-1) << endl;

      cout << "deleted: " << endl;
      for(list< list<QMCWalker>::iterator >::iterator wp=del.begin();
	  wp!=del.end(); ++wp)
	cout << "w= " << (*wp)->getWeight() << " sw = " << fabs((*wp)->getWeight()/aveW-1) << endl;

      cout << "aCandidates(" << aCandidates.size() << "):" << endl;
      for(map<double,list<QMCWalker>::iterator>::iterator it=aCandidates.begin();
	  it != aCandidates.end();
	  it++)
	printf("idx %10.5f w = %20.10f sw= %20.10f\n",(*it).first,(*it).second->getWeight(),fabs((*it).second->getWeight()/aveW-1));
      cout << "dCandidates(" << dCandidates.size() << "):" << endl;
      for(map<double,list<QMCWalker>::iterator>::iterator it=dCandidates.begin();
	  it != dCandidates.end();
	  it++)
	printf("idx %10.5f w = %20.10f sw= %20.10f\n",(*it).first,(*it).second->getWeight(),fabs((*it).second->getWeight()/aveW-1));

      /*
      cout << "New Walker List:" << endl;
      for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end();++wp) 
	cout << wp->ID(true,false);
      */
      //exit(0);
    }  

  for(list< list<QMCWalker>::iterator >::iterator wtd=del.begin();
      wtd!=del.end(); ++wtd)
    wlist.erase( *wtd );
  del.clear();

  wlist.splice( wlist.end(), add );
  
  for(list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end();++wp)
    wp->setWeight(aveW);

}

double QMCRun::getWeights()
{
  // determine the total of all the walker weights
  double total_weights = 0.0;
  
  for( list<QMCWalker>::iterator wp=wlist.begin(); wp!=wlist.end(); ++wp )
    total_weights += wp->getWeight();

  return total_weights;
}

double QMCRun::getPopulationSizeBiasCorrectionFactor()
{
  return populationSizeBiasCorrectionFactor;
}

void QMCRun::toXML(ostream& strm)
{
  strm << "<runType>\n\t" << Input->flags.run_type << endl;
  strm << "</runType>\n";

  if (globalInput.flags.run_type == "diffusion")
    {
      strm << "<populationSizeBiasCorrectionFactor>\n\t"
	   << populationSizeBiasCorrectionFactor
	   << "\n</populationSizeBiasCorrectionFactor>" << endl;
  
      queue<double> temp_CD = correctionDivisor;

      int cDsize = temp_CD.size();

      strm << "<correctionDivisorSize>\n\t" << cDsize
	   << "\n</correctionDivisorSize>" << endl;
      strm << "<correctionDivisor>\n";
      for (int i=0; i<cDsize; i++)
	{
	  strm << temp_CD.front() << endl;
	  temp_CD.pop();
	}
      strm << "</correctionDivisor>" << endl;
    }
  // writes out the properties
  if (Input->flags.use_equilibration_array == 1)
    EquilibrationArray.toXML(strm);
  else 
    {
      Properties.toXML(strm);
      if (Input->flags.future_walking.size() > 1)
	fwProperties.toXML(strm);
    }

  //prints all the walkers
  strm << "<walkers>" << endl;
  
  for(list <QMCWalker>::iterator wp=wlist.begin();wp!=wlist.end();++wp)
    {
      wp->toXML(strm);
    }    
  strm << "</walkers>\n" << endl;
}

bool QMCRun::readXML(istream& strm)
{
  string temp;

  strm >> temp;
  if (temp != "<runType>")
    return false;
  strm >> temp;
  string runType = temp;
  strm >> temp;
  if (temp != "</runType>")
    return false;

  if (runType == "diffusion")
    {
      // read populationSizeBiasCorrectionFactor
      strm >> temp;
      if (temp != "<populationSizeBiasCorrectionFactor>")
	return false;
      strm >> temp;
      populationSizeBiasCorrectionFactor = atof(temp.c_str());
      strm >> temp;
      if (temp != "</populationSizeBiasCorrectionFactor>")
	return false;
  
      strm >> temp;
      if (temp != "<correctionDivisorSize>")
	return false;
      strm >> temp;
      int cDsize = atoi(temp.c_str());
      strm >> temp;
      if (temp != "</correctionDivisorSize>")
	return false;

      strm >> temp;
      if (temp != "<correctionDivisor>")
	return false;
      for (int i=0; i<cDsize; i++)
	{
	  strm >> temp;
	  correctionDivisor.push(atof(temp.c_str()));
	}
      strm >> temp;
      if (temp != "</correctionDivisor>")
	return false;
    }

  // read the properties
  if (Input->flags.use_equilibration_array == 1)
    {
      if (!EquilibrationArray.readXML(strm))
	return false;
    } 
  else 
    {
      if (!Properties.readXML(strm))
	return false;
      if (Input->flags.future_walking.size() > 1)
	if (!fwProperties.readXML(strm))
	  return false;
    }
  
  // read the walkers
  
  strm >> temp;
  if (temp != "<walkers>")
    return false;
  wlist.clear();
  for(int i=0; i<Input->flags.number_of_walkers; i++)
    {
      QMCWalker w;
      w.initialize(Input);
      
      // read a walker
      if (!w.readXML(strm,*QMF))
        return false;
      
      wlist.push_back(w);
    }
    
  strm >> temp;
  if (temp != "</walkers>")
    return false;

  return true;
}

void QMCRun::startTimers()
{
  EquilibrationArray.startTimers();
}

void QMCRun::stopTimers()
{
  EquilibrationArray.stopTimers();
}

Stopwatch * QMCRun::getPropagationStopwatch()
{
  return EquilibrationArray.getPropagationStopwatch();
}

Stopwatch * QMCRun::getEquilibrationStopwatch()
{
  return EquilibrationArray.getEquilibrationStopwatch();
}

QMCProperties * QMCRun::getTimeStepProperties()
{
  return &timeStepProperties;
}

QMCProperties * QMCRun::getProperties()
{
  if (Input->flags.use_equilibration_array == 1)
    return EquilibrationArray.chooseDecorrObject();
  else
    return &Properties;
}

QMCPropertyArrays * QMCRun::getFWTimeStepProperties()
{
  return &fwTimeStepProperties;
}

QMCPropertyArrays * QMCRun::getFWProperties()
{
  /*
    AGA: FIX THIS
  if (Input->flags.use_equilibration_array == 1)
    {
      return EquilibrationArray.chooseDecorrObject();
    } else {
      return &fwProperties;
    }
  */
  return &fwProperties;
}

int QMCRun::getNumberOfWalkers()
{
  return (int)wlist.size();
}

bool QMCRun::step(bool writeConfigs, int iteration)
{
  swTimers(3).start();

  propagateWalkers(writeConfigs,iteration);
  calculatePopulationSizeBiasCorrectionFactor();

  swTimers(2).start();
  calculateObservables();
  swTimers(2).lap(getNumberOfWalkers());

  int step = iteration + Input->flags.equilibration_steps - 1;

  int whichE = -1;
  if(globalInput.flags.one_e_per_iter != 0)
    whichE = step % globalInput.WF.getNumberElectrons();

  if(whichE == -1 || whichE == 0)
    {
      /*
	In the HLR QMC book pg 273, they recommend branching only once
	per N electron cycle.

	Why? Is this important?
      */
      growthRate = getNumberOfWalkers();
      branchWalkers();
      growthRate -= getNumberOfWalkers();
    }

  swTimers(3).lap(getNumberOfWalkers());

  if(getWeights() <= 0.0 || getNumberOfWalkers() == 0)
    {
      cerr << "Error on node " << Input->flags.my_rank << ":" << endl;
      cerr << "       number walkers   = " << getNumberOfWalkers() << endl;
      cerr << "       total weight     = " << getWeights() << endl;
      cerr << "       energy_estimated = " << Input->flags.energy_estimated << endl; 
      cerr << "       energy_trial     = " << Input->flags.energy_trial << endl; 
      cerr << "       energy_original  = " << Input->flags.energy_estimated_original << endl; 

      if(globalInput.flags.nprocs > 4)
	{
	  //Assume that the calculation can be salvaged; only this node was bad.
	  cerr << "       we're going to reinitialize all walkers..." << endl;
	} else {
	  //You're screwed.
	  cerr << "Your calculation is junk! Try lowering dt." << endl;
	  exit(1);
	}

      Input->flags.energy_estimated = Input->flags.energy_estimated_original;
      Input->flags.energy_trial     = Input->flags.energy_estimated_original;
      randomlyInitializeWalkers();
      return false;
    }
  return true;
}

void QMCRun::calculatePopulationSizeBiasCorrectionFactor()
{
  if( Input->flags.correct_population_size_bias &&
      Input->flags.run_type != "variational" )
    {
      /*
	Discussion on Tp:
	Tp needs to decrease with higher population size
	Tp needs to increase with with increasing rms fluctuation of E_local
	Tp needs to be proportional to the autocorrelation time
	UNR93 give no further recommendation except to say that while
	the distribution is unbiased only for infinite Tp, the variance
	goes up with increasing Tp.

	Preliminary tests seem to indicate that 1/dt is best
	(as compared to 0.25/dt, 0.5/dt, 2/dt, 4/dt, and 8/dt)
      */
      const unsigned int Tp = (int)(1.0/Input->flags.dt);

      /*
	if we're going to "absorb the constant part into the definition of G"
	then don't we have to add Input->flags.energy_estimated_original back in somewhere?
	then again, perhaps if walkers at all times are uniformly weighted by this factor, then
	after equilibration, it doesn't really matter... Et - Eo seems to work better.
      */
      double temp = Input->flags.energy_trial - Input->flags.energy_estimated_original;
      //double temp = Input->flags.energy_trial;

      //this used to be just dt, but UNR93 says dt_effective
      temp *= -Input->flags.dt_effective;
      
      temp = exp(temp);
      populationSizeBiasCorrectionFactor *= temp;

      //this is not supposed to be cleared when transitioning from 
      //equilibration to production steps.
      correctionDivisor.push(temp);
      if(correctionDivisor.size() > Tp)
	{
	  populationSizeBiasCorrectionFactor /= correctionDivisor.front();
	  correctionDivisor.pop();
	}

      if (IeeeMath::isNaN(populationSizeBiasCorrectionFactor))
	{
	  cerr << "Error in QMCRun::calculatePopulationSizeBiasCorrectionFactor()" << endl;
	  cerr << "  energy_trial = " << Input->flags.energy_trial << endl;
	  cerr << "  dt_effective = " << Input->flags.dt_effective << endl;
	  cerr << "  multiplier   = " << temp << endl;
	  exit(1);
	}
    }
}

// Updates hartree-fock potential, adding all electrons over all walkers to the
// object.
void QMCRun::updateHFPotential()
{
  int numelecs = Input->WF.getNumberElectrons();

  for (list<QMCWalker>::iterator wp = wlist.begin(); wp != wlist.end(); ++wp)
    {
      Array2D<double> positions = *(wp->getR());
      for (int i = 0; i < numelecs; i++)
	HartreeFock.AddElectron(i,wp->getWeight(),positions(i,0),
				positions(i,1), positions(i, 2));
      HartreeFock.IncrementSample();
    }
}

int QMCRun::getNumTimers()
{
  return QMF->getNumTimers() + 4;
}

void QMCRun::aggregateTimers(Array1D<Stopwatch> & timers,
			     int & idx)
{
  QMF->aggregateTimers(timers,idx);
  for(int i=0; i<swTimers.dim1(); i++)
    timers(idx++).aggregateTimer(swTimers(i));
}
