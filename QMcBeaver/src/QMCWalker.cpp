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

#include "QMCWalker.h"
#include "QMCSurfer.h"
#include "QMCPotential_Energy.h"
#include "math.h"

const double QMCWalker::maxFWAsymp = 1.0;
long int QMCWalker::nextID = 0;
//const int QMCWalker::numAncestors = 10;

QMCWalker::QMCWalker()
{
  TrialWalker    = 0;
  OriginalWalker = 0;

  dW          =  1.0;
  weight      =  1.0;
  age         =  0;
  dR2         =  0.0;
  ageMoved    = -1;
  numWarnings =  0;
  iteration   =  0;
}

void QMCWalker::branchID()
{
  for(int i=numAncestors-1; i > 0; i--)
    {
      genealogy[i]   = genealogy[i-1];
      birthIter[i]   = birthIter[i-1];
      birthWeight[i] = birthWeight[i-1];
      numMoves[i]    = numMoves[i-1];
    }
  birthIter[0]   = iteration;
  birthWeight[0] = weight;
  genealogy[0]   = ++nextID;
  numMoves[0]    = 0;
}

void QMCWalker::newID()
{
  birthIter[0]   = iteration;
  birthWeight[0] = weight;
  genealogy[0]   = ++nextID;
  numMoves[0]    = 0;
  for(int i=1; i < numAncestors; i++)
    {
      genealogy[i]   = -1;
      birthIter[i]   = -1;
      birthWeight[i] = -1;
      numMoves[i]    = -1;
    }
}

QMCWalker::QMCWalker( const QMCWalker & rhs )
{
  TrialWalker    = 0;
  OriginalWalker = 0;

  *this = rhs;
}

QMCWalker::~QMCWalker()
{
  TrialWalker = 0;

  delete OriginalWalker;
  OriginalWalker = 0;  
  
  Input = 0;
  
  R.deallocate();
  
  numFWSteps.clear();

  isCollectingFWResults.deallocate();  
  fwNormalization.deallocate();
  fwR12.deallocate();
  fwR2.deallocate();
  fwiR12.deallocate();
  fwiR.deallocate();
  fwEnergy.deallocate();
  fwKineticEnergy.deallocate();
  fwKineticEnergy_grad.deallocate();
  fwPotentialEnergy.deallocate();

  cs_Energies.deallocate();
  cs_Weights.deallocate();
  p3_xxa.deallocate();
  rp_a.deallocate();

  for(int d1=0; d1<fwNuclearForces.dim1(); d1++)
    for(int d2=0; d2<fwNuclearForces.dim2(); d2++)
      fwNuclearForces(d1,d2).deallocate();
  fwNuclearForces.deallocate();  
}

void QMCWalker::operator=( const QMCWalker & rhs )
{
  // In the equality operator DON'T copy the pointers to the child walkers.
  // If you do, then there is the possibility of big memory leak problems.
  // This is not required because these are just temporary variables used in
  // propagating the electrons.
  
  weight      = rhs.weight;
  age         = rhs.age;
  ageMoved    = rhs.ageMoved;
  numWarnings = rhs.numWarnings;

  for(int i=0; i<numAncestors; i++)
    {
      birthIter[i]   = rhs.birthIter[i];
      birthWeight[i] = rhs.birthWeight[i];
      genealogy[i]   = rhs.genealogy[i];
      numMoves[i]    = rhs.numMoves[i];
    }

  Input = rhs.Input;
  move_accepted         = rhs.move_accepted;
  AcceptanceProbability = rhs.AcceptanceProbability;
  localEnergy           = rhs.localEnergy;
  kineticEnergy         = rhs.kineticEnergy;
  potentialEnergy       = rhs.potentialEnergy;
  neEnergy              = rhs.neEnergy;
  eeEnergy              = rhs.eeEnergy;
  
  fwNormalization       = rhs.fwNormalization;
  fwR12                 = rhs.fwR12;
  fwR2                  = rhs.fwR2;
  fwiR12                = rhs.fwiR12;
  fwiR                  = rhs.fwiR;
  fwEnergy              = rhs.fwEnergy;
  fwKineticEnergy       = rhs.fwKineticEnergy;
  fwKineticEnergy_grad  = rhs.fwKineticEnergy_grad;
  fwPotentialEnergy     = rhs.fwPotentialEnergy;
  isCollectingFWResults = rhs.isCollectingFWResults;
  numFWSteps            = rhs.numFWSteps;
  
  if(Input->flags.nuclear_derivatives != "none"){
    fwNuclearForces.allocate(rhs.fwNuclearForces.dim1(),
                             rhs.fwNuclearForces.dim2());
    
    for (int d1=0; d1<fwNuclearForces.dim1(); d1++)
    {
      for (int d2=0; d2<fwNuclearForces.dim2(); d2++)
      {
        fwNuclearForces(d1,d2).allocate(rhs.fwNuclearForces.get(d1,d2).dim1(),2);
        fwNuclearForces(d1,d2) = rhs.fwNuclearForces.get(d1,d2);
      }
    }
  }

  cs_Energies.allocate(rhs.cs_Energies.dim1());
  cs_Weights.allocate(rhs.cs_Weights.dim1());
  p3_xxa.allocate(rhs.p3_xxa.dim1());
  rp_a.allocate(rhs.rp_a.dim1());
       
  distanceMovedAccepted = rhs.distanceMovedAccepted;
  dR2                   = rhs.dR2;
  R                     = rhs.R;
  walkerData            = rhs.walkerData;
}

void QMCWalker::initializePropagation(QMCWalkerData * &data,
                                      Array2D<double> * &rToCalc,
				      int iteration)
{
  this->iteration = iteration;
  createChildWalkers();

  int step = iteration + Input->flags.equilibration_steps - 1;

  int whichE = -1;
  if(globalInput.flags.one_e_per_iter != 0)
    whichE = step % R.dim1();

  TrialWalker->walkerData.whichE = whichE;
  OriginalWalker->walkerData.whichE = whichE;

  forwardGreensFunction = TrialWalker->moveElectrons();
  data = & TrialWalker->walkerData;
  rToCalc = & TrialWalker->R;
}

void QMCWalker::processPropagation(QMCFunctions & QMF, bool writeConfigs)
{
  if(globalInput.flags.use_surfer == 1)
    {
      static int iteration_to_stop = iteration;

      if(iteration == iteration_to_stop)
	{
	  static QMCSurfer s;
	  cout << "\n**** Surfing walker " << ID(false,false);
	  iteration_to_stop = s.mainMenu(&QMF, iteration, R);
	}
    }


  QMCDouble reverseGreensFunction =
    calculateReverseGreensFunction();
  double GreensFunctionRatio =
    (double)(reverseGreensFunction/forwardGreensFunction);

  if (IeeeMath::isNaN(GreensFunctionRatio)){
    cout << "Error: bad Green's ratio " << GreensFunctionRatio << endl;
    cout << " Forward = " << forwardGreensFunction << endl;
    cout << " Reverse = " << reverseGreensFunction << endl;
    calculateMoveAcceptanceProbability(0.0);
  }
  else
    calculateMoveAcceptanceProbability(GreensFunctionRatio);
  
  acceptOrRejectMove();

  //exit(0);

  if(walkerData.whichE == -1 ||
     walkerData.whichE == 0)
    reweight_walker();

  if(getWeight() == 0)
    return;

  if(move_accepted == false)
    {
      dR2                   = OriginalWalker->dR2;
    }

  if(writeConfigs)
    {
      /*
	If we want to be able to exactly recompose the statistics
	off of a cfgs file, then we'll need to print both configurations,
	since both are used in calculateObservables().
      */
      double p = TrialWalker->AcceptanceProbability;
      double q = 1.0 - p;

      TrialWalker->walkerData.writeConfigs(TrialWalker->R,p);
      OriginalWalker->walkerData.writeConfigs(OriginalWalker->R,q);
    }

  calculateObservables();
  
  if(move_accepted == false)
    {
      R                     = OriginalWalker->R;
      AcceptanceProbability = OriginalWalker->AcceptanceProbability;

      if(globalInput.flags.one_e_per_iter == 0)
	{
	  walkerData            = OriginalWalker->walkerData;
	} else {	 
	  walkerData.updateDistances(R);
	  QMF.evaluate(R,walkerData);
	}
    }
  
  if( TrialWalker->isSingular() )
    {
      cerr << "WARNING: Reinitializing a singular walker!! " << ID(true,false) << endl;
      initializeWalkerPosition(QMF);
    }

  if(writeConfigs)
    {
      /*
	Assuming the AcceptanceProbability is near one, this
	is probably more efficient than the option above.
       */
      //walkerData.writeConfigs(R, AcceptanceProbability);
    }
}

QMCDouble QMCWalker::moveElectrons()
{
  // Make sure that dt is positive!
  if( Input->flags.dt <= 0 )
    {
      cerr << "ERROR: Negative dt value! (" << Input->flags.dt << ")" << endl;
      exit(0);
    }
  QMCDouble gr = QMCDouble();
  
  if(Input->flags.sampling_method == "no_importance_sampling")
    gr = moveElectronsNoImportanceSampling();
  else if(Input->flags.sampling_method == "importance_sampling" )
    gr = moveElectronsImportanceSampling();
  else if(Input->flags.sampling_method == "umrigar93_importance_sampling")
    gr = moveElectronsUmrigar93ImportanceSampling();
  else if(Input->flags.sampling_method == "umrigar93_accelerated_sampling")
    gr = moveElectronsUmrigar93AcceleratedSampling();
  else
    {
      cerr << "ERROR: Improper value for sampling_method set ("
	   << Input->flags.sampling_method << ")!" << endl;
      exit(0);
    }

  walkerData.updateDistances(R);
  return gr;
}

QMCDouble QMCWalker::moveElectronsNoImportanceSampling()
{
  int whichE = walkerData.whichE;

  // Move the electrons of this walker
  double sigma = sqrt(Input->flags.dt);
  dR2 = 0;  
  
  // Don't use the QF  =>  R' = R + gauss_rn
  for(int i=0; i<R.dim1(); i++)
    {
      if(whichE != -1 && i != whichE)
	continue;

      for(int j=0; j<R.dim2(); j++)
	{
	  // Add the randomness to the displacement
	  double drift = sigma*ran.gasdev();

	  // Calculate the square of the magnitude of the displacement
	  dR2 += drift * drift;

	  // Now update the R
	  R(i,j) += drift;
	}
    }
      
  return QMCDouble(1.0);
}

QMCDouble QMCWalker::moveElectronsImportanceSampling()
{
  Array2D<double> Displacement(R.dim1(),R.dim2());
  double tau    = Input->flags.dt;
  double sigma  = sqrt(tau);
  double greens = 0.0;
  bool toMove;
  
  // Use the QF => R' = R + dt * QF + gauss_rn
  // We have already counted for the factor of 2.0
  Displacement  = walkerData.modifiedGradPsiRatio;
  Displacement *= tau;
  dR2 = 0.0;  
  for(int i=0; i<R.dim1(); i++)
    {
      toMove = true;
      if(walkerData.whichE != -1 && i != walkerData.whichE)
	toMove = false;

      if(toMove)
	{
	  for(int j=0; j<R.dim2(); j++)
	    {
	      // Add the randomness to the displacement
	      double drift = sigma*ran.gasdev();
	      
	      // Add the randomness to the displacement
	      Displacement(i,j) += drift;
	      
	      // Calculate the square of the magnitude of the displacement
	      dR2 += Displacement(i,j) * Displacement(i,j);
	      
	      // Now update the R
	      R(i,j) += Displacement(i,j);
	      
	      // Calculate the Green's function for the forward move
	      // The 'Quantum Force' term cancels out
	      greens += drift*drift;
	    }
	} 
      else
	{
	  for(int j=0; j<R.dim2(); j++)
	    //If we're not moving the electron, then the quantum force doesn't cancel
	    greens += Displacement(i,j) * Displacement(i,j);      
	}
    }
    
  //k a^2 exp(c)
  double k = 1.0;
  double a = 2.0*pi*tau;
  double b = -0.5*R.dim1()*R.dim2();
  double c = -greens/(2.0*tau);
  
  return QMCDouble(k,a,b,c);
}

QMCDouble QMCWalker::moveElectronsUmrigar93ImportanceSampling()
{
  int whichE = walkerData.whichE;
  double tau = Input->flags.dt;
  bool toMove;

  dR2 = 0.0;
  Array2D<double> & Displacement = walkerData.modifiedGradPsiRatio;
  
  Array1D<double> zUnitVector(3);
  Array1D<double> radialUnitVector(3);
  Array1D<double> newPosition(3);
  
  QMCDouble GF(1.0);
  QMCDouble GaussianGF;
  QMCDouble SlaterGF;
  QMCDouble OneE;

  for(int electron=0; electron<Input->WF.getNumberElectrons(); electron++)
    {
      toMove = true;
      if(whichE != -1 && electron != whichE)
	toMove = false;

      // Find the nearest nucleus to this electron
      int nearestNucleus = Input->Molecule.findClosestNucleusIndex(walkerData.riA,electron);
      double distanceFromNucleus = walkerData.riA(electron,nearestNucleus);
      
      for(int i=0; i<3; i++)
	zUnitVector(i) = walkerData.riA_uvec(electron,nearestNucleus,i);
      
      // Decompose the modified QF into components in the z and radial
      // directions
      
      double zComponentQF = 0.0;
      
      for(int i=0; i<3; i++)
	zComponentQF += zUnitVector(i)*Displacement(electron,i);
        
      double radialComponentQF = 0.0;
      
      for(int i=0; i<3; i++)
        {
          radialUnitVector(i) = Displacement(electron,i) -
                                zComponentQF * zUnitVector(i);
                                
          radialComponentQF += radialUnitVector(i)*radialUnitVector(i);
        }
        
      radialComponentQF = sqrt(radialComponentQF);
      
      if (radialComponentQF > 1e-15)
	radialUnitVector /= radialComponentQF;
      
      // Calculate the gaussian drift components in cylindrical coordinates
      double zCoordinate = max(distanceFromNucleus + zComponentQF * tau,0.0);
                               
      double radialCoordinate = 2.0 * radialComponentQF * tau *
                           zCoordinate / ( distanceFromNucleus + zCoordinate );
                                
      // Calculate the hydrogen like exponential factor
      
      double expParam = Input->Molecule.Z(nearestNucleus);
      expParam = expParam*expParam + 1.0/tau;
      expParam = sqrt(expParam);
      
      // Calculate the probability of moving the electron with respect
      // to a gaussian type distribution with a QF drift or a hydrogenic
      // atom type slater function centered on the closest nucleus
      
      /*
	double probabilitySlaterTypeMove = 0.5 *
	erfc( (distanceFromNucleus + zComponentQF * tau) / sqrt(2.0 * tau) );
	/*/
      double probabilitySlaterTypeMove = 0.5 *
      MathFunctions::erfc( (distanceFromNucleus + zComponentQF * tau) / sqrt(2.0 * tau) );
      //*/

      if (probabilitySlaterTypeMove > 1.0)
	{
	  cerr << "Warning: probabilitySlaterTypeMove = ";
	  cerr << probabilitySlaterTypeMove << endl;
	  probabilitySlaterTypeMove = 1.0;
	}
      if (probabilitySlaterTypeMove < 0.0)
	{
	  cerr << "Warning: probabilitySlaterTypeMove = ";
	  cerr << probabilitySlaterTypeMove << endl;
	  probabilitySlaterTypeMove = 0.0;
	}
                                                              
      double probabilityGaussianTypeMove = 1.0 - probabilitySlaterTypeMove;
      
      if(toMove)
	{
	  // Randomly decide which electron moving method to use
	  if( probabilityGaussianTypeMove > ran.unidev() )
	    {
	      // Gaussian Type Move
	      
	      // Particle is moved in the direction
	      // Rnuc + radialCoordinate * radialUnitVector +
	      //   zCoordinate * zUnitVector +
	      //   gaussian random number with standard deviation sqrt(tau)
	      for(int i=0; i<3; i++)
		newPosition(i) = Input->Molecule.Atom_Positions(nearestNucleus,i)
		  + radialCoordinate * radialUnitVector(i)
		  + zCoordinate * zUnitVector(i)
		  + sqrt(tau)*ran.gasdev();
	    }
	  else
	    {
	      // Slater Type Move
	      
	      for(int i=0; i<3; i++)
		newPosition(i) = Input->Molecule.Atom_Positions(nearestNucleus,i);
	      
	      // add random part
	      
	      double r = ran.randomDistribution1()/(2.0*expParam);
	      double phi = 2*pi*ran.unidev();
	      double theta = ran.sindev();
	      
	      newPosition(0) += r*sin(theta)*cos(phi);
	      newPosition(1) += r*sin(theta)*sin(phi);
	      newPosition(2) += r*cos(theta);
	    }
	} else {
	  //We aren't moving this electron
	  //So set it to be the original position
	  for(int i=0; i<3; i++)
	    newPosition(i) = R(electron,i);
	}

      double temp;
      // Update the distance attempted to move squared and move the electron.
      for(int i=0; i<3; i++)
        {
          temp = newPosition(i) - R(electron,i);
          dR2 += temp*temp;
          R(electron,i) = newPosition(i);
        }

      // Update the greens function
      double distance1Sq = 0.0;
      
      for(int i=0; i<3; i++)
        {
          temp = newPosition(i) -
	    ( Input->Molecule.Atom_Positions(nearestNucleus,i) +
	      radialCoordinate * radialUnitVector(i) + zCoordinate * zUnitVector(i) );
	  
          distance1Sq += temp*temp;
        }
      
      double ga = 2*pi*tau;
      double gb = -1.5;
      double gc = -distance1Sq/(2*tau);
      
      GaussianGF=QMCDouble(probabilityGaussianTypeMove,ga,gb,gc);

      /*
      if(probabilitySlaterTypeMove < 1e-20)
	{
	  GF *= GaussianGF;
	  continue;	  
	}
      */
      double distance2Sq = 0.0;
      
      for(int i=0; i<3; i++)
        {
          temp = newPosition(i) -
	    Input->Molecule.Atom_Positions(nearestNucleus,i);
	  
          distance2Sq += temp*temp;
        }
      
      double sk = probabilitySlaterTypeMove/pi;
      double sa = expParam;
      double sb = 3.0;
      double sc = -2.0*expParam*sqrt(distance2Sq);
      
      SlaterGF = QMCDouble(sk,sa,sb,sc);
      
      // The addition and multiplication operations here have been causing a
      // lot of problems.
      OneE = SlaterGF + GaussianGF;
      
      GF *= OneE;
    }
  return GF;
}

QMCDouble QMCWalker::moveElectronsUmrigar93AcceleratedSampling()
{
  int whichE = walkerData.whichE;
  double cos_tm = cos(globalInput.flags.accel_tm);
  double delta = globalInput.flags.accel_delta;
  bool toMove;

  dR2 = 0.0;
  Array2D<double> & Fk = walkerData.modifiedGradPsiRatio;
  //Array2D<double> & Fk = walkerData.gradPsiRatio;
  
  Array1D<double> riA_hat(3);
  Array1D<double> Fx(3);
  Array1D<double> newPosition(3);
  
  QMCDouble GF(1.0);
  for(int electron=0; electron<Input->WF.getNumberElectrons(); electron++)
    {
      toMove = true;
      if(whichE != -1 && electron != whichE)
	toMove = false;
      if(!toMove) continue;

      // Find the nearest nucleus to this electron
      int nearestNucleus = Input->Molecule.findClosestNucleusIndex(walkerData.riA,electron);
      double Z = Input->Molecule.Z(nearestNucleus);
      double riA = walkerData.riA(electron,nearestNucleus);
      
      for(int xyz=0; xyz<3; xyz++)
	riA_hat(xyz) = walkerData.riA_uvec(electron,nearestNucleus,xyz);      

      /*
	First, calculate the new radius. To do this, we need to fit
	parameters a and zeta so that the U function matches log derivatives.
	Frk = rki . Fk
      */
      double Frk = 0.0;      
      for(int xyz=0; xyz<3; xyz++)
	Frk += riA_hat(xyz)*Fk(electron,xyz);

      //Get zeta and a
      double zeta, a, I;
      MathFunctions::fitU(Z,Frk,riA,delta,zeta,a,I);

      //Now that we have a and zeta, we can sample the function sqrt(rf) |U(rf)|
      double rkf;
      if(toMove) rkf = ran.randomDistribution2(a,zeta,riA/delta,riA*delta);
      else	 rkf = riA;
      double rav = (riA + rkf)/2.0;

      /*
	Obtain the new theta
      */
      double tM, cos_tM;
      cos_tM = cos_tm - (1.0 + cos_tm)/(1.0 + Z*Z*rav*rav);
      tM = acos(cos_tM);

      double tkf;
      if(toMove){
	tkf = ran.sindev(tM);	
      } else {
	tkf = 0.0;
      }
      double s_tkf = sin(tkf);
      double c_tkf = cos(tkf);

      /*
	Obtain the new phi
	Fxk = | Fk - Frk rk |
      */
      double Fxk = 0.0;
      for(int xyz=0; xyz<3; xyz++)
        {
          Fx(xyz) = Fk(electron,xyz) - Frk*riA_hat(xyz);
          Fxk += Fx(xyz)*Fx(xyz);
        }
      Fxk = sqrt(Fxk);

      double Fpk = Fxk*rav*s_tkf;
      double pkf;
      if(toMove) pkf = ran.randomDistribution3(Fpk);
      else       pkf = 0.0;

      // Our new position is rotated theta degrees in the direction of Fx
      for(int xyz=0; xyz<3; xyz++)
	newPosition(xyz) = c_tkf*riA_hat(xyz) + s_tkf*Fx(xyz)/Fxk;

      // Rotate our new position phi degrees around riA_hat
      newPosition.rotate(riA_hat,pkf);
      newPosition *= rkf;

      double temp;
      for(int xyz=0; xyz<3; xyz++)
        {
	  //Move the electron
          temp = newPosition(xyz) - R(electron,xyz);
          R(electron,xyz) = newPosition(xyz);
          dR2 += temp*temp;
	}

      /*
	Finally, calculate the transition matrix terms
	S(Rf|Ri) is the current going to Rf from Ri
	S(Rf|Ri) = | Phi(Rf|Ri) | / rkf^1.5 (eqn 24)
	Phi(Rf|Ri) = W V U
      */
      double S = 1.0;
      S *= (1.0 + a*rkf)*exp(-zeta*rkf);
      S *= 1.0 + min(Fpk,1.0)*cos(pkf);
      S *= pow(rkf,-1.5);
      S = fabs(S);

      I *= (1.0 - cos_tM);
      I *= 2*pi;

      //T(Rf|Ri) = S(Rf|Ri)/I(Ri) is probability of going from Ri to Rf 
      if(fabs(I) > 1e-50)
	GF *= QMCDouble(S/I);

      if(GF.isNotValid()){
	printf("move electron=%i age=%i\n",electron,age);
	printf("rki=%20.10f a=%20.10f zeta=%20.10f i=%20.10f o=%20.10f tM=%20.10e\n",riA,a,zeta,riA/delta,riA*delta,tM);
	printf("Frk=%20.10f Fxk=%20.10f Fpk=%20.10f\n",Frk,Fxk,Fpk);
	printf("rkf=%20.10f tkf=%20.10f pkf=%20.10f\n",rkf,tkf,pkf); 
	printf("S=%20.10e I=%20.10e T=%20.10e\n",S,I,(double)(GF));
      }
    }
  return GF;
}

QMCDouble QMCWalker::calculateReverseGreensFunction()
{
  // Make sure that dt is positive!
  if( Input->flags.dt <= 0 )
    {
      cerr << "ERROR: Negative dt value! (" << Input->flags.dt << ")" << endl;
      exit(0);
    }
    
  if(Input->flags.sampling_method == "no_importance_sampling")
      return calculateReverseGreensFunctionNoImportanceSampling();
  else if(Input->flags.sampling_method == "importance_sampling" )
      return calculateReverseGreensFunctionImportanceSampling();
  else if(Input->flags.sampling_method == "umrigar93_importance_sampling")
      return calculateReverseGreensFunctionUmrigar93ImportanceSampling();
  else if(Input->flags.sampling_method == "umrigar93_accelerated_sampling")
      return calculateReverseGreensFunctionUmrigar93AcceleratedSampling();
  else
    {
      cerr << "ERROR: Improper value for sampling_method set ("
	   << Input->flags.sampling_method << ")!" << endl;
      exit(0);
    }
    
  // should never reach this
  
  QMCDouble TRASH = QMCDouble();
  return TRASH;
}

QMCDouble \
QMCWalker::calculateReverseGreensFunctionNoImportanceSampling()
{
  return QMCDouble(1.0);
}

QMCDouble \
QMCWalker::calculateReverseGreensFunctionImportanceSampling()
{
  double tau = Input->flags.dt;  
  double greens = 0.0;

  for(int i=0; i<R.dim1(); i++)
    {      
      for(int j=0; j<R.dim2(); j++)
	{
	  double temp = OriginalWalker->R(i,j) - TrialWalker->R(i,j) -
	    tau*TrialWalker->walkerData.modifiedGradPsiRatio(i,j);
	  
	  greens += temp*temp;
	}
    }

  // k * a^b * exp(c)
  double k = 1.0;
  double a = 2.0*pi*tau;
  double b = -0.5*R.dim1()*R.dim2();
  double c = -greens/(2.0*tau);

  return QMCDouble(k,a,b,c);
}

QMCDouble \
QMCWalker::calculateReverseGreensFunctionUmrigar93ImportanceSampling()
{
  int whichE = walkerData.whichE;
  bool moved;

  double tau = Input->flags.dt;
  Array1D<double> zUnitVector(3);
  Array1D<double> radialUnitVector(3);

  QMCDouble GF(1.0);
  QMCDouble GaussianGF;
  QMCDouble SlaterGF;
  QMCDouble OneE;

  for(int electron=0; electron<Input->WF.getNumberElectrons(); electron++)
    {      
      moved = true;
      if(whichE != -1 && electron != whichE)
	moved = false;

      // Find the nearest nucleus to this electron
      int nearestNucleus = Input->Molecule.findClosestNucleusIndex(TrialWalker->walkerData.riA,electron);
      double distanceFromNucleus = walkerData.riA(electron,nearestNucleus);
      
      for(int i=0; i<3; i++)
	zUnitVector(i) = walkerData.riA_uvec(electron,nearestNucleus,i);
      
      // Decompose the modified QF into components in the z and radial
      // directions
      
      double zComponentQF = 0.0;
      
      for(int i=0; i<3; i++)
          zComponentQF += zUnitVector(i)*
                          TrialWalker->walkerData.modifiedGradPsiRatio(electron,i);
        
      double radialComponentQF = 0.0;
      
      for(int i=0; i<3; i++)
        {
          radialUnitVector(i) =
            TrialWalker->walkerData.modifiedGradPsiRatio(electron,i) -
            zComponentQF * zUnitVector(i);
            
          radialComponentQF += radialUnitVector(i)*radialUnitVector(i);
        }
        
      radialComponentQF = sqrt(radialComponentQF);
      
      if (radialComponentQF > 1e-15)
	radialUnitVector /= radialComponentQF;
      
      // Calculate the gaussian drift components in cylindrical coordinates
      double zCoordinate = max(distanceFromNucleus + zComponentQF * tau,0.0);
                               
      double radialCoordinate = 2.0 * radialComponentQF * tau *
                           zCoordinate / ( distanceFromNucleus + zCoordinate );
                                
      // Calculate the hydrogen like exponential factor
      
      double expParam = Input->Molecule.Z(nearestNucleus);
      expParam = expParam*expParam + 1.0/tau;
      expParam = sqrt(expParam);
      
      // Calculate the probability of moving the electron with respect
      // to a gaussian type distribution with a QF drift or a hydrogenic
      // atom type slater function centered on the closest nucleus
      
      /*
      double probabilitySlaterTypeMove = 0.5 *
	erfc( (distanceFromNucleus + zComponentQF * tau) / sqrt(2.0 * tau) );
	/*/
      double probabilitySlaterTypeMove = 0.5 *
      MathFunctions::erfc( (distanceFromNucleus + zComponentQF * tau) / sqrt(2.0 * tau) );
      //*/

      if (probabilitySlaterTypeMove > 1.0)
	{
	  cerr << "Warning: probabilitySlaterTypeMove = ";
	  cerr << probabilitySlaterTypeMove << endl;
	  probabilitySlaterTypeMove = 1.0;
	}
      if (probabilitySlaterTypeMove < 0.0)
	{
	  cerr << "Warning: probabilitySlaterTypeMove = ";
	  cerr << probabilitySlaterTypeMove << endl;
	  probabilitySlaterTypeMove = 0.0;
	}

      double probabilityGaussianTypeMove = 1.0 - probabilitySlaterTypeMove;
      
      // Update the greens function
      
      double distance1Sq = 0.0;
      double temp;
      
      for(int i=0; i<3; i++)
        {
          temp = OriginalWalker->R(electron,i) -
	    ( Input->Molecule.Atom_Positions(nearestNucleus,i) +
	      radialCoordinate * radialUnitVector(i) + zCoordinate * zUnitVector(i) );
                          
          distance1Sq += temp*temp;
        }
        
      double ga = 2*pi*tau;
      double gb = -1.5;
      double gc = -distance1Sq/(2*tau);
      
      GaussianGF=QMCDouble(probabilityGaussianTypeMove,ga,gb,gc);

      /*
      if(probabilitySlaterTypeMove < 1e-20)
	{
	  GF *= GaussianGF;
	  continue;
	}
      */

      double distance2Sq = 0.0;
      
      for(int i=0; i<3; i++)
        {
          temp = OriginalWalker->R(electron,i) - \
                        Input->Molecule.Atom_Positions(nearestNucleus,i);
                        
          distance2Sq += temp*temp;
        }
        
      double sk = probabilitySlaterTypeMove/pi;
      double sa = expParam;
      double sb = 3.0;
      double sc = -2.0*expParam*sqrt(distance2Sq);
      
      SlaterGF = QMCDouble(sk,sa,sb,sc);

      // The addition and multiplication operations here have been causing a
      // a lot of problems.
      OneE = SlaterGF + GaussianGF;

      GF *= OneE;
    }
  return GF;
}

QMCDouble QMCWalker::calculateReverseGreensFunctionUmrigar93AcceleratedSampling()
{
  int whichE = walkerData.whichE;
  bool moved;

  double cos_tm = cos(globalInput.flags.accel_tm);
  double delta  = globalInput.flags.accel_delta;

  Array1D<double> riA_hat(3);
  Array1D<double> rfA_hat(3);

  Array1D<double> Fx(3);
  Array2D<double> & Fk = TrialWalker->walkerData.modifiedGradPsiRatio;
  //Array2D<double> & Fk = TrialWalker->walkerData.gradPsiRatio;

  QMCDouble GF(1.0);
  for(int electron=0; electron<Input->WF.getNumberElectrons(); electron++)
    {      
      moved = true;
      if(whichE != -1 && electron != whichE)
	moved = false;
      if(!moved) continue;

      // Find the nearest nucleus to this electron
      int nearestNucleus = Input->Molecule.findClosestNucleusIndex(TrialWalker->walkerData.riA,electron);
      double Z = Input->Molecule.Z(nearestNucleus);
      double rkf = OriginalWalker->walkerData.riA(electron,nearestNucleus);
      double riA = TrialWalker->walkerData.riA(electron,nearestNucleus);
      double rav = (riA + rkf)/2.0;
      for(int xyz=0; xyz<3; xyz++)
	{
	  riA_hat(xyz) = TrialWalker->walkerData.riA_uvec(electron,nearestNucleus,xyz);
	  rfA_hat(xyz) = OriginalWalker->walkerData.riA_uvec(electron,nearestNucleus,xyz);
	}      

      /*
	First, calculate the new radius. To do this, we need to fit
	parameters a and zeta so that the U function matches log derivatives.
      */
      double Frk = 0.0;
      for(int xyz=0; xyz<3; xyz++)
	Frk += riA_hat(xyz)*Fk(electron,xyz);

      //Get zeta and a
      double zeta, a, I;
      MathFunctions::fitU(Z,Frk,riA,delta,zeta,a,I);

      if(rkf > riA*delta || rkf < riA/delta){
	printf("Warning: rkf is outside range! rkf=%20.10f l=%20.10f u=%20.10f\n",
	       rkf,riA/delta,riA*delta);
	GF *= QMCDouble(0.0);
	return GF;
      }

      /*
	Obtain the new theta, zenith angle
      */
      double tM, cos_tM;
      cos_tM = cos_tm - (1.0 + cos_tm)/(1.0 + Z*Z*rav*rav);
      tM = acos(cos_tM);

      // This should be exactly the same as before, if nearestNucleus is the same
      double ri_dot_rf = riA_hat * rfA_hat;
      double tkf;
      if(fabs(ri_dot_rf-1.0) < 1e-10) tkf = 0.0;
      else tkf = acos(riA_hat * rfA_hat);

      if(tkf > tM){
	printf("Warning: tkf is outside range! tkf = %20.10f tM=%20.10f\n",tkf,tM);
	GF *= QMCDouble(0.0);
	return GF;
      }

      double s_tkf = sin(tkf);

      /*
	Obtain the new phi, azimuth angle
      */
      double Fxk = 0.0;
      for(int i=0; i<3; i++)
        {
          Fx(i) = Fk(electron,i) - Frk*riA_hat(i);
          Fxk += Fx(i)*Fx(i);
        }
      Fxk = sqrt(Fxk);
      double Fpk = Fxk*rav*s_tkf;
      double pkf;
      double rf_dot_Fx = rfA_hat * Fx;
      if(fabs(rf_dot_Fx) < 1e-10) pkf = 0.5*pi;
      else pkf = acos((rfA_hat * Fx)/(s_tkf*Fxk));

      /*
      //It doesn't matter whether it's pkf or -pkf
      //but this will help verify that we have the correct coordinates
      double c_tkf = cos(tkf);
      if((rfA_hat.cross(Fx) * riA_hat) > 0.0) pkf *= -1;
      Array1D<double> newPosition(3);
      for(int xyz=0; xyz<3; xyz++)
	newPosition(xyz) = c_tkf*riA_hat(xyz) + s_tkf*Fx(xyz)/Fxk;
      newPosition.rotate(riA_hat,pkf);
      newPosition *= rkf;
      for(int xyz=0;xyz<3;xyz++)
	printf("%i pkf %20.10f diff %20.10f \n",xyz,pkf,
	       newPosition(xyz)-rfA_hat(xyz)*rkf);
      //*/
      //printf("rkf=%20.10f tkf=%20.10f pkf=%20.10f\n",rkf,tkf,pkf); 

      double S = 1.0;
      S *= (1.0 + a*rkf)*exp(-zeta*rkf);
      S *= 1.0 + min(Fpk,1.0)*cos(pkf);
      S *= pow(rkf,-1.5);
      S = fabs(S);

      I *= (1.0 - cos_tM);
      I *= 2*pi;

      if( fabs(I) < 1e-50)
	GF *= 0.0;
      else
	GF *= QMCDouble(S/I);

      if(GF.isNotValid()){
	printf("moved electron=%i age=%i\n",electron,age);
	printf("a=%20.10f zeta=%20.10f i=%20.10f o=%20.10f cos_tM=%20.10f\n",a,zeta,riA/delta,riA*delta,cos_tM);  
	printf("Frk=%20.10f Fxk=%20.10f Fpk=%20.10f\n",Frk,Fxk,Fpk);
	printf("rkf=%20.10f tkf=%20.10f pkf=%20.10f\n",rkf,tkf,pkf); 
	printf("S=%20.10e I=%20.10e T=%20.10e\n",S,I,(double)(GF));
      }
    }
  return GF;
}

void QMCWalker::reweight_walker()
{
  dW = 0.0;
  double S_trial = 0.0;
  double S_original = 0.0;
  double trialEnergy    = TrialWalker->walkerData.localEnergy;
  double originalEnergy = OriginalWalker->walkerData.localEnergy;
  double p = TrialWalker->getAcceptanceProbability();
  double q = 1.0 - p;      
      
  // determine the weighting factor dW so that the new weight = weight*dW
  
  bool weightIsNaN = false;

  if( Input->flags.run_type == "variational")
    {
      // Keep weights constant for VMC
      dW = 1.0;
      return;
    }

  if( IeeeMath::isNaN(trialEnergy) )
    {
      cerr << "WARNING: trial energy = ";
      cerr << TrialWalker->walkerData.localEnergy << "!" << endl;
      weightIsNaN = true;
    }
  else if( Input->flags.energy_modification_type == "none" )
    {
      S_trial    = Input->flags.energy_trial - trialEnergy;
      S_original = Input->flags.energy_trial - originalEnergy;
    }
  else
    {
      Array2D<double> * tMGPR =
	& TrialWalker->walkerData.modifiedGradPsiRatio;
      Array2D<double> * tGPR = & TrialWalker->walkerData.gradPsiRatio;
      Array2D<double> * oMGPR =
	& OriginalWalker->walkerData.modifiedGradPsiRatio;
      Array2D<double> * oGPR = & OriginalWalker->walkerData.gradPsiRatio;

      //Limit[ Vbar/V , R->nucleus ] = 1
      //Limit[ Vbar/V , R->node ]    = 0
      double lengthGradTrialModified =
	sqrt((*tMGPR).dotAllElectrons(*tMGPR));
      double lengthGradTrialUnmodified =
	sqrt((*tGPR).dotAllElectrons(*tGPR));
      double lengthGradOriginalModified =
	sqrt((*oMGPR).dotAllElectrons(*oMGPR));
      double lengthGradOriginalUnmodified =
	sqrt((*oGPR).dotAllElectrons(*oGPR));
      
      if(Input->flags.energy_modification_type=="modified_umrigar93")
	{
	  S_trial = (Input->flags.energy_trial - trialEnergy)*
	    (lengthGradTrialModified/lengthGradTrialUnmodified);
	  
	  S_original = (Input->flags.energy_trial - originalEnergy) *
	    (lengthGradOriginalModified/lengthGradOriginalUnmodified);
	}
      else if( Input->flags.energy_modification_type == "umrigar93" )
	{
	  S_trial =
	    (Input->flags.energy_trial - Input->flags.energy_estimated)
	    +(Input->flags.energy_estimated - trialEnergy)*
	    (lengthGradTrialModified/lengthGradTrialUnmodified);
	  
	  S_original =
	    (Input->flags.energy_trial - Input->flags.energy_estimated)
	    +(Input->flags.energy_estimated - originalEnergy) *
	    (lengthGradOriginalModified/lengthGradOriginalUnmodified);
	}
      else
	{
	  cerr << "ERROR: unknown energy modification method!" << endl;
	  exit(0);
	}
    }

  double reweight = 0.0;
  if( Input->flags.walker_reweighting_method ==
      "simple_symmetric" )
    {
      // This is the "classical" reweighting factor most people use.
      // umrigar says this has larger timestep and statistical errors than
      // umrigar93_probability_weighted

      // Umrigar, Nightingale, and Runge JCP 99(4) 2865; 1993, eq 8
      reweight = 0.5*(S_trial+S_original);
    }
  else if( Input->flags.walker_reweighting_method ==
	   "simple_asymmetric" )
    {
      // supposedly between simple_symmetric and
      // umrigar93_probability_weighted in terms of its timestep error and
      // statistical performance

      // Umrigar, Nightingale, and Runge JCP 99(4) 2865; 1993, eq 23
      if( move_accepted )
	reweight = 0.5*(S_trial+S_original);
      else
	reweight = S_original;
    }
  else if( Input->flags.walker_reweighting_method ==
	   "umrigar93_probability_weighted" )
    {
      // Umrigar claims this has a small time step error and a small
      // statistical error compared to simple_asymmetric and
      // simple_symmetric

      // Umrigar, Nightingale, and Runge JCP 99(4) 2865; 1993, eq 25
      reweight = (p*0.5*(S_original+S_trial) + q*S_original);
    }
  else
    {
      cerr << "ERROR: unknown reweighting method!" << endl;
      exit(0);
    }

  dW = exp(reweight*Input->flags.dt_effective);
  setWeight( getWeight() * dW );

  //if(getWeight() <= 0.0)
  //  return;

  /*
    Below this line are the statistical manipulations.
  */
  //return;

  double rel_diff;
  if(move_accepted)
    {
      rel_diff = fabs( (TrialWalker->walkerData.localEnergy -
			Input->flags.energy_estimated_original)/
		       Input->flags.energy_estimated_original);
    }
  else
    {
      rel_diff = fabs( (OriginalWalker->walkerData.localEnergy -
			Input->flags.energy_estimated_original)/
		       Input->flags.energy_estimated_original);
    }

  if(iteration < -100)
    {
      /*
	We're equilibrating when iteration < 0.

	If we're equilibrating, we can be more aggressive about
	handling bad walkers.
      */
      
      // we probably don't need warnings for the first couple steps
      int steps = iteration + Input->flags.equilibration_steps;
      // the number of steps until warnings are printed
      int min_steps = 5;

      if(getWeight() > 2.0*Input->flags.branching_threshold)
	{
	  if(steps > min_steps)
	    {
	      cerr << "WARNING: Deleting heavy walker " << ID(move_accepted,false);
	      cerr.flush();
	    }
	  setWeight(0.0);
	  return;
	}

      if(getWeight() < 0.1)
	{
	  if(steps > min_steps)
	    {
	      cerr << "WARNING: Deleting light walker " << ID(move_accepted,false);
	      cerr.flush();
	    }
	  setWeight(0.0);
	  return;
	}

      if(dW > 1.5 || dW < 0.5 || IeeeMath::isNaN(dW))
	{
	  if(steps > min_steps)
	    {
	      if(dW > 2.0 || dW < 0.25 || IeeeMath::isNaN(dW))
		{
		  cerr << "WARNING: Deleting walker with bad dW " << ID(move_accepted,false);
		  cerr << "       p = " << p << "; q = " << q;
		  cerr << "       energy_trial     = " << Input->flags.energy_trial << endl;
		  cerr << "       energy_estimated = " << Input->flags.energy_estimated << endl;
		  cerr << "       S_trial          = " << S_trial << endl;
		  cerr << "       S_original       = " << S_original << endl;
		  cerr << "TrialWalker:" << endl << TrialWalker->walkerData;
		  cerr << "OriginalWalker:" << endl << OriginalWalker->walkerData;
		  cerr << endl;
		}
	      else
		{
		  //cerr << "WARNING: Deleting fast growth walker " << ID(move_accepted);
		}
	      cerr.flush();
	    }
	  setWeight(0.0);
	  return;
	}

      double rel_cutoff = Input->flags.rel_cutoff;
      if(rel_diff > rel_cutoff && steps > min_steps)
	{
	  /*
	    A lot of walkers when initialized will have bad energies. So we
	    want to give them a couple of chances to move first.

	    Since rel_diff uses the SCF energy, it will bias the energy
	    toward zero, so we want a tradeoff between biasing the energy
	    and not allowing horrible energies into the average.

	    rel_diff > 1 corresponds to a positive energy, clearly unacceptable.
	    Is there any argument suggesting that a Local Energy can never be +'ve?
	   */
	  if(rel_diff > rel_cutoff && steps > 2*min_steps){
	    cerr << "WARNING: Deleting walker with bad energy " << ID(move_accepted,false);
	    cerr.flush();
	  }
	  setWeight(0.0);
	  return;
	}

      int age_cutoff = 30;
      if(age >= age_cutoff)
	{
	  // I don't see how a warning here could be very useful...
	  if(steps > age_cutoff && false)
	    {
	      cerr << "WARNING: Deleting aged walker " << ID(move_accepted,false);
	      cerr.flush();
	    }
	  setWeight(0.0);
	  return;
	}
    }
  else
    {
      /*
	We're not equillibrating anymore, so we want to be more
	careful. We'll only delete walkers we know will mess up
	the calculation.

	How aggressive can we be?
      */
      if(dW > 4.0 || IeeeMath::isNaN(dW))
	{
	  stringstream strm;	  
	  strm << "ERROR: Deleting fast growth walker " << ID(move_accepted,true);
	  int wi = 25;
	  int pr = 15;

	  strm.width(10);
	  strm << " p = "
	       << setw(wi) << setprecision(pr)
	       << fixed << p;
	  strm.width(10);
	  strm << "  q = "
	       << setw(wi) << setprecision(pr)
	       << fixed << q;
	  strm.width(10);
	  strm << endl;

	  strm.width(10);
	  strm << "Trial E = "
	       << setw(wi) << setprecision(pr)
	       << fixed << Input->flags.energy_trial;
	  strm.width(10);
	  strm << "  Est E = "
	       << setw(wi) << setprecision(pr)
	       << fixed << Input->flags.energy_estimated;
	  strm.width(10);
	  strm << endl;

	  strm.width(10);
	  strm << "reweight= "
	       << setw(wi) << setprecision(pr)
	       << scientific << reweight;
	  strm.width(10);
	  strm << " Eff dt = "
	       << setw(wi) << setprecision(pr)
	       << scientific << Input->flags.dt_effective;
	  strm.width(10);
	  strm << endl;

	  strm.width(10);
	  strm << "S_trial = "
	       << setw(wi) << setprecision(pr)
	       << scientific << S_trial;
	  strm.width(10);
	  strm << " S_orig = "
	       << setw(wi) << setprecision(pr)
	       << scientific << S_original;
	  strm.width(10);
	  strm << endl;


	  strm << "TrialWalker:" << endl << TrialWalker->walkerData;
	  strm << "OriginalWalker:" << endl << OriginalWalker->walkerData;
	  strm << endl;
	  cerr << strm.str();
	  cerr.flush();
	  setWeight(0.0);
	  return;
	}

      if(getWeight() > 50.0)
	{
	  cerr << "ERROR: Deleting heavy walker " << ID(move_accepted,true);
	  cerr.flush();
	  setWeight(0.0);
	  return;
	}

      if(rel_diff > globalInput.flags.rel_cutoff)
	{
	  cerr << "ERROR: Deleting walker with bad energy " << ID(move_accepted,!true);
	  cerr.flush();
	  setWeight(0.0);
	  return;
	}

    }
}

bool QMCWalker::branchRecommended()
{
  if(globalInput.flags.limit_branching == 0)
    return true;

  /*
    This function will be queried before a branch, meaning the walker is a
    candidate for branching. Since branching is purely an
    efficiency issue, we can do whatever we want and shouldn't have to worry about
    ruining detailed balance.

    We're equilibrating when iteration < 0, so we can be more aggressive in penalizing
    bad walkers... On the other hand, since we're deleting a lot of walkers during
    equilibration, maybe we want to make branching easier to replenish the loses.

    The difficulty is trying to get some warnings, but not too many. If we didn't
    delete all the bad walkers during equilibration, then "badness" can propagate
    through the calculation.
    
    By the time this function is called, OriginalWalker == TrialWalker
  */
  /*
    Some ideas to try:
    1) No branching until you've moved at least once.
    2) No branching if LE - TrialEnergy is outside some range (relative or abs?).
    3) A walker can't contribute to the energy average if its weight is too large
    4) Walkers can only be deleted if they haven't moved (and then raise the threshold?)
    5) 
   */
  bool shouldRecommend = true;
  string shouldWarn = "";
  //int aged = age - 10;

  if(age > 4)
    {
      shouldRecommend = false;
      if(age >= Input->flags.old_walker_acceptance_parameter && age%10 == 0)
	shouldWarn = "lazy";
    }
  if(dW > 1.5)
    {
      shouldRecommend = false;
      if(iteration%10 == 0 && getWeight() > 2.0*Input->flags.branching_threshold)
	shouldWarn = "dW";
      if(dW > 2.0)
	shouldWarn = "dW";
    }
  if(getWeight() > 2.0*Input->flags.branching_threshold)
    {
      shouldRecommend = false;
      if(iteration%50 == 0 &&
	 getWeight() > 2.0*Input->flags.branching_threshold &&
	 dW > 1.0)
	shouldWarn = "W";
    }

  double virial = -TrialWalker->walkerData.potentialEnergy/TrialWalker->walkerData.kineticEnergy;

  // Relatively large virial ratios do not appear to correlate with bad E_L.
  // It's hard to know how this would be a good indicator given rel_diff
  if(fabs(virial) < 1e-3)
    {
      shouldRecommend = false;
      //only warn when we arrive at a bad spot
      if(age == 0)
	shouldWarn = "virial";
    }

  double rel_diff = fabs( (TrialWalker->walkerData.localEnergy -
			   Input->flags.energy_estimated_original)/
			  Input->flags.energy_estimated_original);

  //should I use rel_diff > globalInput.flags.rel_cutoff here?
  if(rel_diff > globalInput.flags.rel_cutoff)
    {
      shouldRecommend = false;
      if(age == 0)
	shouldWarn = "rel E";
    }

  if(shouldWarn != "" && !shouldRecommend && iteration > 0)
    {
      if(shouldWarn == "lazy")
	{
	  /*
	  cerr << "WARNING: Not recommending (reason = " << shouldWarn << ") a branch for walker";
	  cerr << " age = " << setw(4) << age;
	  cerr << " (" << genealogy[0] << "::r" << Input->flags.my_rank << ")" << endl;
	  */
	}
      else
	{
	  cerr << "WARNING: Not recommending (reason = " << shouldWarn << ") a branch for walker";
	  cerr << ID(move_accepted,false);
	}
      cerr.flush();
    }

  return shouldRecommend;
}

string QMCWalker::ID(bool showTrial, bool verbose)
{
  /*
    The output only depends upon move_accepted if
    ID() is called before the processPropagation function completes
    since the last thing processPropagation does is 
    syncronize the Trial and Original walkerData.
   */
  QMCWalkerData * wd;
  if(showTrial)
    {
      wd = & TrialWalker->walkerData;
    }
  else
    {
      wd = & OriginalWalker->walkerData;
    }
  int w = 25;
  int p = 15;

  double virial = 0;
  if(fabs(wd->kineticEnergy) > 1e-30)
    virial = -wd->potentialEnergy/wd->kineticEnergy;

  /*
  double rel_diff = fabs( (wd->localEnergy -
		    Input->flags.energy_estimated_original)/
		   Input->flags.energy_estimated_original);
  */
  stringstream id;

  if(verbose){
    id << "(r" << Input->flags.my_rank << ")" << endl;
    id << setw(15) << "ID";
    id << setw(15) << "Born";
    id << setw(15) << "Age Branch";
    id << setw(15) << "Num Moves";
    id << setw(15) << "Weight" << endl;
    id << setw(15) << genealogy[0];
    id << setw(15) << birthIter[0];
    id << setw(15) << iteration-birthIter[0];
    id << setw(15) << numMoves[0];
    id << setw(15) << setprecision(5) << birthWeight[0] << endl;
    for(int i=1; i<numAncestors; i++)
      {
	if(genealogy[i] == -1)
	  break;
	id << setw(15) << genealogy[i];
	id << setw(15) << birthIter[i];
	id << setw(15) << birthIter[i-1]-birthIter[i];
	id << setw(15) << numMoves[i];
	id << setw(15) << setprecision(5) << birthWeight[i] << endl;
      }
  } else {
    id << "(" << genealogy[0];
    long numPrint = numAncestors < 5 ? numAncestors : 5;
    for(int i=1; i<numPrint; i++)
      {
	if(genealogy[i] == -1)
	  break;
	id << "<" << genealogy[i];
      }
    id << "::r" << Input->flags.my_rank << ")" << endl;
  }

  id << *wd;

  if(globalInput.flags.run_type != "variational")
    {
      id.width(10);
      id << "weight= "
	 << setw(w) << setprecision(p)
	 << scientific << getWeight();
      id.width(10);
      id << "   dW = "
	 << setw(w) << setprecision(p);
      if(fabs(dW) > 100)
	id << scientific << dW;
      else
	id << fixed << dW;
      id.width(10);
      id << endl;
    }

  id.width(10);
  id << " iter = "
     << setw(w) << setprecision(p)
     << fixed << iteration;
  id.width(10);
  id << "  age = "
     << setw(w) << setprecision(p)
     << scientific << age;
  id.width(10);
  id << endl;

  return id.str();
}

void QMCWalker::calculateMoveAcceptanceProbability(double GreensRatio)
{
  // This tells us the probability of accepting or rejecting a proposed move  
  double PsiRatio = TrialWalker->walkerData.psi/OriginalWalker->walkerData.psi;

  // calculate the probability of accepting the trial move
  double p = PsiRatio * PsiRatio * GreensRatio;

  /*
  printf("OP=%20.10e TP=%20.10e PsiRatio^2=%20.10e GR=%20.10e p=%20.10e\n",
	 (double)OriginalWalker->walkerData.psi,
	 (double)TrialWalker->walkerData.psi,
	 PsiRatio*PsiRatio,GreensRatio,p);
  //*/
  if( !(IeeeMath::isNaN(p)))
    {
      // 1.1 ^ 1000 ~= 1e41
      if(age - Input->flags.old_walker_acceptance_parameter > 1000)
	{
	  //setWeight(0.0);
	  p = 1.0;
	} 
      else if(age - Input->flags.old_walker_acceptance_parameter > 500)
	{
	  p = 1.0;
	} 
      else if(age - Input->flags.old_walker_acceptance_parameter > 0)
	{
	  p *= pow(1.1,age-Input->flags.old_walker_acceptance_parameter);
	}
      
    } else {
      // if the aratio is NaN then reject the move
      numWarnings++;
      cerr << "WARNING: Rejecting trial walker with NaN p! " << ID(true,false);
      cerr << "         PsiRatio    = " << PsiRatio << endl;
      cerr << "         GreensRatio = " << GreensRatio << endl;
      p = 0.0;
      // The energies are assigned zero because the later multiplication by 
      // p=0 doesn't result in 0 contribution.
      TrialWalker->walkerData.zero();
    }
    
  // if the trial position is singular reject the move
  if( TrialWalker->isSingular() )
    {
      numWarnings++;
      cerr << "WARNING: Rejecting singular trial walker! " << ID(true,false);
      p = 0.0;
      TrialWalker->walkerData.zero();
    }

  if(numWarnings >= 10){
    /*
      I've had problems with exploding files. It's usually an indicator of a
      programming error somewhere.
    */
    cerr << "ERROR: There have been " << numWarnings << " warnings, so I'm quitting." << endl;

#ifdef PARALLEL
    MPI_Abort(MPI_COMM_WORLD,1);
#endif
    exit(0);
  }
    
  // Fixed Node Condition
  if( PsiRatio < 0 && Input->flags.run_type == "diffusion" )
    {
      // Fixed Node Condition
      p = 0.0;
    }
    
  // Limit the probability to 1
  if( p > 1.0)
    {
      p = 1.0;
    }

  // Set the probabilities of taking the trial move or keeping the current
  // move
  TrialWalker->setAcceptanceProbability(p);
}

void QMCWalker::acceptOrRejectMove()
{
  if( TrialWalker->getAcceptanceProbability() > ran.unidev() )
    {
      // accept the move
      ageMoved = age;
      age = 0;
      move_accepted = true;
      numMoves[0]++;
    }
  else
    {
      // reject the move
      //age measures the number of steps for which the move was not accepted
      age++;
      move_accepted = false;
    }

  //int aged = age - Input->flags.old_walker_acceptance_parameter;
  if( iteration > 0)
    {
      if(ageMoved > Input->flags.old_walker_acceptance_parameter + 15)
	{
	  //cerr << "WARNING: Old walker moved after " << ageMoved << " iterations " << ID(false,true);
	  cerr << "WARNING: Old walker moved after " << ageMoved << " iterations." << endl;
	  cerr << "Original Walker ";
	  cerr << ID(false,true);
	  cerr << "Trial Walker ";
	  cerr << ID(true,false);
	  //cerr << "WARNING: Old walker " moved after " << ageMoved << " iterations." << endl;
	}
      /*
      if(aged >= 50 && aged%50 == 0)
	{
	  cerr << "WARNING: Aged walker " << ID();
	}
      */
    }
}

void QMCWalker::createChildWalkers()
{
  /*
    We only ever need 2 walker objects, so AGA modified the code
    so that we don't use 3 walker objects anymore. In order to preserve
    the rest of the code, the TrialWalker pointer just points to "this" now,
    this results in fewer memory copies (since the move is likely to be accepted)
    than if OriginalWalker pointed to "this".

    Correspondingly, I eliminated redundant data.
  */
  if( OriginalWalker == 0 )
    {
      OriginalWalker = new QMCWalker();
    }

  //We create a copy of the data we need
  OriginalWalker->R                     = R;
  OriginalWalker->dR2                   = dR2;
  OriginalWalker->AcceptanceProbability = AcceptanceProbability;

  if(globalInput.flags.one_e_per_iter == 0)
    {
      OriginalWalker->walkerData            = walkerData;
    } else {
      //when moving one electron at a time, a full copy is expensive
      OriginalWalker->walkerData.partialCopy(walkerData);
    }

  //And use a psuedonym
  TrialWalker = this;
}

void QMCWalker::initialize(QMCInput *INPUT)
{
  Input = INPUT;

  int numFW = Input->flags.future_walking.size();
  int numElectrons  = Input->WF.getNumberElectrons();
  int numDimensions;
  int numNucForceDim1 = -1;
  int numNucForceDim2 = -1;
  
  if(Input->flags.trial_function_type == "restricted" ||
     Input->flags.trial_function_type == "unrestricted"){     
    numDimensions = 3;
  } else {
    //We're hijacking the Nbasisfunc parameter input
    numDimensions = Input->flags.Nbasisfunc;
  }

  if(Input->flags.nuclear_derivatives != "none")
    {
      if(Input->flags.nuclear_derivatives != "bin_force_density")
	{
	  numNucForceDim1 = Input->Molecule.getNumberAtoms();
	  numNucForceDim2 = 3;
	} else {
	  numNucForceDim1 = QMCNuclearForces::getNumBins();
	  numNucForceDim2 = 1;
	}
      
      fwNuclearForces.allocate(numNucForceDim1,numNucForceDim2);
      
      for(int d1=0; d1<fwNuclearForces.dim1(); d1++)
	for(int d2=0; d2<fwNuclearForces.dim2(); d2++)
	  fwNuclearForces(d1,d2).allocate(numFW,2);
    }

  R.allocate(numElectrons,numDimensions);

  walkerData.initialize(numDimensions,
			numNucForceDim1,numNucForceDim2);

  numFWSteps.resize(numFW);
  
  isCollectingFWResults.allocate(numFW,2);
  fwNormalization.allocate(numFW,2);
  fwR12.allocate(numFW,2);
  fwR2.allocate(numFW,2);
  fwiR12.allocate(numFW,2);
  fwiR.allocate(numFW,2);
  fwEnergy.allocate(numFW,2);
  fwKineticEnergy.allocate(numFW,2);
  fwKineticEnergy_grad.allocate(numFW,2);
  fwPotentialEnergy.allocate(numFW,2);
  resetFutureWalking();

  //initialize acceptance probability
  setAcceptanceProbability(0.0);
}

void QMCWalker::calculateElectronDensities(double max_pair_distance, double dr,
                          Array1D<double> &pll_spin, Array1D<double> &opp_spin,
		                     Array1D< Array1D<double> > &alpha_density,
                                      Array1D< Array1D<double> > &beta_density)
{
  int nalpha = Input->WF.getNumberElectrons(true);
  int nbeta = Input->WF.getNumberElectrons(false);
  
  double dist = 0.0;
  int index = 0;

  // We calculate the distance between each same spin pair and record them in 
  // the histogram.
  if (nalpha > 1)
    {
      for (int i=0; i<nalpha-1; i++)
        for (int j=i+1; j<nalpha; j++)
          {
            dist = sqrt( (R(i,0) - R(j,0)) * (R(i,0) - R(j,0)) +
                         (R(i,1) - R(j,1)) * (R(i,1) - R(j,1)) +
                         (R(i,2) - R(j,2)) * (R(i,2) - R(j,2)) );
	    if (dist < max_pair_distance)
	      {
		index = int(dist/dr);
		pll_spin(index) += weight;
	      }
          }
    }
    
  if (nbeta > 1)
    {
      for (int i=nalpha; i<nalpha+nbeta-1; i++)
        for (int j=i+1; j<nalpha+nbeta; j++)
          {
            dist = sqrt( (R(i,0) - R(j,0)) * (R(i,0) - R(j,0)) +
                         (R(i,1) - R(j,1)) * (R(i,1) - R(j,1)) +
                         (R(i,2) - R(j,2)) * (R(i,2) - R(j,2)) );
	    if (dist < max_pair_distance)
	      {
		index = int(dist/dr);
		pll_spin(index) += weight;
	      }
          }
    }

  // We calculate the distance between each opposite spin pair and record them
  // in the histogram.
  if (nalpha > 0 && nbeta > 0)
    {
      for (int i=0; i<nalpha; i++)
        for (int j=nalpha; j<nalpha+nbeta; j++)
          {
            dist = sqrt( (R(i,0) - R(j,0)) * (R(i,0) - R(j,0)) +
                         (R(i,1) - R(j,1)) * (R(i,1) - R(j,1)) +
                         (R(i,2) - R(j,2)) * (R(i,2) - R(j,2)) );
	    if (dist < max_pair_distance)
	      {
		index = int(dist/dr);
		opp_spin(index) += weight;
	      }
          }
    }

  for (int i=0; i<Input->flags.Natoms; i++)
    {
      string NucleusType = Input->Molecule.Atom_Labels(i);

      int nuc_index = -1;
      for (int j=0; j<Input->Molecule.NucleiTypes.dim1(); j++)
	if (Input->Molecule.NucleiTypes(j) == NucleusType)
	  {
	    nuc_index = j;
	    break;
	  }

      double nuc_x = Input->Molecule.Atom_Positions(i,0);
      double nuc_y = Input->Molecule.Atom_Positions(i,1);
      double nuc_z = Input->Molecule.Atom_Positions(i,2);
      
      for (int k=0; k<nalpha; k++)
	{
	  dist = sqrt( (R(k,0) - nuc_x) * (R(k,0) - nuc_x) + 
		       (R(k,1) - nuc_y) * (R(k,1) - nuc_y) +
		       (R(k,2) - nuc_z) * (R(k,2) - nuc_z) );
	  if (dist < max_pair_distance)
	    {
	      index = int(dist/dr);
	      (alpha_density(nuc_index))(index) += weight;
	    }
	}

      for (int l=nalpha; l<nalpha+nbeta; l++)
	{
	  dist = sqrt( (R(l,0) - nuc_x) * (R(l,0) - nuc_x) +
		       (R(l,1) - nuc_y) * (R(l,1) - nuc_y) +
		       (R(l,2) - nuc_z) * (R(l,2) - nuc_z) );
	  if (dist < max_pair_distance)
	    {
	      index = int(dist/dr);
	      (beta_density(nuc_index))(index) += weight;
	    }
	}
    }
}

void QMCWalker::calculatePllCorrelationDiagram(int coord, double min, 
		    double max, Array1D< Array1D<double> > &CorrelationDiagram)
{
  int nalpha = Input->WF.getNumberElectrons(true);
  int nbeta = Input->WF.getNumberElectrons(false);
  
  double dr = (max-min)/CorrelationDiagram.dim1();
  int index1 = -1;
  int index2 = -1;

  for (int i=0; i<nalpha-1; i++)
    for (int j=i+1; j<nalpha; j++)
      if ( (R(i,coord) >= min) && (R(i,coord) <=max) )
	if ( (R(j,coord) >= min) && (R(j,coord) <= max) )
	  {
	    index1 = int((R(i,coord)-min)/dr);
	    index2 = int((R(j,coord)-min)/dr);
	    (CorrelationDiagram(index1))(index2) += weight;
	  }

  for (int i=nalpha; i<nalpha+nbeta-1; i++)
    for (int j=i+1; j<nalpha+nbeta; j++)
      if ( (R(i,coord) >= min) && (R(i,coord) <=max) )
	if ( (R(j,coord) >= min) && (R(j,coord) <= max) )
	  {
	    index1 = int((R(i,coord)-min)/dr);
	    index2 = int((R(j,coord)-min)/dr);
	    (CorrelationDiagram(index1))(index2) += weight;
	  }
}

void QMCWalker::calculateOppCorrelationDiagram(int coord, double min, 
		    double max, Array1D< Array1D<double> > &CorrelationDiagram)
{
  int nalpha = Input->WF.getNumberElectrons(true);
  int nbeta = Input->WF.getNumberElectrons(false);
  
  double dr = (max-min)/CorrelationDiagram.dim1();
  int index1 = 0;
  int index2 = 0;

  for (int i=0; i<nalpha; i++)
    for (int j=nalpha; j<nalpha+nbeta; j++)
      if ( (R(i,coord) >= min) && (R(i,coord) <=max) )
	if ( (R(j,coord) >= min) && (R(j,coord) <= max) )
	  {
	    index1 = int((R(i,coord)-min)/dr);
	    index2 = int((R(j,coord)-min)/dr);
	    (CorrelationDiagram(index1))(index2) += weight;
	  }
}

void QMCWalker::toXML(ostream& strm)
{
  strm << "<QMCWalker>" << endl;
  strm << "\t<Position>" <<endl;
  for(int ep=0; ep<R.dim1(); ep++)
    {
      strm << "\t\t";
      for(int j=0;j<R.dim2();j++)
	strm << R(ep,j) << "\t";
      strm << endl;
    }
  strm << "\t</Position>" << endl;

  strm << "\t<Weight>\n\t\t" << getWeight() << "\n\t</Weight>" << endl;
  strm << "\t<Age>\n\t\t" << getAge() << "\n\t</Age>" << endl;
  strm << "</QMCWalker>\n" << endl;
}

bool QMCWalker::readXML(istream& strm, QMCFunctions & QMF)
{
  string temp;
  strm >> temp;
  
  if (temp != "<QMCWalker>")
    return false;

  // Read position
  strm >> temp;
  if (temp != "<Position>")
    return false;
  
  for(int ep=0; ep<R.dim1(); ep++)
    for(int j=0;j<R.dim2();j++)
      strm >> R(ep,j);
    
  strm >> temp;
  if (temp != "</Position>")
    return false;
  
  // Read weight
  strm >> temp;
  if (temp != "<Weight>")
    return false;
  strm >> temp;
  weight = atof(temp.c_str());
  strm >> temp;
  if (temp != "</Weight>")
    return false;

  // Read age
  strm >> temp;
  if (temp != "<Age>")
    return false;
  strm >> temp;
  age = atoi(temp.c_str());
  strm >> temp;
  if (temp != "</Age>")
    return false;

  strm >> temp;
  if (temp != "</QMCWalker>")
    return false;

  walkerData.updateDistances(R);
  QMF.evaluate(R,walkerData);
  newID();

  return true;
}

void QMCWalker::initializeWalkerPosition(QMCFunctions & QMF)
{
  QMCInitializeWalker * IW =
    QMCInitializeWalkerFactory::initializeWalkerFactory(Input,
        Input->flags.walker_initialization_method);
  
  setR(IW->initializeWalkerPosition());
  QMF.evaluate(R,walkerData);
  
  int initilization_try = 1;
  while( isSingular() )
    {
      cerr << "Regenerating Walker..." << endl;
      
      if( initilization_try > 10 )
        {
          cerr << "ERROR: 10 consecutive singular configurations while "
          << "trying to initilize walker!" << endl;
          exit(0);
        }

      setR(IW->initializeWalkerPosition());
      QMF.evaluate(R,walkerData);
      initilization_try++;
    }
    
  delete IW;
  IW = 0;
}

double QMCWalker::getWeight()
{
  return weight;
}

void QMCWalker::setWeight(double val)
{
  weight = val;
}

int QMCWalker::getAge()
{
  return age;
}

double QMCWalker::getAcceptanceProbability()
{
  return AcceptanceProbability;
}

void QMCWalker::setAcceptanceProbability(double val)
{
  AcceptanceProbability = val;
}

Array2D<double> * QMCWalker::getR()
{
  return &R;
}

bool QMCWalker::setR(Array2D<double> temp_R)
{
  bool ok = true;
  for(int i=0; i<temp_R.dim1(); i++)
    for(int j=0; j<temp_R.dim2(); j++)
      if(IeeeMath::isNaN(temp_R(i,j)) || temp_R(i,j) == 0 || fabs(temp_R(i,j)) > 500.0 )
	{
	  //cout << "Error: setting electron position (" << i << "," << j << ") = " << temp_R(i,j) << endl;
	  ok = false;
	}
  if(ok)
    {
      R = temp_R;
      walkerData.updateDistances(R);
    }
  return ok;
}

QMCWalkerData* QMCWalker::getWalkerData()
{
  return &walkerData;
}

void QMCWalker::calculateObservables()
{
  double p = TrialWalker->AcceptanceProbability;
  double q = 1.0 - p;
  
  // Calculate the Energy ...
  localEnergy = p * TrialWalker->walkerData.localEnergy +
                q * OriginalWalker->walkerData.localEnergy;

  // Calculate the kinetic energy...
  kineticEnergy = p * TrialWalker->walkerData.kineticEnergy +
                  q * OriginalWalker->walkerData.kineticEnergy;
  
  // Calculate the potential energy
  potentialEnergy = p * TrialWalker->walkerData.potentialEnergy +
                    q * OriginalWalker->walkerData.potentialEnergy;
  
  // Calculate the ne and ee potential energy
  neEnergy = p * TrialWalker->walkerData.neEnergy + 
             q * OriginalWalker->walkerData.neEnergy;
  eeEnergy = p * TrialWalker->walkerData.eeEnergy + 
             q * OriginalWalker->walkerData.eeEnergy;


  if(globalInput.cs_Parameters.dim1() > 1)
    {
      cs_Energies.allocate(TrialWalker->walkerData.cs_LocalEnergy.dim1());
      cs_Weights.allocate(TrialWalker->walkerData.cs_Weights.dim1());
      
      if(OriginalWalker->walkerData.cs_LocalEnergy.dim1() ==
	 TrialWalker->walkerData.cs_LocalEnergy.dim1())
	{
	  for(int cs=0; cs<cs_Energies.dim1(); cs++)
	    {
	      cs_Energies(cs) = p * TrialWalker->walkerData.cs_LocalEnergy(cs) +
		q * OriginalWalker->walkerData.cs_LocalEnergy(cs);
	      cs_Weights(cs) = p * TrialWalker->walkerData.cs_Weights(cs) +
		q * OriginalWalker->walkerData.cs_Weights(cs);
	      //cs_Weights(cs) = TrialWalker->walkerData.cs_Weights(cs);
	    }
	} else {
	  cs_Energies = 0.0;
	  cs_Weights = 0.0;
	}
    }
  
  if(globalInput.flags.calculate_Derivatives == 1)
    {
      p3_xxa.allocate(TrialWalker->walkerData.p3_xxa.dim1());
      rp_a.allocate(TrialWalker->walkerData.rp_a.dim1());
      
      if(iteration == 1)
	{
	  /*
	    If we weren't calculating derivatives while equilibrating, then
	    we need this special case in order to be consistent.
	  */
	  for(int ai=0; ai<rp_a.dim1(); ai++)
	    {
	      p3_xxa(ai) = TrialWalker->walkerData.p3_xxa(ai);
	      
	      rp_a(ai) = TrialWalker->walkerData.rp_a(ai);
	    }
	} else {
	  for(int ai=0; ai<rp_a.dim1(); ai++)
	    {     
	      p3_xxa(ai) = p * TrialWalker->walkerData.p3_xxa(ai) +
		q * OriginalWalker->walkerData.p3_xxa(ai);
	      
	      rp_a(ai) = p * TrialWalker->walkerData.rp_a(ai) +
		q * OriginalWalker->walkerData.rp_a(ai);
	    }
	}
    }
  
  // Calculate the DistanceMovedAccepted this is the average distance
  // moved on a step
  distanceMovedAccepted = p * dR2;
  
  if (Input->flags.calculate_bf_density == 1)
    {
      for (int i=0; i<Input->WF.getNumberBasisFunctions(); i++)
        walkerData.chiDensity(i) =
          p * TrialWalker->walkerData.chiDensity(i) +
          q * OriginalWalker->walkerData.chiDensity(i);
    }
  
  if(Input->flags.nuclear_derivatives != "none")
    {
      for (int d1=0; d1<walkerData.nuclearDerivatives.dim1(); d1++)
        for (int d2=0; d2<walkerData.nuclearDerivatives.dim2(); d2++)
          walkerData.nuclearDerivatives(d1,d2) =
            p * TrialWalker->walkerData.nuclearDerivatives(d1,d2) +
            q * OriginalWalker->walkerData.nuclearDerivatives(d1,d2);
    }
  
  //eventually, i'll remove the r2 and r12 measurement. but for Helium, it's
  //nice to use when verifying FW
  double calcR12_T=0, calcR12_O=0;
  double calcR1_T=0, calcR1_O=0;
  double calcR2_T=0, calcR2_O=0;
  r2 = 0.0;
  for(int i=0; i<R.dim2(); i++)
    {
      //assuming the nucleus is at the origin
      //also, we're measuring the same thing for both electrons
      //and averaging
      calcR1_T += TrialWalker->R(0,i)*TrialWalker->R(0,i);
      calcR2_T += TrialWalker->R(1,i)*TrialWalker->R(1,i);
      calcR1_O += OriginalWalker->R(0,i)*OriginalWalker->R(0,i);
      calcR2_O += OriginalWalker->R(1,i)*OriginalWalker->R(1,i);
      
      double tempT, tempO;
      tempT = TrialWalker->R(0,i)    - TrialWalker->R(1,i);
      tempO = OriginalWalker->R(0,i) - OriginalWalker->R(1,i);
      
      calcR12_T += tempT*tempT;
      calcR12_O += tempO*tempO;
    }
  r12  = p*sqrt(calcR12_T) + q*sqrt(calcR12_O);
  ir12 = p/sqrt(calcR12_T) + q/sqrt(calcR12_O);
  ir = (p*(1.0/sqrt(calcR1_T) + 1.0/sqrt(calcR2_T))
     +  q*(1.0/sqrt(calcR1_O) + 1.0/sqrt(calcR2_O)))/2.0;
  r2 = (p*(calcR1_T + calcR2_T)
     +  q*(calcR1_O + calcR2_O))/2.0;

  /*
    This section of code probably measure anything useful
  */
  double kineticEnergy_grad_O = 0.0;
  double kineticEnergy_grad_T = 0.0;
  kineticEnergy_grad_O = OriginalWalker->walkerData.gradPsiRatio.dotAllElectrons(OriginalWalker->walkerData.gradPsiRatio);
  kineticEnergy_grad_T = TrialWalker->walkerData.gradPsiRatio.dotAllElectrons(TrialWalker->walkerData.gradPsiRatio);
  kineticEnergy_grad = p*kineticEnergy_grad_T*TrialWalker->walkerData.psi
                     + q*kineticEnergy_grad_O*OriginalWalker->walkerData.psi;

  //This is the forward walking portion of the calculation as described in:
  //J. Casulleras and J. Boronat, Phys. Rev. B 52, 3654 (1995) aka "CB95"
  //this is what CB95 suggest
  //double aWeight = getWeight();
  //double pWeight = getWeight()/OriginalWalker->getWeight();
  
  //this produces the same results as CB95; the difference cancels off in the end
  double aWeight = 1.0;
  double pWeight = 1.0;
    
  //This happens in both Formula 15 and 16 from CB95
  for(int i=0; i<isCollectingFWResults.dim1(); i++)
    {
      for(int j=0; j<isCollectingFWResults.dim2(); j++)
	{
	  
	  if(Input->flags.future_walking[i] == 0)
	    continue;
	  
	  if( isCollectingFWResults(i,j) != DONE)
	    {
	      fwNormalization(i,j)     *= pWeight;
	      fwR12(i,j)               *= pWeight;
	      fwR2(i,j)                *= pWeight;
	      fwiR12(i,j)              *= pWeight;
	      fwiR(i,j)                *= pWeight;
	      fwEnergy(i,j)            *= pWeight;
	      fwKineticEnergy(i,j)     *= pWeight;
	      fwKineticEnergy_grad(i,j)*= pWeight;
	      fwPotentialEnergy(i,j)   *= pWeight;

	      /*
		//I probably commented out this section because
		//It doesn't work with harmonic oscillators...
	      if(Input->flags.nuclear_derivatives != "none")
		{
		  for(int d1=0; d1<fwNuclearForces.dim1(); d1++)
		    for(int d2=0; d2<fwNuclearForces.dim2(); d2++)
		      (fwNuclearForces(d1,d2))(i,j) *= pWeight;
		}
	      */
	    }
	  
	  if(isCollectingFWResults(i,j) == ACCUM )
	    {
	      //We are collecting results (Formula 15)
	      fwNormalization(i,j)   += aWeight;
	      fwR12(i,j)             += aWeight * r12;
	      fwR2(i,j)              += aWeight * r2;
	      fwiR12(i,j)            += aWeight * ir12;
	      fwiR(i,j)              += aWeight * ir;
	      fwEnergy(i,j)          += aWeight * localEnergy;
	      fwKineticEnergy(i,j)   += aWeight * kineticEnergy;
	      fwKineticEnergy_grad(i,j) += aWeight * kineticEnergy_grad;
	      fwPotentialEnergy(i,j) += aWeight * potentialEnergy;

	      /*
	      if(Input->flags.nuclear_derivatives != "none")
		{
		  for(int d1=0; d1<fwNuclearForces.dim1(); d1++)
		    for(int d2=0; d2<fwNuclearForces.dim2(); d2++)
		      (fwNuclearForces(d1,d2))(i,j) += aWeight * walkerData.nuclearDerivatives(d1,d2);
		}
	      */
	    }
	}
    }
}

void QMCWalker::calculateObservables( QMCProperties & props )
{
  if(getWeight() <= 0.0)
    return;

  if(walkerData.whichE != -1 &&
     walkerData.whichE != 0)
    {
      /*
	If we're moving electrons individually, we still only want
	to collect statistics once all electrons have moved.

	Is this necessary? If we are able to calculate serial
	correlation effectively, does it really matter? Or is there
	another issue?
      */
      return;
    }


  // Add the data from this walker to the accumulating properties
  
  /*
    What we want to measure is how bad the bad ones are.
    So, we want to set this cutoff high enough to cut out
    all the good to decent walkers.
   */
  if(iteration > 0 && ageMoved > Input->flags.old_walker_acceptance_parameter){
    props.walkerAge.newSample(ageMoved,1.0);
  }
  ageMoved = -1;


  props.weightChange.newSample(dW,1.0);
    
  double rel_diff = fabs( (localEnergy -
			   Input->flags.energy_estimated_original)/
			  Input->flags.energy_estimated_original);
  if(rel_diff < Input->flags.rel_cutoff)
    {
      // Calculate the Energy ...
      props.energy.newSample( localEnergy, getWeight() );
      
      // Calculate the Kinetic Energy  
      props.kineticEnergy.newSample( kineticEnergy, getWeight() );
      
      // Calculate the Potential Energy  
      props.potentialEnergy.newSample( potentialEnergy, getWeight() );
      
      // Calculate the ne and ee potential energy
      props.neEnergy.newSample( neEnergy, getWeight() );
      props.eeEnergy.newSample( eeEnergy, getWeight() );
    }
  else
    {
      if(age == 0)
	{
	  bool shouldWarn = true;
	  if(shouldWarn && iteration + Input->flags.equilibration_steps > 5){
	    cerr << "ERROR: Not including walker in average " << ID(move_accepted,false);
	    cerr << "   localEnergy = " << localEnergy << endl;
	    cerr.flush();
	  }
	}
    }

  // Calculate the Acceptance Probability ...
  props.acceptanceProbability.newSample( getAcceptanceProbability(), getWeight() );
  
  // Calculate the DistanceMovedAccepted this is the average distance
  // moved on a step
  props.distanceMovedAccepted.newSample( distanceMovedAccepted, getWeight() );
  
  // Calculate the DistanceMovedTrial this is the average step length for
  // a trial move
  props.distanceMovedTrial.newSample( dR2, getWeight() );
  
  // Calculate the log of the weight
  props.logWeights.newSample(getWeight(),1);
}

void QMCWalker::calculateObservables( QMCPropertyArrays & fwProps )
{

  double rel_diff = fabs( (localEnergy -
			   Input->flags.energy_estimated_original)/
			  Input->flags.energy_estimated_original);

  double fWeight = getWeight();

  // Add the data from this walker to the accumulating properties  
  for(int i=0; i<isCollectingFWResults.dim1(); i++)
    {
      numFWSteps[i]++;
      
      //This collects the non-forward walking results
      if(Input->flags.future_walking[i] == 0)
	{
	  // Calculate the nuclear forces      
	  if(Input->flags.nuclear_derivatives != "none")
	    for (int d1=0; d1<walkerData.nuclearDerivatives.dim1(); d1++)
	      for (int d2=0; d2<walkerData.nuclearDerivatives.dim2(); d2++)
		(fwProps.nuclearForces(i))(d1,d2).newSample( walkerData.nuclearDerivatives(d1,d2), fWeight );
	  
	  (fwProps.props[FW_It])(i).newSample(1.0, fWeight);
	  (fwProps.props[FW_R12])(i).newSample(r12, fWeight);
	  (fwProps.props[FW_R2])(i).newSample(r2, fWeight);
	  (fwProps.props[FW_iR12])(i).newSample(ir12, fWeight);
	  (fwProps.props[FW_iR])(i).newSample(ir, fWeight);
	  (fwProps.props[FW_TE])(i).newSample(localEnergy, fWeight);
	  (fwProps.props[FW_KE])(i).newSample(kineticEnergy, fWeight);
	  (fwProps.props[FW_KEg])(i).newSample(kineticEnergy_grad, fWeight);
	  (fwProps.props[FW_PE])(i).newSample(potentialEnergy, fWeight);

	  (fwProps.props[FW_R12_2])(i).newSample(r12*r12, fWeight);
	  (fwProps.props[FW_R2_2])(i).newSample(r2*r2, fWeight);
	  (fwProps.props[FW_iR12_2])(i).newSample(ir12*ir12, fWeight);
	  (fwProps.props[FW_iR_2])(i).newSample(ir*ir, fWeight);
	  (fwProps.props[FW_TE_2])(i).newSample(localEnergy*localEnergy, fWeight);
	  (fwProps.props[FW_KE_2])(i).newSample(kineticEnergy*kineticEnergy, fWeight);
	  (fwProps.props[FW_KEg_2])(i).newSample(kineticEnergy_grad*kineticEnergy_grad, fWeight);
	  (fwProps.props[FW_PE_2])(i).newSample(potentialEnergy*potentialEnergy, fWeight);
	  continue;
	}
      
      //this should be >=, but with > we can allow 1.0 for maxFWAsymp
      //so we actually do 1 more step than was input
      if(numFWSteps[i] > Input->flags.future_walking[i])
	{
	  int whichIsDone = isCollectingFWResults(i,0) == DONE ? 0 : 1;
	  int otherStage = (whichIsDone+1)%2;
	  isCollectingFWResults(i,whichIsDone) = ACCUM;
	  isCollectingFWResults(i,otherStage)  = ASYMP;
	  numFWSteps[i] = 0;  
	}
      else if(numFWSteps[i] == int(maxFWAsymp*Input->flags.future_walking[i]+0.5))
	{
	  //CB95 do not suggest how long to propagate the FW steps into the next
	  //block, so maxFWAsymp is the percentage of the next block to use.
	  //I've found that a value of maxFWAsymp=1 works best.

	  int whichIsDone = isCollectingFWResults(i,0) == ACCUM ? 1 : 0;
	  
	  if(isCollectingFWResults(i,whichIsDone) == DONE)
	    break;
	  isCollectingFWResults(i,whichIsDone) = DONE;
	  
	  double norm = 1.0/fwNormalization(i,whichIsDone);
	  (fwProps.props[FW_It])(i).newSample(fwNormalization(i,whichIsDone)/numFWSteps[i],fWeight);
	  (fwProps.props[FW_R12])(i).newSample(fwR12(i,whichIsDone)*norm,fWeight);
	  (fwProps.props[FW_R2])(i).newSample(fwR2(i,whichIsDone)*norm,fWeight);
	  (fwProps.props[FW_iR12])(i).newSample(fwiR12(i,whichIsDone)*norm,fWeight);
	  (fwProps.props[FW_iR])(i).newSample(fwiR(i,whichIsDone)*norm,fWeight);
	  (fwProps.props[FW_TE])(i).newSample(fwEnergy(i,whichIsDone)*norm,fWeight);
	  (fwProps.props[FW_KE])(i).newSample(fwKineticEnergy(i,whichIsDone)*norm,fWeight);
	  (fwProps.props[FW_KEg])(i).newSample(fwKineticEnergy_grad(i,whichIsDone)*norm,fWeight);
	  (fwProps.props[FW_PE])(i).newSample(fwPotentialEnergy(i,whichIsDone)*norm,fWeight);

	  (fwProps.props[FW_R12_2])(i).newSample(fwR12(i,whichIsDone)*norm*r12,fWeight);
	  (fwProps.props[FW_R2_2])(i).newSample(fwR2(i,whichIsDone)*norm*r2,fWeight);
	  (fwProps.props[FW_iR12_2])(i).newSample(fwiR12(i,whichIsDone)*norm*ir12,fWeight);
	  (fwProps.props[FW_iR_2])(i).newSample(fwiR(i,whichIsDone)*norm*ir,fWeight);
	  (fwProps.props[FW_TE_2])(i).newSample(fwEnergy(i,whichIsDone)*norm*localEnergy,fWeight);
	  (fwProps.props[FW_KE_2])(i).newSample(fwKineticEnergy(i,whichIsDone)*norm*kineticEnergy,fWeight);
	  (fwProps.props[FW_KEg_2])(i).newSample(fwKineticEnergy_grad(i,whichIsDone)*norm*kineticEnergy,fWeight);
	  (fwProps.props[FW_PE_2])(i).newSample(fwPotentialEnergy(i,whichIsDone)*norm*potentialEnergy,fWeight);
	  
	  // Calculate the nuclear forces      
	  if(Input->flags.nuclear_derivatives != "none")//what about binning?
	    for (int d1=0; d1<walkerData.nuclearDerivatives.dim1(); d1++)
	      for (int d2=0; d2<walkerData.nuclearDerivatives.dim2(); d2++)
		(fwProps.nuclearForces(i))(d1,d2).newSample( (fwNuclearForces(d1,d2))(i,whichIsDone)*norm, fWeight );
	  
	  resetFutureWalking(i,whichIsDone);
	} 
    }

  //if rel diff is really bad, then we have no reason to collect statistics
  if(rel_diff > Input->flags.rel_cutoff)
    return;
  
  if(globalInput.cs_Parameters.dim1() > 1)
    {	
      for(int cs=0; cs<cs_Energies.dim1(); cs++)
	{
	  fwProps.cs_Energies(cs).newSample( cs_Energies(cs), getWeight() * cs_Weights(cs));
	}
    }
  else
    {
      cs_Energies.deallocate();
      cs_Weights.deallocate();
    }

  // Calculate the basis function densities
  if (Input->flags.calculate_bf_density == 1)
    for (int i=0; i<Input->WF.getNumberBasisFunctions(); i++)
      fwProps.chiDensity(i).newSample( walkerData.chiDensity(i), getWeight() );
}

void QMCWalker::calculateDerivatives( QMCPropertyArrays & props )
{
  //if rel diff is really bad, then we have no reason to collect statistics
  double rel_diff = fabs( (localEnergy -
			   Input->flags.energy_estimated_original)/
			  Input->flags.energy_estimated_original);
  if(rel_diff > Input->flags.rel_cutoff)
    return;

  if(globalInput.flags.calculate_Derivatives == 1)
    {
      props.numDerHessSamples++;

      int numAI = rp_a.dim1();
      for(int ai=0; ai<numAI; ai++)
	{
	  /*
	    This isn't exactly the same since localEnergy is a weighted average
	    of original and trial walkers...
	  */
	  double e2 = localEnergy * localEnergy;
	  props.der(ai,0) += rp_a(ai);
	  props.der(ai,1) += rp_a(ai) * localEnergy;
	  props.der(ai,2) += rp_a(ai) * e2;
	  props.der(ai,3) += p3_xxa(ai);
	  props.der(ai,4) += p3_xxa(ai) * localEnergy;
	}
      
      if(globalInput.flags.optimize_Psi_criteria == "analytical_energy_variance")
	{
	  for(int ai=0; ai<numAI; ai++)
	    for(int aj=0; aj<=ai; aj++)
	      (props.hess(0))(ai,aj) += p3_xxa(ai) * p3_xxa(aj);
	  
	} else if(globalInput.flags.optimize_Psi_criteria == "generalized_eigenvector") {

	  for(int ai=0; ai<numAI; ai++)
	    for(int aj=0; aj<numAI; aj++)
	      {
		(props.hess(0))(ai,aj) += rp_a(ai) * rp_a(aj) * localEnergy;
		(props.hess(1))(ai,aj) += rp_a(ai) * rp_a(aj);
		(props.hess(2))(ai,aj) += rp_a(ai) * p3_xxa(aj); //hessian is not symmmetric 
	      }
	}
    } else {
      rp_a.deallocate();
      p3_xxa.deallocate();
    }
}

void QMCWalker::resetFutureWalking(int whichBlock, int whichIsDone)
{
  if(whichIsDone > 1 || whichIsDone < 0)
  {
    cerr << "Error: There are only 2 FW stages!\n";
    exit(1);
  }
  
  fwNormalization(whichBlock,whichIsDone)     = 0;
  fwR12(whichBlock,whichIsDone)               = 0;
  fwR2(whichBlock,whichIsDone)                = 0;
  fwiR12(whichBlock,whichIsDone)              = 0;
  fwiR(whichBlock,whichIsDone)                = 0;
  fwEnergy(whichBlock,whichIsDone)            = 0;
  fwKineticEnergy(whichBlock,whichIsDone)     = 0;
  fwKineticEnergy_grad(whichBlock,whichIsDone)     = 0;
  fwPotentialEnergy(whichBlock,whichIsDone)   = 0;
      
  if(Input->flags.nuclear_derivatives != "none")
    for (int d1=0; d1<fwNuclearForces.dim1(); d1++)
      for (int d2=0; d2<fwNuclearForces.dim2(); d2++)
        (fwNuclearForces(d1,d2))(whichBlock,whichIsDone) = 0;
}

void QMCWalker::resetFutureWalking()
{
  numFWSteps = 0;
  for(int i=0; i<isCollectingFWResults.dim1(); i++)
  {
    isCollectingFWResults(i,0) = ACCUM;
    isCollectingFWResults(i,1) = DONE;
    for(int j=0; j<isCollectingFWResults.dim2(); j++)
      resetFutureWalking(i,j);
  }
}


bool QMCWalker::isSingular()
{
  return walkerData.isSingular();
}

double QMCWalker::getLocalEnergyEstimator()
{
  return localEnergy;
}

