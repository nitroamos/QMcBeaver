/*
  Surfing the waves! wavefunctions, that is...
*/
#include "QMCSurfer.h"

using namespace std;

QMCSurfer::QMCSurfer()
{
  QMF = 0;
  walkerData.initialize(&globalInput,3,-1,-1);
}

QMCSurfer::~QMCSurfer()
{
  QMF = 0;
  R.deallocate();
}

int QMCSurfer::mainMenu(QMCFunctions * useQMF, int iteration,
			Array2D<double> & newR)
{
  /*
  int entry = min(globalInput.flags.max_time_steps -
		  globalInput.flags.equilibration_steps - 10,
		  (unsigned long int)1000);
  */
  R = newR;
  QMF = useQMF;
  static bool done = false;
  static int iteration_to_stop = iteration + 100 - 1;
  cout << "iteration_to_stop = " << iteration_to_stop << endl;
  if( iteration >= iteration_to_stop)
    {
      iteration_to_stop += 1000;
      done = true;
      cout << "We're at iteration " << iteration << ", next stop will be at " << iteration_to_stop << endl;

      char ok = 'e';
      while(ok != 'n' && ok != 'N' && ok != 'a' && ok != 'A')
	{
	  string mystr;
	  ok = 'n';
	  cout << "[e]xit [n]ext [s]can [p]ositions [a]nother walker [i]terations [" << ok << "]: ";
	  getline(cin,mystr);
	  if(mystr != "") ok = (mystr.c_str())[0];

	  if(ok == 'e' || ok == 'E')
	    {
	      exit(0);
	    } else if(ok == 'p' || ok == 'P')
	      {
		interparticleDistanceMatrix();
	      } else if(ok == 'n' || ok == 'N')
		{
		  //just so we don't fall to the error
		} else if(ok == 's' || ok == 'S')
		  {
		    surfaceExplorer();
		  } else if(ok == 'a' || ok == 'A')
		    {
		      iteration_to_stop -= 1000;
		      done = false;
		    } else if(ok == 'i' || ok == 'I')
		      {
			int skip = 10000;
			cout << "Iterations to skip [" << skip << "]: ";
			getline(cin,mystr);
			if(mystr != "") stringstream(mystr) >> skip;
			//subract the 1000 we already added
			iteration_to_stop -= 1000;
			iteration_to_stop += skip;
			cout << "Next stop will be at iteration " << iteration_to_stop << endl;
		      } else {
			cout << "Error: selected \'" << ok << "\'" << endl; 
		      }
	}
    } else {
      if( (iteration + globalInput.flags.equilibration_steps) % 1000 != 0)
	done = false;
    }
  return iteration_to_stop;
}

void QMCSurfer::interparticleDistanceMatrix()
{
  int numE = R.dim1();
  int numZ = globalInput.Molecule.getNumberAtoms();
  int numA = globalInput.WF.getNumberAlphaElectrons();
  double dist_hilight = 0.3;
  for(int i=0; i<numA; i++){
    for(int j=0; j<i; j++){
      double x = R(i,0) - R(j,0);
      double y = R(i,1) - R(j,1);
      double z = R(i,2) - R(j,2);
      double dist = sqrt(x*x + y*y + z*z);
      if( (i-numA+0.5)*(j-numA+0.5) > 0.0)
	{
	  cout << "\033[0;35m";
	} else {
	  cout << "\033[0;31m";
	}
      if(dist < dist_hilight) cout << "\033[0;36m";
      if(j%5 == 0) cout << endl;
      printf("(%2i,%2i %24.16f)  ",i,j,dist);
    }
  }
  for(int i=numA; i<numE; i++){
    for(int j=numA; j<i; j++){
      double x = R(i,0) - R(j,0);
      double y = R(i,1) - R(j,1);
      double z = R(i,2) - R(j,2);
      double dist = sqrt(x*x + y*y + z*z);
      if( (i-numA+0.5)*(j-numA+0.5) > 0.0)
	{
	  cout << "\033[0;35m";
	} else {
	  cout << "\033[0;31m";
	}
      if((j-numA)%5 == 0) cout << endl;
      if(dist < dist_hilight) cout << "\033[0;36m";
      printf("(%2i,%2i %24.16f)  ",i,j,dist);
    }
  }
  for(int i=0; i<numA; i++){
    for(int j=numA; j<numE; j++){
      double x = R(i,0) - R(j,0);
      double y = R(i,1) - R(j,1);
      double z = R(i,2) - R(j,2);
      double dist = sqrt(x*x + y*y + z*z);
      if( (i-numA+0.5)*(j-numA+0.5) > 0.0)
	{
	  cout << "\033[0;35m";
	} else {
	  cout << "\033[0;31m";
	}
      if((j-numA)%5 == 0) cout << endl;
      if(dist < dist_hilight) cout << "\033[0;36m";
      printf("(%2i,%2i %24.16f)  ",i,j,dist);
    }
  }
  cout << "\033[0m" << endl; 
  for(int i=0; i<numZ; i++){
    for(int j=0; j<numE; j++){
      double x = globalInput.Molecule.Atom_Positions(i,0) - R(j,0);
      double y = globalInput.Molecule.Atom_Positions(i,1) - R(j,1);
      double z = globalInput.Molecule.Atom_Positions(i,2) - R(j,2);
      double dist = sqrt(x*x + y*y + z*z);
      if(j%5 == 0) cout << endl;
      if(dist < dist_hilight)
	cout << "\033[0;36m";
      else
	cout << "\033[0m";
      printf("(%2s,%2i %24.16f)  ",globalInput.Molecule.Atom_Labels(i).data(),j,dist);
    }
  }
  cout << "\033[0m" << endl; 
}

void QMCSurfer::equipotentialSurface()
{
  int numE = R.dim1();
  int numZ = globalInput.Molecule.getNumberAtoms();

  for(int i=0; i<numE; i++){
    double potential = 0.0;
    for(int j=0; j<numE; j++){
      if(i == j) continue;
      double x = R(i,0) - R(j,0);
      double y = R(i,1) - R(j,1);
      double z = R(i,2) - R(j,2);
      double dist = sqrt(x*x + y*y + z*z);
      potential += 1.0/dist;
    }
    for(int j=0; j<numZ; j++){
      double x = globalInput.Molecule.Atom_Positions(j,0) - R(i,0);
      double y = globalInput.Molecule.Atom_Positions(j,1) - R(i,1);
      double z = globalInput.Molecule.Atom_Positions(j,2) - R(i,2);
      double dist = sqrt(x*x + y*y + z*z);
      potential -= globalInput.Molecule.Z(j) / dist;
    }
    printf("(%2i) %20.10f    ",i,potential);
    if((i+1)%5 == 0 || i == numE - 1) printf("\n");
  }
}

void QMCSurfer::scanEnergies(int moveE, int nucStart, int nucStop, int numSteps,
			     double startFrac, double stopFrac, double perturb)
{
  double sq2 = sqrt(2.0);
  double pi = acos(-1.0);
  double r=0, theta=0, phi=0;
  double x, y, z;
  double dx=-1, dy=-1, dz=-1;
  int numA = globalInput.Molecule.getNumberAtoms();
  //this ensures we don't actually land on the nucleus
  perturb = 1e-10;

  if(nucStart - numA > R.dim1()){
    cout << "Start index of " << nucStart << " is outside available range." << endl;
    return;
  }      
  if(nucStop - numA > R.dim1()){
    cout << "Stop index of " << nucStop << " is outside available range." << endl;
    return;
  }      

  if(nucStart < numA){
    x = globalInput.Molecule.Atom_Positions(nucStart,0);
    y = globalInput.Molecule.Atom_Positions(nucStart,1);
    z = globalInput.Molecule.Atom_Positions(nucStart,2);
  } else {
    x = R(nucStart-numA,0);
    y = R(nucStart-numA,1);
    z = R(nucStart-numA,2);
  }

  if(nucStart != nucStop && nucStop != -2){
    if(nucStart < numA){
      cout << "Moving electron " << moveE << " from atom " << nucStart << "\n " << globalInput.Molecule.Atom_Labels(nucStart);
    } else {
      cout << "Moving electron " << moveE << " from electron\ne" << (nucStart-numA);
    }
    cout << " at " << x << ", " << y << ", " << z;

    if(nucStop < numA){
      cout << "\nto atom " << nucStop << "\n " << globalInput.Molecule.Atom_Labels(nucStop) << " at ";
      dx = globalInput.Molecule.Atom_Positions(nucStop,0);
      dy = globalInput.Molecule.Atom_Positions(nucStop,1);
      dz = globalInput.Molecule.Atom_Positions(nucStop,2);
    } else {
      cout << "\nto electron\ne" << (nucStop-numA);
      dx = R(nucStop-numA,0);
      dy = R(nucStop-numA,1);
      dz = R(nucStop-numA,2);
    }
    cout << " at " << dx << ", " << dy << ", " << dz << endl;
  } else {

    if(nucStart < numA){
      cout << "Moving electron " << moveE << " relative to atom " << nucStart << "\n" << globalInput.Molecule.Atom_Labels(nucStart);
    } else {
      cout << "Moving electron " << moveE << " relative to electron\ne" << (nucStart-numA);
    }
    cout << " at " << x << ", " << y << ", " << z << endl;

    //This is if we want Electron 2 somewhere specific
    /*

    theta = 0.0;
    phi   = pi/2.0;
    r     = 2.0;
    
    R(1,0) = r*cos(theta)*sin(phi);
    R(1,1) = r*sin(theta)*sin(phi);
    R(1,2) = r*cos(phi);

    x = 0.0 - x;
    y = r - y;
    z = 0.0 - z;
    */
    //printf("Electron 2: x = %20.10e y = %20.10e z = %20.10e\n",R(1,0),R(1,1),R(1,2));
  }

  double delta;
  if(nucStop != -2)
    {
      dx = dx - x;
      dy = dy - y;
      dz = dz - z;
      
      //i'd rather have one more step actually printed than
      //have delta be a strange number (numsteps) or (numsteps-1)
      delta = (stopFrac - startFrac)/(numSteps);
      printf("scan distance = %10.5e (dx = %10.5f dy = %10.5f dz = %10.5f) delta = %10.5f\n",
	     sqrt(dx*dx+dy*dy+dz*dz),dx,dy,dz,delta);
      printf(" %10s %10s","step","dist","x","y","z");
    } else {
      r = startFrac;
      phi = stopFrac;

      delta = 2.0*pi/(numSteps);
      startFrac = 0.0;
      stopFrac  = 2.0*pi;
      printf("circumference = %10.5e, phi = %10.5f, dr = %10.5e\n",stopFrac*r,phi,delta);
      printf(" %10s %10s","theta","%","x","y","z");
    }

  printf(" %20s %20s %20s","TE","KE","PE");
  printf(" %20s %20s %20s %20s %20s","NE","EE","d(EE)/d(dist)","Virial","Psi");
  printf("\n");
  double ee_start = 0;
  for(double frac=startFrac; (frac - stopFrac) < delta; frac += delta)
    {
      double dist;
      if(nucStop != -2)
	{
	  //perturb doesn't move us off the vector
	  R(moveE,0) = x + (frac+perturb)*dx;
	  R(moveE,1) = y + (frac+perturb)*dy;
	  R(moveE,2) = z + (frac+perturb)*dz;

	  dist = (dx*dx + dy*dy + dz*dz);
	  dist = sqrt(dist)*(frac+perturb);
	} else {
	  theta = frac;

	  R(moveE,0) = x + r*cos(theta)*sin(phi);
	  R(moveE,1) = y + r*sin(theta)*sin(phi);
	  R(moveE,2) = z + r*cos(phi); 
	  
	  dist = r*frac;
	}

      /*
	 There is a minor bug in QMCBasisfunction such that
	 if an electron is at exactly 0 for any of it's coordinates,
	 it divides by zero. So we'll shift it very slightly.
      */
      for(int xyz=0; xyz<3; xyz++)
	if(R(moveE,xyz) == 0.0)
	  R(moveE,xyz) = perturb;

      /*      
      interparticleDistanceMatrix();
      for(int electron=0; electron<globalInput.WF.getNumberElectrons(); electron++)
	printf("el=%i: (%24.16e , %24.16e , %24.16e)\n",
	       electron,
	       R(electron,0),
	       R(electron,1),
	       R(electron,2));
	  //*/
      //equipotentialSurface();

      QMF->evaluate(R,walkerData);
      if((double)walkerData.psi > 0)
	{
	  //numbers 31 <= i <= 37 are permitted
	  cout << "\033[0;32m";
	} else {
	  cout << "\033[0;33m";
	}

      double sum = 0.0;
      int numE = globalInput.WF.getNumberElectrons();

      /*
      for(int i=0; i<numE; i++)
	{
	  double ke = -0.5*walkerData.laplacianE_PsiRatio(i);
	  double pe = walkerData.potentialE_PsiRatio(i);
	  if(i%4 == 0 && i > 0) cout << endl;
	  if(fabs(ke) > 1e3)
	    printf("(%i %15.8e ",i,ke);
	  else
	    printf("(%i %15.10f ",i,ke);

	  if(fabs(pe) > 1e3)
	    printf("%15.8e ",pe);
	  else
	    printf("%15.10f ",pe);

	  //this total does not include the atom-atom coulomb repulsion
	  if(fabs(pe+ke) > 1e3)
	    printf("%15.8e) ",pe+ke);
	  else
	    printf("%15.10f) ",pe+ke);
	  sum += ke + pe;
	}

      //cout << endl << "Sum = " << sum << endl;
      cout << endl;
      */      
      if(frac >= 0)
	{
	  printf(" %10.8f %10.5e",
		 frac,
		 dist);
	} else {
	  printf(" %10.7f %10.4e",
		 frac,
		 dist);
	}

      if(fabs(walkerData.localEnergy) >= 1e4)
	{
	  printf(" %20.10e",
		 walkerData.localEnergy);
	} else {
	  printf(" %20.14f",
		 walkerData.localEnergy);
	}
      if(fabs(walkerData.kineticEnergy) >= 1e4 ||
	 fabs(walkerData.potentialEnergy) >= 1e4)
	{
	  printf(" %20.10e %20.10e %20.10e",
		 walkerData.kineticEnergy,
		 walkerData.potentialEnergy,
		 walkerData.neEnergy);
	} else {
	  printf(" %20.14f %20.14f %20.14f",
		 walkerData.kineticEnergy,
		 walkerData.potentialEnergy,
		 walkerData.neEnergy);
	}

      if(ee_start == 0) ee_start = walkerData.eeEnergy;
      double deeddist = (walkerData.eeEnergy - ee_start)/(dist*delta);
      if(deeddist >= 0)
	{
	  //numbers 31 <= i <= 37 are permitted
	  /*
	    30 = black
	    31 = red
	    32 = dark green
	    33 = orange
	    34 = blue
	    35 = purple
	    36 = light blue
	    37 = normal terminal (green?)
	  */
	  cout << "\033[0;35m";
	} else {
	  cout << "\033[0;31m";
	}

      if(fabs(walkerData.eeEnergy) >= 1e4)
	{
	  printf(" %20.10e %20.10e",
		 walkerData.eeEnergy,
		 deeddist);
	} else {
	  printf(" %20.14f %20.10e",
		 walkerData.eeEnergy,
		 deeddist);
	}

      if((double)walkerData.psi > 0)
	{
	  //numbers 31 <= i <= 37 are permitted
	  cout << "\033[0;32m";
	} else {
	  cout << "\033[0;33m";
	}

      printf(" %20.10f %20.10e",
	     (-1.0*walkerData.potentialEnergy/walkerData.kineticEnergy),
	     (double)walkerData.psi);


      ee_start = walkerData.eeEnergy;
      printf("\n");
      cout << "\033[0m";
    }
  cout.flush();
}

void QMCSurfer::surfaceExplorer()
{
  int heaviest = 0;
  for(int i=0; i<globalInput.Molecule.getNumberAtoms(); i++)
    {
      if(globalInput.Molecule.Z(heaviest) < globalInput.Molecule.Z(i))
	heaviest = i;
    }
  static double startFrac = 0.0;
  static double stopFrac = 1.0;
  static int moveE = 0;
  static int startAtom = heaviest;
  static int stopAtom = -2;
  static int steps = 50;
  static double defR = 0.0000001;
  static double defP = acos(-1.0)/2.0;

  //char ok = 'y';
  char ok = 'n';

  while(ok == 'n' || ok == 'N')
    {
      string mystr;

      cout << "Enter electron to move [" << moveE << "]: ";
      getline(cin,mystr);
      if(mystr != "") stringstream(mystr) >> moveE;

      //cout << "Atom indices are 0 -> numAtoms-1; electron indicies start at numAtoms; -1 for all" << endl;
      cout << "Enter starting particle [" << startAtom << "]: ";
      getline(cin,mystr);
      if(mystr != "") stringstream(mystr) >> startAtom;
      
      //cout << "New option: -2 for ring" << endl;
      cout << "Enter stopping particle [" << stopAtom << "]: ";
      getline(cin,mystr);
      if(mystr != "") stringstream(mystr) >> stopAtom;
      
      cout << "Selected electron " << moveE << " to move from " << startAtom << " to " << stopAtom << endl; 
      if(startAtom <= -2 || moveE < 0 || moveE >= globalInput.WF.getNumberElectrons())
	{
	  cerr << "Bad values, try again\n";
	} else {
	  ok = 'y';
	  //cout << " ok? ";
	  //cin >> ok;
	}
    }
  ok = 'n';
  while(ok == 'n' || ok == 'N')
    {
      string mystr;
      cout << "Enter steps [" << steps << "]: ";
      getline(cin,mystr);
      if(mystr != "") stringstream(mystr) >> steps;

      if(stopAtom != -2)
	{ 
	  cout << "Enter starting fraction [" << startFrac << "]: ";
	  getline(cin,mystr);
	  if(mystr != "") stringstream(mystr) >> startFrac;
	  cout << "Enter stopping fraction [" << stopFrac << "]: ";
	  getline(cin,mystr);
	  if(mystr != "") stringstream(mystr) >> stopFrac;
	  cout << "Selected " << startFrac << " to " << stopFrac << " in " << steps << " steps" << endl; 
	} else {
	  startFrac = defR;
	  stopFrac = defP;
	  cout << "Enter radius [" << startFrac << "]: ";
	  getline(cin,mystr);
	  if(mystr != "") stringstream(mystr) >> startFrac;
	  cout << "Enter phi [" << stopFrac << "]: ";
	  getline(cin,mystr);
	  if(mystr != "") stringstream(mystr) >> stopFrac;
	  cout << "Selected radius " << startFrac << ", phi = " << stopFrac << " in " << steps << " steps" << endl; 
	}

      if(stopFrac == startFrac || (stopAtom == -2 && startFrac < 0.0))
	{
	  cerr << "Bad values, try again\n";
	} else {
	  ok = 'y';
	  if(stopAtom == -2)
	    {
	      defR = startFrac;
	      defP = stopFrac;
	    }
	  //cout << " ok? ";
	  //cin >> ok;
	}
    }

  for(int i=0; i<globalInput.Molecule.getNumberAtoms(); i++)
    {
      for(int j=0; j<globalInput.Molecule.getNumberAtoms(); j++)
	{
	  bool run = false;
	  if( startAtom == i && stopAtom == j)
	    run = true;
	  
	  if(startAtom >= globalInput.Molecule.getNumberAtoms() || stopAtom  >= globalInput.Molecule.getNumberAtoms())
	    {
	      if(i==0 && j==0){
		if(startAtom >= globalInput.Molecule.getNumberAtoms() &&
		   stopAtom  >= globalInput.Molecule.getNumberAtoms())
		  run = true;
		if( startAtom == -2 || stopAtom == -2)
		  run = true;		
	      }
	      if(i==j){
		if(startAtom == -1 || stopAtom == -1)
		  run = true;
		if(startAtom == i || stopAtom == j)
		  run = true;
	      }
	    }

	  if(startAtom == stopAtom && startAtom >= 0)
	    run = true;

	  if(startAtom == -1 && stopAtom == j && i != j)
	    run = true;

	  if(stopAtom  == -1 && startAtom == i && i != j)
	    run = true;

	  if(stopAtom  == -1 && startAtom == -1 && i < j)
	    run = true;

	  if(stopAtom  == -2 && startAtom == i && i == j)
	    run = true;

	  if(stopAtom  == -2 && startAtom == -1 && i == j)
	    run = true;

	  if(run){
	    int indexi = i, indexj = j;
	    if(stopAtom == -2) indexj = -2;
	    if(startAtom >= globalInput.Molecule.getNumberAtoms()) indexi = startAtom;
	    if(stopAtom >= globalInput.Molecule.getNumberAtoms()) indexj = stopAtom;
	    Array2D<double> temp = R;
	    scanEnergies(moveE,indexi,indexj,steps,startFrac,stopFrac,0);		
	    R = temp;
	    QMF->evaluate(R,walkerData);
	  }
	}
    }
}