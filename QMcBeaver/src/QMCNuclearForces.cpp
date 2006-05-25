/*
  Copyright (c) Amos G. Anderson 2006
  Distributed under GNU general public license (GPL)
  No guarantee or warantee regarding usability or stability is expressed or implied.
  nitroamos@gmail.com
*/

#include "QMCNuclearForces.h"

QMCNuclearForces::QMCNuclearForces()
{}

QMCNuclearForces::~QMCNuclearForces()
{
  waveValuesHFSpline.deallocate();
  numPolyBasisF.deallocate();
  nucleusContributions.deallocate();
  alphaOrbitals.deallocate();
  betaOrbitals.deallocate();
  
  for(int i=0; i<radialPoints.dim1(); i++)
    radialPoints(i).deallocate();
  radialPoints.deallocate();

  for(int i=0; i<waveValuesHF.dim1(); i++)
  for(int j=0; j<waveValuesHF.dim2(); j++)
    waveValuesHF(i,j).deallocate();
  waveValuesHF.deallocate();

  for(int i=0; i<basisCoeffs.dim1(); i++)
  for(int j=0; j<basisCoeffs.dim2(); j++)
    basisCoeffs(i,j).deallocate();
  basisCoeffs.deallocate();
}

void QMCNuclearForces::initialize(QMCInput *Ip)
{
  if(Ip->flags.nuclear_derivatives == "none")
    return;
    
  Input = Ip;
  BF = & Input->BF;
  WF = & Input->WF;
  
  int constNumBF = 5;
  fittingWeightM = 0.0;
  int numNuc = Input->Molecule.getNumberAtoms();
  numSamplesPerArea = 1e2;
  
  numPolyBasisF.allocate(numNuc);
  basisCoeffs.allocate(numNuc,3);
  radialPoints.allocate(numNuc);
  waveValuesHF.allocate(numNuc,4);//1 s + 3 p
  waveValuesHFSpline.allocate(numNuc,3);

  for(int i=0; i<numNuc; i++)
  {
    numPolyBasisF(i) = constNumBF;
    for(int q=0; q<3; q++)
      basisCoeffs(i,q).allocate(numPolyBasisF(i));
    radialPoints(i).allocate(1);
    (radialPoints(i))(0) = 1.0;
  }

  calculateNuclearContributions();

  for(int nuc=0; nuc<numNuc; nuc++)
  {
    calcCoefficients(nuc);
    
    if(false){
      cout << "\n Coeff. " << Input->Molecule.Atom_Labels(nuc) << ": ";
      for(int j=0; j<basisCoeffs(nuc,0).dim1(); j++)
        printf(" %15.7e",(basisCoeffs(nuc,0))(j));
      
      cout << endl;
      printf("%15s   %15s %15s\n","rAi","functionValue","temperTerm");
      for(int rIndex = 0; rIndex<radialPoints(nuc).dim1(); rIndex++)
      {
        int numBF = numPolyBasisF(nuc);
        int q = 0;
        double rAi = (radialPoints(nuc))(rIndex);
        double temperTerm = (basisCoeffs(nuc,0))(numBF-1);
        for(int i=numBF-2; i >= 0; i--)
        {
          temperTerm = temperTerm*rAi + (basisCoeffs(nuc,q))(i);
        }
        
        if(Input->flags.nuclear_derivatives == "slaterpoly_hellmann_feynman")
        {
          spliner.initializeWithFunctionValues(radialPoints(nuc),
                                               waveValuesHF(nuc,q+1), 0.0, 0.0);
          spliner.evaluate(rAi);
          temperTerm *= spliner.getFunctionValue()*pow(rAi,fittingWeightM);                
        } else {
          temperTerm *= pow(rAi,1.0 + fittingWeightM);
        }
        
        printf("%15.7e   %15.7e %15.7e\n",rAi,spliner.getFunctionValue(),temperTerm);
      }
    }
    //exit(0);
  }
  
  for(int nuc=0; nuc<1; nuc++){
    if(false){
      printf("%15s %15s %15s %15s %15s\n","rad","s","px","py","pz");
      for(int j=0; j<radialPoints(nuc).dim1(); j++)
        printf("%15.7f %15.7e %15.7e %15.7e %15.7e\n",
          (radialPoints(nuc))(j),
          (waveValuesHF(nuc,0))(j),
          (waveValuesHF(nuc,1))(j),
          (waveValuesHF(nuc,2))(j),
          (waveValuesHF(nuc,3))(j));
    }
  }
  //double x = ran1(&Input->flags.iseed);
  //exit(0);
}

void QMCNuclearForces::evaluate(QMCWalkerData &walkerData, Array2D<double> &xData)
{
  Array1D<QMCWalkerData *> oneWalker = Array1D<QMCWalkerData *>(1);
  Array1D<Array2D<double> * > oneX = Array1D<Array2D<double> * >(1);
  oneWalker(0) = & walkerData;
  oneX(0) = & xData;
  evaluate(oneWalker,oneX,1);
  oneWalker.deallocate();
  oneX.deallocate();
}

/*
  This is the function that would be repeatedly called each QMC
  iteration
*/
void QMCNuclearForces::evaluate(Array1D<QMCWalkerData *> &walkerData, 
		Array1D<Array2D<double> * > &xData, int num)
{
  int numNuc = Input->Molecule.getNumberAtoms();
  int numEle = xData(0)->dim1();
  
  /*
    If we are doing a binning procedure, then we don't
    need the nucleus-nucleus contributions to the force.
    
    The same data structure is used for both force and bin output.
  */
  for(int walker=0; walker<num; walker++)
    if(Input->flags.nuclear_derivatives == "bin_force_density"){
      walkerData(walker)->nuclearDerivatives = 0;
    } else {
      walkerData(walker)->nuclearDerivatives = nucleusContributions;    
    }
    
  for(int walker=0; walker<num; walker++)
  {
    for(int nucA=0; nucA<numNuc; nucA++)
    {
      for(int eleI=0; eleI<numEle; eleI++)
      {
        
        double rAi = 0;
        for(int q=0; q<3; q++)
        {
          rAi +=
          ( (*xData(walker))(eleI,q) -
            Input->Molecule.Atom_Positions(nucA,q) ) *
          ( (*xData(walker))(eleI,q) -
            Input->Molecule.Atom_Positions(nucA,q) );
        }
        rAi = sqrt(rAi);
        
        double nuclearForce = 0;
        for(int q=0; q<3; q++)
        {
          nuclearForce = (*xData(walker))(eleI,q) -
                         Input->Molecule.Atom_Positions(nucA,q);
          nuclearForce *= Input->Molecule.Z(nucA);
          nuclearForce *= -1.0 / (rAi * rAi * rAi);
          
          /*
            If we are below some cutoff, then we probably want
            to multiply the nuclear force by some tempering term
            which is designed to cut the s-wave contributions out,
            leaving only the p-wave.
            This tempering term is \frac{ \rho_z(r) * z }{ \rho(r) * r }
            which needs to go to zero as r goes to zero.
          */
          int numBF = basisCoeffs(nucA,q).dim1();

	  double radialCutoff = (radialPoints(nucA))(radialPoints(nucA).dim1()-1);
          if(rAi < radialCutoff)
          {
            if(Input->flags.nuclear_derivatives == "bare_hellmann_feynman")
            {
              //do nothing
            }
            
            if(Input->flags.nuclear_derivatives == "poly_hellmann_feynman" ||
               Input->flags.nuclear_derivatives == "slaterpoly_hellmann_feynman")
            {

	      //This evaluates the polynomial (Eq. 7 from CCZ05)
	      //where the m contributions have been factored out
              double temperTerm = (basisCoeffs(nucA,q))(numBF-1);
              for(int i=numBF-2; i >= 0; i--)
              {
                temperTerm = temperTerm*rAi + (basisCoeffs(nucA,q))(i);
              }
              
              if(Input->flags.nuclear_derivatives == "slaterpoly_hellmann_feynman")
              {
		//This adds in the Slater term (Eq. 8 from CCZ05)
                waveValuesHFSpline(nucA,q).evaluate(rAi);
                temperTerm *= waveValuesHFSpline(nucA,q).getFunctionValue()
		  *pow(rAi,fittingWeightM);                
              } else {
		//Now we add the m contribution to Eq. 7
                temperTerm *= pow(rAi,1.0 + fittingWeightM);
              }

	      //This completes Eq. 6
              nuclearForce *= temperTerm;
            }

          }//end if
          
          if(Input->flags.nuclear_derivatives == "bin_force_density")
          {
            //let's only bin the 1st nucleus' z-axis
	    if(q == 2 && nucA == 0)
	      {
		int binIndex = (int)( min(rAi/radialCutoff,1.0)*( getNumBins() - 1) );
		walkerData(walker)->nuclearDerivatives(binIndex,1) += 1;
	      }
          } else {
            walkerData(walker)->nuclearDerivatives(nucA,q) += nuclearForce;
          }
        }//end q
      }//end eleI
    }//end nucA
  }//end walker
  
  if(false)
    {
      cout << "\ncalculated force:\n";
      for(int w=0; w<num; w++){
	for(int i=0; i<4; i++)
	  printf("%s: %17.7f %17.7f %17.7f\n",
		 Input->Molecule.Atom_Labels(i).c_str(),
		 walkerData(w)->nuclearDerivatives(i,0),
		 walkerData(w)->nuclearDerivatives(i,1),
		 walkerData(w)->nuclearDerivatives(i,2));
	cout << endl;
      }
    }
  //exit(0);
}

void QMCNuclearForces::calcCoefficients(int whichNucleus)
{
  if(Input->flags.nuclear_derivatives == "none" ||
     Input->flags.nuclear_derivatives == "bare_hellmann_feynman" ||
     Input->flags.nuclear_derivatives == "bin_force_density")
  {
    //nothing to be done
    return;
  }

  double radialCutoff = (radialPoints(whichNucleus))(radialPoints(whichNucleus).dim1()-1);
  int numBF = numPolyBasisF(whichNucleus);
  Array2D<double> hilbertMatrix  = Array2D<double> (numBF,numBF);
  Array2D<double> inverseMatrix  = Array2D<double> (numBF,numBF);
  Array2D<double> residualVector = Array2D<double> (numBF, 1);
  Array2D<double> coeffs         = Array2D<double> (numBF, 1);
  double det=0; bool isOK=true;
  
  if(Input->flags.nuclear_derivatives == "poly_hellmann_feynman")
  {
    /*
      This section calculates the coefficients as detailed in formula 8
      of the paper. The coefficients are the same for all three coordinate
      directions.
      
      The coefficients calculated (formula 8) are not the same as those shown
      in formula 4. These are the result of some manipulations.
       
      Note: Fortran and C indices are different. This code was indexed to
      look like the Fortran loops.
    */
    for(int j=1; j<=numBF; j++)
    {
      /*
        calculate the residual vector
        basis functions are: \phi_i(r) = r^i (where i starts indexing at 1)
        Integral is Rj = \int_0^R dr \phi_j(r)
          = \frac{ R^{j + 1} }{j + 1}
      */
      residualVector(j-1,0) = pow(radialCutoff, j + 1.0) / ( j + 1.0 );
      for(int k=1; k<=numBF; k++)
      {
        /*
          calculate the Hilbert matrix
          Integral is Sij = \int_0^R dr r^m \phi_i(r) \phi_j(r)
            = \frac{ R^{ m + k + j + 1} }{ m + k + j + 1 }
        */
        hilbertMatrix(j-1,k-1) = pow(radialCutoff, fittingWeightM + k + j + 1.0) /
                                 ( fittingWeightM + k + j + 1.0 );
      }
    }
    
    //determinant_and_inverse(hilbertMatrix,inverseMatrix,det,&isOK);
    hilbertMatrix.determinant_and_inverse(inverseMatrix,det,&isOK);

    if(!isOK)
    {
      cerr << "ERROR: Hilbert matrix is not ok\n";
      exit(0);
    }
    
    // \bm{c} = \bm{ S^{-1} } \bm{ R }
    //coeffs = inverseMatrix.multiply(residualVector);

    inverseMatrix.gemm(residualVector,coeffs,false);

    for(int i=0; i<numBF; i++)
      for(int coord=0; coord<3; coord++)
        (basisCoeffs(whichNucleus,coord))(i) = coeffs(i,0);
    
  }
  else if(Input->flags.nuclear_derivatives == "slaterpoly_hellmann_feynman")
  {
    /*
      This section calculates the coefficients for formula 9
      of the paper.
      
      These formulas have also been manipulated, but we can use
      the same procedure detailed in the supplement, the formulas
      are shown in the comments.
       
      Note: I went ahead with C like indices here
    */
    int numPT = 100;
    waveMemorization(whichNucleus, numPT, radialCutoff);

    Array1D<double> radialValues = Array1D<double>(numPT);

    for(int coord=0; coord<3; coord++)
    {
      /*
        calculate the residual vector
        basis functions are: \phi_i(r) = f_z(r) r^i (where i starts indexing at 0)
        Integral is Ri = \int_0^R dr \phi_i(r)
      */
      for(int powerI=0; powerI<numBF; powerI++)
      {
        radialValues = waveValuesHF(whichNucleus,coord+1);
        for(int rIndex=0; rIndex<numPT; rIndex++)
        {
          radialValues(rIndex) *= fastPower((radialPoints(whichNucleus))(rIndex),powerI);
        }

        //it seems reasonable to choose 0.0 for both end point derivatives
        spliner.initializeWithFunctionValues(radialPoints(whichNucleus), radialValues, 0.0, 0.0);
        residualVector(powerI,0) = spliner.integrate(0, numPT);
      }
      
      /*
        calculate the Hilbert matrix
        Integral is Sij = \int_0^R dr r^m \phi_i(r) \phi_j(r) 
      */
      for(int powerI=0; powerI<numBF; powerI++)
      {
        for(int powerJ=0; powerJ<=powerI; powerJ++)
        {
          radialValues = waveValuesHF(whichNucleus,coord+1);
          for(int rIndex=0; rIndex<numPT; rIndex++)
          {
            radialValues(rIndex) *= radialValues(rIndex) *
                                    pow((radialPoints(whichNucleus))(rIndex),
                                    powerI + powerJ + fittingWeightM);
          }
         
          spliner.initializeWithFunctionValues(radialPoints(whichNucleus), radialValues, 0.0, 0.0);
          hilbertMatrix(powerI,powerJ) = spliner.integrate(0,numPT);
          /*
          //This code will calculate the fit to f_z(r)
          hilbertMatrix(powerI,powerJ) = pow(rMaxNuc(whichNucleus), fittingWeightM + powerI + powerJ + 1.0) /
                                   ( fittingWeightM + powerI + powerJ + 1.0 );
          */
          
          //the Hilbert matrix is symmetric
          if(powerI != powerJ) hilbertMatrix(powerJ,powerI) = hilbertMatrix(powerI,powerJ);
        }
      }

      hilbertMatrix.determinant_and_inverse(inverseMatrix,det,&isOK);
      //determinant_and_inverse(hilbertMatrix,inverseMatrix,det,&isOK);
      
      if(!isOK)
      {
        cerr << "ERROR: Hilbert matrix was not inverted.\n";
        exit(1);
      }

      inverseMatrix.gemm(residualVector,coeffs,false); 

      for(int i=0; i<numBF; i++)
        (basisCoeffs(whichNucleus,coord))(i) = coeffs(i,0);

    }// end for coord
  }
  else
  {
    cerr << "ERROR: \"" << Input->flags.nuclear_derivatives << 
      "\" is an unknown nuclear force method.\n";
    cerr << "The choices are:\n" <<
      "\tnone\n" <<
      "\tbin_force_density\n" <<      
      "\tbare_hellmann_feynman\n" <<
      "\tpoly_hellmann_feynman\n" <<
      "\tslaterpoly_hellmann_feynman\n" << endl;
    exit(1);
  }

  coeffs.deallocate();  
  hilbertMatrix.deallocate();
  inverseMatrix.deallocate();
  residualVector.deallocate();
}

void QMCNuclearForces::waveMemorization(int whichNucleus, int numKnots, double radialCutoff)
{
  if(numKnots <= 1)
    {
      cout << "Error: we need more than 1 shell to memorize the wavefunction!\n";
      exit(1);
    }
  double pi = acos(-1.0);
  double dr = radialCutoff/(numKnots-1);
  
  //we will reallocate radialPoints and the original value will now be the last value
  radialPoints(whichNucleus).allocate(numKnots);
  for(int i=0; i<4; i++)
    {
      waveValuesHF(whichNucleus,i).allocate(numKnots);
      waveValuesHF(whichNucleus,i) = 0;
    }
  
  Array1D<double> center = Array1D<double>(3);
  for(int i=0; i<3; i++)
    center(i) = Input->Molecule.Atom_Positions(whichNucleus,i);
  
  Array2D<double> cube = Array2D<double>(8,3);
  Array1D<double> densities = Array1D<double>(8);
  for(int irad=0; irad<numKnots; irad++)
    {
      double rad = irad*dr;
      (radialPoints(whichNucleus))(irad) = rad;
      int numSamples = (int)(rad*rad*numSamplesPerArea + 1.0);
      //printf("%5i ",numSamples);
      double aveSurfacePerPoint = 8.0*numSamples/(4.0*pi);
      
      for(int sample = 0; sample < numSamples; sample++)
	{
	  generateCube(cube,rad);
	  randomlyRotate(cube,1.0);
	  
	  for(int i=0; i<cube.dim1(); i++)
	    for(int coord=0; coord<3; coord++)
	      cube(i,coord) += center(coord);
	  
	  getDensities(cube,densities);
	  
	  /*
	    For spherical wavefunctions, we will get 8 times the contribution.
	    For p-waves, the 8 pieces will cancel each other out to the degree
	    that the wave is 3-fold symmetric.
	  */
	  for(int i=0; i<densities.dim1(); i++)
	    { 
	      (waveValuesHF(whichNucleus,0))(irad) += densities(i);          
	      for(int j=0; j<3; j++)
		(waveValuesHF(whichNucleus,j+1))(irad) += densities(i)*(cube(i,j)-center(j))/rad;
	    }
	  
	  if(irad == 0)
	    for(int j=0; j<3; j++)
	      (waveValuesHF(whichNucleus,j+1))(irad) = 0;
	}//samples
      
      for(int i=0; i<4; i++)
	(waveValuesHF(whichNucleus,i))(irad) /= aveSurfacePerPoint;
      
    }//distances

  //This is the spline we'll use to preserve our measurements.
  for(int q=0; q<3; q++)
    waveValuesHFSpline(whichNucleus,q).initializeWithFunctionValues(radialPoints(whichNucleus),
								    waveValuesHF(whichNucleus,q+1),0.0,0.0);
  
  center.deallocate();
  cube.deallocate();
  densities.deallocate();
}

void QMCNuclearForces::getDensities(Array2D<double> & X, Array1D<double> & densities)
{
  int nBasisFun = WF->getNumberBasisFunctions();
  int nOrbitals = WF->getNumberOrbitals();
  Array2D<double> Chi = Array2D<double>(X.dim1(),nBasisFun);
  alphaOrbitals.allocate(X.dim1(),nOrbitals);
  betaOrbitals.allocate(X.dim1(),nOrbitals);  
  
  BF->evaluateBasisFunctions(X,Chi);
  //alphaOrbitals = Chi * Input->WF.AlphaCoeffs;
  Chi.gemm(Input->WF.AlphaCoeffs,alphaOrbitals,true);
  if (Input->flags.trial_function_type == "restricted")
  {
    betaOrbitals = alphaOrbitals;
  }
  else if (Input->flags.trial_function_type == "unrestricted")
  {
    //betaOrbitals = Chi * Input->WF.BetaCoeffs;
    Chi.gemm(Input->WF.AlphaCoeffs,betaOrbitals,true);
  }

  densities = 0;
  for(int i=0; i<X.dim1(); i++)
  {
    for(int j=0; j<WF->getNumberDeterminants(); j++)
    {
      double detSum = 0.0;
      for(int k=0; k<nOrbitals; k++)
      {
        if(Input->WF.AlphaOccupation(j,k) == 1)
          detSum += alphaOrbitals(i,k)*alphaOrbitals(i,k);
        if(Input->WF.BetaOccupation(j,k) == 1)
          detSum += betaOrbitals(i,k)*betaOrbitals(i,k);
      }
      densities(i) += detSum*WF->CI_coeffs(j);
    }
  }
  Chi.deallocate();
}

void QMCNuclearForces::printPoints(Array2D<double> & points)
{
  for(int i=0; i<points.dim1(); i++)
  {
    printf("%10.7f  %10.7f  %10.7f\n",points(i,0),points(i,1),points(i,2));
  }
  cout << endl;
}

void QMCNuclearForces::printCubeLengths(Array2D<double> & points)
{
  double dist;
  for(int i=0; i<points.dim1(); i++)
  {
    dist = sqrt( points(i,0)*points(i,0)
               + points(i,1)*points(i,1)
               + points(i,2)*points(i,2) );
    printf("ORG %3i %10.7f\n",i,dist);
      
    for(int j=i+1; j<points.dim1(); j++)
    {
      dist = sqrt( (points(i,0)-points(j,0))*(points(i,0)-points(j,0))
                 + (points(i,1)-points(j,1))*(points(i,1)-points(j,1))
                 + (points(i,2)-points(j,2))*(points(i,2)-points(j,2)) );
      printf("%3i %3i %10.7f\n",i,j,dist);
    }
  }
  cout << endl;
}

void QMCNuclearForces::generateCube(Array2D<double> & cube, double length)
{
  //8 vertices, 3 coordinates per vertex
  cube.allocate(8,3);
  double pos = length/sqrt(3.0);
  cube(0,0) =  pos;  cube(0,1) =  pos;  cube(0,2) =  pos;
  cube(1,0) =  pos;  cube(1,1) =  pos;  cube(1,2) = -pos;
  cube(2,0) =  pos;  cube(2,1) = -pos;  cube(2,2) =  pos;
  cube(3,0) =  pos;  cube(3,1) = -pos;  cube(3,2) = -pos;
  cube(4,0) = -pos;  cube(4,1) =  pos;  cube(4,2) =  pos;
  cube(5,0) = -pos;  cube(5,1) =  pos;  cube(5,2) = -pos;
  cube(6,0) = -pos;  cube(6,1) = -pos;  cube(6,2) =  pos;
  cube(7,0) = -pos;  cube(7,1) = -pos;  cube(7,2) = -pos;
}

void QMCNuclearForces::randomlyRotate(Array2D<double> & points, double scale)
{
  //this is the so-called "y-convention" or "zyz" rotation
  Array2D<double> rotmat = Array2D<double>(3,3);
  double pi = acos(-1.0);
  //phi (aka alpha) goes from 0 to 2pi
  double phi = 2.0*pi*ran1(&Input->flags.iseed);
  //theta (aka beta) ranges between 0 to pi, but different distribution
  double theta = acos(2.0*(ran1(&Input->flags.iseed) - 0.5));
  //psi (aka gamma) goes from 0 to 2pi
  double psi = 2.0*pi*ran1(&Input->flags.iseed);
  
  double sinpsi, cospsi, sinthe, costhe, sinphi, cosphi;
  sinphi = sin(phi);   cosphi = cos(phi);
  sinpsi = sin(psi);   cospsi = cos(psi);
  sinthe = sin(theta); costhe = cos(theta);

  //(row, column)
  rotmat(0,0) =  cosphi*costhe*cospsi - sinphi*sinpsi;
  rotmat(0,1) =  sinphi*costhe*cospsi + cosphi*sinpsi;
  rotmat(0,2) = -sinthe*cospsi;
  rotmat(1,0) = -cosphi*costhe*sinpsi - sinphi*cospsi;
  rotmat(1,1) = -sinphi*costhe*sinpsi + cosphi*cospsi;
  rotmat(1,2) =  sinthe*sinpsi;
  rotmat(2,0) =  cosphi*sinthe;
  rotmat(2,1) =  sinphi*sinthe;
  rotmat(2,2) =  costhe;

  double x, y, z;
  for(int i=0; i<points.dim1(); i++)
  {
    x = points(i,0);
    y = points(i,1);
    z = points(i,2);
    for(int j=0; j<3; j++)
    {
      points(i,j) = x*rotmat(j,0)
                  + y*rotmat(j,1)
                  + z*rotmat(j,2);
      points(i,j) *= scale;
    }
  }
  rotmat.deallocate();
}

void QMCNuclearForces::calculateNuclearContributions()
{
  int numNuc = Input->Molecule.getNumberAtoms();
  nucleusContributions.allocate(numNuc,3);
  nucleusContributions = 0;

  for(int nucA=0; nucA<numNuc; nucA++)
    {
      for(int nucB=0; nucB<nucA; nucB++)
	{ 
	  double rAB = 0;
	  for(int q=0; q<3; q++)
	    {
	      rAB +=
		(Input->Molecule.Atom_Positions(nucA,q) -
		 Input->Molecule.Atom_Positions(nucB,q) ) *
		(Input->Molecule.Atom_Positions(nucA,q) -
		 Input->Molecule.Atom_Positions(nucB,q) );
	    }
	  rAB = sqrt(rAB);
	  
	  double nuclearForce = 0;
	  for(int q=0; q<3; q++)
	    {
	      nuclearForce =
		Input->Molecule.Atom_Positions(nucB,q) -
		Input->Molecule.Atom_Positions(nucA,q);
	      nuclearForce *= 
		Input->Molecule.Z(nucA) *
		Input->Molecule.Z(nucB);
	      nuclearForce *=
		1.0 / (rAB * rAB * rAB);
	      
	      nucleusContributions(nucA,q) += nuclearForce;
	    }//end q  
	}//end nucB
    }//end nucA
}
  



