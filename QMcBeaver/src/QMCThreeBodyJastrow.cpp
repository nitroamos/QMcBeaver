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

#include "QMCThreeBodyJastrow.h"
#include "MathFunctions.h"

QMCThreeBodyJastrow::QMCThreeBodyJastrow()
{
}

QMCThreeBodyJastrow::~QMCThreeBodyJastrow()
{
}

void QMCThreeBodyJastrow::initialize(QMCInput * input)
{
}

double QMCThreeBodyJastrow::get_p3_xxa_ln(int ai)
{
  return p3_xxa(ai);
}

Array2D<double> * QMCThreeBodyJastrow::get_p2_xa_ln(int ai)
{
  return &p2_xa(ai);
}

double QMCThreeBodyJastrow::get_p_a_ln(int ai)
{
  return p_a(ai);
}

void QMCThreeBodyJastrow::evaluate(QMCJastrowParameters & JP, 
				   QMCWalkerData * wData, Array2D<double> & X)
{
  wd = wData;
  if(wd->whichE == -1)
    updateAll(JP,X);
  else
    updateOne(JP,X);
  wd = 0;
}

void QMCThreeBodyJastrow::updateOne(QMCJastrowParameters & JP, Array2D<double> & X)
{
  int nalpha = globalInput.WF.getNumberElectrons(true);
  int nbeta = globalInput.WF.getNumberElectrons(false);
  int E = wd->whichE;
  bool isAlpha = E < nalpha ? true : false;

  for(int E2=0; E2<wd->UijA.dim1(); E2++)
    {
      if(E == E2) continue;
      int EH = max(E,E2);
      int EL = min(E,E2);

      wd->U -= wd->UijA(EH,EL);
      wd->UijA(EH,EL) = 0.0;
      
      wd->U_xx -= wd->UijA_xx(EH,EL);
      wd->UijA_xx(EH,EL) = 0.0;

      for(int xyz=0; xyz<3; xyz++)
	{
	  wd->U_x(EH,xyz)  -= (wd->UijA_x1(xyz))(EH,EL); 
	  (wd->UijA_x1(xyz))(EH,EL) = 0.0;
	  wd->U_x(EL,xyz) -= (wd->UijA_x2(xyz))(EH,EL); 
	  (wd->UijA_x2(xyz))(EH,EL) = 0.0;
	}
    }

  Array1D<string> * NucleiTypes = JP.getNucleiTypes();

  EupEdnNuclear = JP.getElectronUpElectronDownNuclearParameters();
  for(int nuc=0; nuc<EupEdnNuclear->dim1(); nuc++)
    (*EupEdnNuclear)(nuc).zeroOutDerivatives();

  EupEupNuclear = JP.getElectronUpElectronUpNuclearParameters();
  for(int nuc=0; nuc<EupEupNuclear->dim1(); nuc++)
    (*EupEupNuclear)(nuc).zeroOutDerivatives();
  
  EdnEdnNuclear = JP.getElectronDownElectronDownNuclearParameters();
  for(int nuc=0; nuc<EdnEdnNuclear->dim1(); nuc++)
    (*EdnEdnNuclear)(nuc).zeroOutDerivatives();

  Array1D<double> xyz(3);
  for(int Nuclei=0; Nuclei<globalInput.Molecule.getNumberAtoms(); Nuclei++)
    {
      // Find the number of the current nucleus in the nuclei list
      int NucleiType = -1;
      for( int i=0; i<NucleiTypes->dim1(); i++ )
	if( globalInput.Molecule.Atom_Labels(Nuclei) == (*NucleiTypes)(i) )
	  {
	    NucleiType = i;
	    break;
	  } 
      bool used;
      
      // Now we do all opposite spin pairs with this nucleus	 
      if (nalpha > 0 && nbeta > 0)
	if(isAlpha)
	  {
	    for (int EB=nalpha; EB<X.dim1(); EB++)
	      {
		for(int i=0; i<3; i++) xyz(i)  = wd->riA_uvec(EB,Nuclei,i);
		used = (*EupEdnNuclear)(NucleiType).getThreeBodyCorrelationFunction()->setElectron(true, xyz, wd->riA(EB,Nuclei));
		if(!used) continue;		
		collectForPair(EB, E, Nuclei,& (*EupEdnNuclear)(NucleiType));
	      }
	  } else {
	    for (int EA=0; EA<nalpha; EA++)
	      {
		for(int i=0; i<3; i++) xyz(i)  = wd->riA_uvec(E,Nuclei,i);
		used = (*EupEdnNuclear)(NucleiType).getThreeBodyCorrelationFunction()->setElectron(true, xyz, wd->riA(E,Nuclei));
		if(!used) continue;
		collectForPair(E, EA, Nuclei, & (*EupEdnNuclear)(NucleiType));
	      }
	  }

      // First we do all pairs of two alphas with this nucleus

      if (nalpha > 1)
	if(isAlpha)
	  {
	    for(int EA=0; EA<E; EA++)
	      {
		for(int i=0; i<3; i++) xyz(i)  = wd->riA_uvec(E,Nuclei,i);
		used = (*EupEupNuclear)(NucleiType).getThreeBodyCorrelationFunction()->setElectron(true, xyz, wd->riA(E,Nuclei));
		if(!used) continue;
		collectForPair(E, EA, Nuclei, & (*EupEupNuclear)(NucleiType));
	      }
	    for(int EA=E+1; EA<nalpha; EA++)
	      {
		for(int i=0; i<3; i++) xyz(i)  = wd->riA_uvec(EA,Nuclei,i);
		used = (*EupEupNuclear)(NucleiType).getThreeBodyCorrelationFunction()->setElectron(true, xyz, wd->riA(EA,Nuclei));
		if(!used) continue;		
		collectForPair(EA, E, Nuclei, & (*EupEupNuclear)(NucleiType));
	      }
	  }
      
      // Finally we do all beta spin pairs with this nucleus
      if (nbeta > 1)
	if(!isAlpha)
	  {
	    for(int EB=nalpha; EB<E; EB++)
	      {
		for(int i=0; i<3; i++) xyz(i)  = wd->riA_uvec(E,Nuclei,i);
		used = (*EdnEdnNuclear)(NucleiType).getThreeBodyCorrelationFunction()->setElectron(true, xyz, wd->riA(E,Nuclei));
		if(!used) continue;
		collectForPair(E, EB, Nuclei, & (*EdnEdnNuclear)(NucleiType));
	      }
	    for(int EB=E+1; EB<X.dim1(); EB++)
	      {
		for(int i=0; i<3; i++) xyz(i)  = wd->riA_uvec(EB,Nuclei,i);
		used = (*EdnEdnNuclear)(NucleiType).getThreeBodyCorrelationFunction()->setElectron(true, xyz, wd->riA(EB,Nuclei));
		if(!used) continue;
		collectForPair(EB, E, Nuclei, & (*EdnEdnNuclear)(NucleiType));
	      }
	  }
    }

  for(int E2=0; E2<wd->UijA.dim1(); E2++)
    {
      if(E == E2) continue;
      int EH = max(E,E2);
      int EL = min(E,E2);

      wd->U    += wd->UijA(EH,EL);
      wd->U_xx += wd->UijA_xx(EH,EL);
      for(int xyz=0; xyz<3; xyz++)
	{
	  wd->U_x(EH,xyz) += (wd->UijA_x1(xyz))(EH,EL); 
	  wd->U_x(EL,xyz) += (wd->UijA_x2(xyz))(EH,EL); 
	}
    }  

  xyz.deallocate();
}

double QMCThreeBodyJastrow::jastrowOnGrid(QMCJastrowParameters & JP,
					  int E,
					  Array2D<double> & R,
					  Array2D<double> & grid,
					  Array1D<double> & integrand)
{
  int nalpha = globalInput.WF.getNumberElectrons(true);
  int nbeta = globalInput.WF.getNumberElectrons(false);
  bool isAlpha = E < nalpha ? true : false;

  double denom = 0.0;
  Array1D<double> sumU(integrand.dim1());
  sumU = 0.0;

  Array1D<string> * NucleiTypes = JP.getNucleiTypes();

  EupEdnNuclear = JP.getElectronUpElectronDownNuclearParameters();
  EupEupNuclear = JP.getElectronUpElectronUpNuclearParameters();
  EdnEdnNuclear = JP.getElectronDownElectronDownNuclearParameters();

  for(int Nuclei=0; Nuclei<globalInput.Molecule.getNumberAtoms(); Nuclei++)
    {
      // Find the number of the current nucleus in the nuclei list
      int NucleiType = -1;
      for( int i=0; i<NucleiTypes->dim1(); i++ )
	if( globalInput.Molecule.Atom_Labels(Nuclei) == (*NucleiTypes)(i) )
	  {
	    NucleiType = i;
	    break;
	  } 

      QMCThreeBodyCorrelationFunction *U_Function;

      double r1 = MathFunctions::rij(R,globalInput.Molecule.Atom_Positions,E,Nuclei);
      double r1g[grid.dim1()];
      for(int gp=0; gp<grid.dim1(); gp++)
	r1g[gp] = MathFunctions::rij(grid,globalInput.Molecule.Atom_Positions,gp,Nuclei);

      // Now we do all opposite spin pairs with this nucleus	 
      if((*EupEdnNuclear)(NucleiType).isUsed()){
	U_Function = (*EupEdnNuclear)(NucleiType).getThreeBodyCorrelationFunction();
	if (nalpha > 0 && nbeta > 0)
	  if(isAlpha){
	    for (int EB=nalpha; EB<R.dim1(); EB++){
	      double r2 = MathFunctions::rij(R,globalInput.Molecule.Atom_Positions,EB,Nuclei);
	      denom += U_Function->getFunctionValue(MathFunctions::rij(R,EB,E),r1,r2);
	      for(int gp=0; gp<grid.dim1(); gp++)
		sumU(gp) += U_Function->getFunctionValue(MathFunctions::rij(grid,R,gp,EB),r1g[gp],r2);
	    }
	  } else {
	    for (int EA=0; EA<nalpha; EA++){
	      double r2 = MathFunctions::rij(R,globalInput.Molecule.Atom_Positions,EA,Nuclei);
	      denom += U_Function->getFunctionValue(MathFunctions::rij(R,EA,E),r1,r2);
	      for(int gp=0; gp<grid.dim1(); gp++)
		sumU(gp) += U_Function->getFunctionValue(MathFunctions::rij(grid,R,gp,EA),r1g[gp],r2);
	    }
	  }
      }

      // First we do all pairs of two alphas with this nucleus
      if((*EupEupNuclear)(NucleiType).isUsed()){
	U_Function = (*EupEupNuclear)(NucleiType).getThreeBodyCorrelationFunction();
	if (nalpha > 1)
	  if(isAlpha)
	    {
	      for(int EA=0; EA<nalpha; EA++)
		{
		  if(EA == E) continue;
		  double r2  = MathFunctions::rij(R,globalInput.Molecule.Atom_Positions,EA,Nuclei);
		  denom += U_Function->getFunctionValue(MathFunctions::rij(R,EA,E),r1,r2);
		  for(int gp=0; gp<grid.dim1(); gp++)
		    sumU(gp) += U_Function->getFunctionValue(MathFunctions::rij(grid,R,gp,EA),r1g[gp],r2);
		}
	    }
      }      

      // Finally we do all beta spin pairs with this nucleus
      if((*EdnEdnNuclear)(NucleiType).isUsed()){
	U_Function = (*EdnEdnNuclear)(NucleiType).getThreeBodyCorrelationFunction();
	if (nbeta > 1)
	  if(!isAlpha)
	    {
	      for(int EB=nalpha; EB<R.dim1(); EB++)
		{
		  if(EB == E) continue;
		  double r2  = MathFunctions::rij(R,globalInput.Molecule.Atom_Positions,EB,Nuclei);
		  denom += U_Function->getFunctionValue(MathFunctions::rij(R,EB,E),r1,r2);
		  for(int gp=0; gp<grid.dim1(); gp++)
		    sumU(gp) += U_Function->getFunctionValue(MathFunctions::rij(grid,R,gp,EB),r1g[gp],r2);
		}
	    }
      }
    }

  for(int gp=0; gp<grid.dim1(); gp++)
    integrand(gp) *= exp(sumU(gp));
  return exp(denom);
}

void QMCThreeBodyJastrow::updateAll(QMCJastrowParameters & JP, Array2D<double> & X)
{
  int nalpha = globalInput.WF.getNumberElectrons(true);
  int nbeta = globalInput.WF.getNumberElectrons(false);

  for(int E1=0; E1<wd->UijA.dim1(); E1++)
    for(int E2=0; E2<E1; E2++)
      {
	wd->U    -= wd->UijA(E1,E2);
	wd->U_xx -= wd->UijA_xx(E1,E2);
	for(int xyz=0; xyz<3; xyz++)
	  {
	    wd->U_x(E1,xyz) -= (wd->UijA_x1(xyz))(E1,E2); 
	    wd->U_x(E2,xyz) -= (wd->UijA_x2(xyz))(E1,E2); 
	  }
      }

  wd->UijA    = 0;
  wd->UijA_xx = 0;
  for (int i=0; i<3; i++)
    {
      wd->UijA_x1(i) = 0.0;
      wd->UijA_x2(i) = 0.0;
    }
  
  // Get values from JP that will be needed during the calc

  Array1D<string> * NucleiTypes = JP.getNucleiTypes();

  EupEdnNuclear = JP.getElectronUpElectronDownNuclearParameters();
  for(int nuc=0; nuc<EupEdnNuclear->dim1(); nuc++)
    (*EupEdnNuclear)(nuc).zeroOutDerivatives();

  EupEupNuclear = JP.getElectronUpElectronUpNuclearParameters();
  for(int nuc=0; nuc<EupEupNuclear->dim1(); nuc++)
    (*EupEupNuclear)(nuc).zeroOutDerivatives();
  
  EdnEdnNuclear = JP.getElectronDownElectronDownNuclearParameters();
  for(int nuc=0; nuc<EdnEdnNuclear->dim1(); nuc++)
    (*EdnEdnNuclear)(nuc).zeroOutDerivatives();

  Array1D<double> xyz(3);
  for(int Nuclei=0; Nuclei<globalInput.Molecule.getNumberAtoms(); Nuclei++)
    {
      // Find the number of the current nucleus in the nuclei list
      int NucleiType = -1;
      for( int i=0; i<NucleiTypes->dim1(); i++ )
	if( globalInput.Molecule.Atom_Labels(Nuclei) == (*NucleiTypes)(i) )
	  {
	    NucleiType = i;
	    break;
	  } 
      
      // Now we do all opposite spin pairs with this nucleus	 
      if (nalpha > 0 && nbeta > 0)
	for (int E2=nalpha; E2<nalpha+nbeta; E2++)
	  {
	    bool used;
	    for(int i=0; i<3; i++)
	      xyz(i)  = wd->riA_uvec(E2,Nuclei,i);
	    used = (*EupEdnNuclear)(NucleiType).getThreeBodyCorrelationFunction()->setElectron(true, xyz, wd->riA(E2,Nuclei));
	    if(!used) continue;
	    
	    for(int E1=0; E1<nalpha; E1++)
	      collectForPair(E2, E1, Nuclei,
			     & (*EupEdnNuclear)(NucleiType));
	  }

      // First we do all pairs of two alphas with this nucleus

      if (nalpha > 1)
	for(int E1=1; E1<nalpha; E1++)
	  {
	    bool used;
	    for(int i=0; i<3; i++)
	      xyz(i)  = wd->riA_uvec(E1,Nuclei,i);
	    used = (*EupEupNuclear)(NucleiType).getThreeBodyCorrelationFunction()->setElectron(true, xyz, wd->riA(E1,Nuclei));
	    if(!used) continue;

	    for(int E2=0; E2<E1; E2++)	      
	      collectForPair(E1, E2, Nuclei,
			     & (*EupEupNuclear)(NucleiType));
	  }

      // Finally we do all beta spin pairs with this nucleus
      if (nbeta > 1)
	for(int E1=nalpha+1; E1<X.dim1(); E1++)
	  {
	    bool used;
	    for(int i=0; i<3; i++)
	      xyz(i)  = wd->riA_uvec(E1,Nuclei,i);
	    used = (*EdnEdnNuclear)(NucleiType).getThreeBodyCorrelationFunction()->setElectron(true, xyz, wd->riA(E1,Nuclei));
	    if(!used) continue;

	    for(int E2=nalpha; E2<E1; E2++)
	      collectForPair(E1, E2, Nuclei,
			     & (*EdnEdnNuclear)(NucleiType));
	  }
    }

  for(int E1=0; E1<wd->UijA.dim1(); E1++)
    for(int E2=0; E2<E1; E2++)
      {
	wd->U    += wd->UijA(E1,E2);
	wd->U_xx += wd->UijA_xx(E1,E2);
	for(int xyz=0; xyz<3; xyz++)
	  {
	    wd->U_x(E1,xyz) += (wd->UijA_x1(xyz))(E1,E2); 
	    wd->U_x(E2,xyz) += (wd->UijA_x2(xyz))(E1,E2); 
	  }
      }

  xyz.deallocate();
  packageDerivatives();
}

inline void QMCThreeBodyJastrow::collectForPair(int E1, 
						int E2,
						int Nuclei,
						QMCThreeBodyCorrelationFunctionParameters * paramset)
{
  if(!paramset->isUsed()) return;

  QMCThreeBodyCorrelationFunction *U_Function = paramset->getThreeBodyCorrelationFunction();

  Array1D<double> xyz2(3);
  Array1D<double> xyz12(3);
  
  for(int i=0; i<3; i++)
    {
      xyz2(i)  = wd->riA_uvec(E2,Nuclei,i);
      xyz12(i) = wd->rij_uvec(E1, E2,i);
    }

  bool used;
  used = U_Function->setElectron(false,xyz2, wd->riA(E2,Nuclei));
  if(!used) return;

  U_Function->evaluate(xyz12, wd->rij(E1,E2));

  if(IeeeMath::isNaN( U_Function->getLaplacianValue() ) && false ){
    printf("(%2i,%2i) sum %20.10e cur %20.10e psi %20.10e r12 %20.10e r1 %20.10e r2 %20.10e \n",E1,E2,
	   wd->UijA_xx(E1,E2),
	   U_Function->getLaplacianValue(),
	   U_Function->getFunctionValue(),
	   wd->rij(E1,E2),wd->riA(E1,Nuclei),wd->riA(E2,Nuclei)
	   );
    U_Function->print(cerr);
    //exit(0);
  }

  wd->UijA(E1,E2)    += U_Function->getFunctionValue();
  wd->UijA_xx(E1,E2) += U_Function->getLaplacianValue();

  Array1D<double> * grad1 = U_Function->getElectron1Gradient();
  Array1D<double> * grad2 = U_Function->getElectron2Gradient();
  for (int i=0; i<3; i++)
    {
      (wd->UijA_x1(i))(E1,E2) += (*grad1)(i);
      (wd->UijA_x2(i))(E1,E2) += (*grad2)(i);
    }

  if(globalInput.flags.calculate_Derivatives == 1 &&
     globalInput.flags.optimize_NEE_Jastrows == 1)
    {
      for(int ai=0; ai<paramset->getNumberOfTotalParameters()+1; ai++)
	{
	  paramset->pt_a(ai)    += U_Function->get_p_a(ai);
	  paramset->pt3_xxa(ai) += U_Function->get_p3_xxa(ai);

	  for(int i=0; i<3; i++)
	    {
	      (paramset->pt2_xa(ai))(E1,i) += U_Function->get_p2_xa(true,i,ai); 
	      (paramset->pt2_xa(ai))(E2,i) += U_Function->get_p2_xa(false,i,ai); 
	    }
	}
    }
}

void QMCThreeBodyJastrow::packageDerivatives()
{
  if(globalInput.flags.calculate_Derivatives != 1 ||
     globalInput.flags.optimize_NEE_Jastrows == 0)
    return;

  int numNEE =
    globalInput.JP.getNumberNEupEdnParameters() +
    globalInput.JP.getNumberNEupEupParameters() + 
    globalInput.JP.getNumberNEdnEdnParameters();  

  p_a.allocate(numNEE);
  p2_xa.allocate(numNEE);
  p3_xxa.allocate(numNEE);
  for(int i=0; i<p2_xa.dim1(); i++)
    {
      p2_xa(i).allocate(globalInput.WF.getNumberElectrons(),3);
      p2_xa(i) = 0.0;
    }
  p_a = 0.0;
  p3_xxa = 0.0;

  int numNuc = EupEdnNuclear->dim1();
  int ai = 0;
  for(int nuc=0; nuc<numNuc; nuc++)
    {
      QMCThreeBodyCorrelationFunctionParameters * accum = & (*EupEdnNuclear)(nuc);

      switch(globalInput.flags.link_NEE_Jastrows)
	{
	case 0:
	  accum->totalDerivativesToFree(ai,p_a,p2_xa,p3_xxa);
	  accum = & (*EupEupNuclear)(nuc);
	  accum->totalDerivativesToFree(ai,p_a,p2_xa,p3_xxa);
	  accum = & (*EdnEdnNuclear)(nuc);
	  accum->totalDerivativesToFree(ai,p_a,p2_xa,p3_xxa);
	  break;

	case 1:
	  accum->totalDerivativesToFree(ai,p_a,p2_xa,p3_xxa);
	  accum = & (*EupEupNuclear)(nuc);
	  accum->totalDerivativesAccumulate((*EdnEdnNuclear)(nuc));
	  accum->totalDerivativesToFree(ai,p_a,p2_xa,p3_xxa);
	  break;

	case 2:
	  accum->totalDerivativesAccumulate((*EupEupNuclear)(nuc));
	  accum->totalDerivativesAccumulate((*EdnEdnNuclear)(nuc));
	  accum->totalDerivativesToFree(ai,p_a,p2_xa,p3_xxa);
	  break;
	}
    }
}

