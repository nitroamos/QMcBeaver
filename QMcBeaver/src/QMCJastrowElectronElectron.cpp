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

#include "QMCJastrowElectronElectron.h"
#include "MathFunctions.h"

QMCJastrowElectronElectron::QMCJastrowElectronElectron(){}

QMCJastrowElectronElectron::~QMCJastrowElectronElectron()
{
  for(int i=0; i<p2_xa.dim1(); i++)
    p2_xa(i).deallocate();

  p_a.deallocate();
  p2_xa.deallocate();
  p3_xxa.deallocate();
}

void QMCJastrowElectronElectron::initialize(QMCInput * input)
{
  wd = 0;

  int numEE = globalInput.JP.getNumberEEParameters();
  p_a.allocate(numEE);
  p2_xa.allocate(numEE);
  p3_xxa.allocate(numEE);

  for(int i=0; i<p2_xa.dim1(); i++)
    p2_xa(i).allocate(globalInput.WF.getNumberElectrons(),3);
}

double QMCJastrowElectronElectron::get_p3_xxa_ln(int ai)
{
  return p3_xxa(ai);
}

Array2D<double> * QMCJastrowElectronElectron::get_p2_xa_ln(int ai)
{
  return &p2_xa(ai);
}

double QMCJastrowElectronElectron::get_p_a_ln(int ai)
{
  return p_a(ai);
}

void QMCJastrowElectronElectron::evaluate(QMCJastrowParameters & JP,
					  QMCWalkerData * wData,
					  Array2D<double> & X)
{
  wd = wData;
  if(wd->whichE != -1)
    updateOne(JP,X);
  else
    updateAll(JP,X);
  wd = 0;
}

void QMCJastrowElectronElectron::updateOne(QMCJastrowParameters & JP,
					   Array2D<double> & X)
{
  int E = wd->whichE;

  int nalpha = globalInput.WF.getNumberElectrons(true);
  int nbeta  = globalInput.WF.getNumberElectrons(false);
  bool isAlpha = E < nalpha ? true : false;
  
  // Get values from JP that will be needed during the calc

  QMCCorrelationFunctionParameters * EupEdn = 0;

  if(nalpha > 0 && nbeta > 0)
    EupEdn = JP.getElectronUpElectronDownParameters();
  
  QMCCorrelationFunctionParameters * EupEup = 0;
  
  if(nalpha > 1)
    EupEup = JP.getElectronUpElectronUpParameters();
  
  QMCCorrelationFunctionParameters * EdnEdn = 0;
  
  if(nbeta > 1 )
    EdnEdn = JP.getElectronDownElectronDownParameters();

  QMCCorrelationFunction *U_Function = 0;

  int index = 0;
  int numP  = 0;

  if(EupEdn != 0)
    {
      U_Function = EupEdn->getCorrelationFunction();
      numP = EupEdn->getTotalNumberOfParameters();

      if(isAlpha)
	{
	  for(int EB=nalpha; EB<X.dim1(); EB++)
	    collectForPair(EB,E,U_Function,X,index,numP);
	} else {
	  for(int EA=0; EA<nalpha; EA++)
	    collectForPair(E,EA,U_Function,X,index,numP);
	}
    }
  
  index += numP;
  
  if(EupEup != 0)
    {
      numP  = EupEup->getTotalNumberOfParameters();
      U_Function = EupEup->getCorrelationFunction();

      if(isAlpha)
	{
	  for(int EA=0; EA<E; EA++)
	    collectForPair(E,EA,U_Function,X,index,numP);
	  for(int EA=E+1; EA<nalpha; EA++)
	    collectForPair(EA,E,U_Function,X,index,numP);
	}
    } else {
      numP = 0;
    }

  if(globalInput.flags.link_Jastrow_parameters == 0)
    index += numP;

  if(EdnEdn != 0)
    {
      numP  = EdnEdn->getTotalNumberOfParameters();
      U_Function = EdnEdn->getCorrelationFunction();

      if(!isAlpha)
	{
	  for(int EB=nalpha; EB<E; EB++)
	    collectForPair(E,EB,U_Function,X,index,numP);      
	  for(int EB=E+1; EB<X.dim1(); EB++)
	    collectForPair(EB,E,U_Function,X,index,numP);
	}
    } else {
      numP  = 0;
    }
}

void QMCJastrowElectronElectron::updateAll(QMCJastrowParameters & JP,
					   Array2D<double> & X)
{
  // initialize the results
  p_a              = 0.0;
  p3_xxa           = 0.0;
  for(int ai=0; ai<p2_xa.dim1(); ai++)
    p2_xa(ai)      = 0.0;

  int nalpha = globalInput.WF.getNumberElectrons(true);
  int nbeta  = globalInput.WF.getNumberElectrons(false);

  // Get values from JP that will be needed during the calc

  QMCCorrelationFunctionParameters * EupEdn = 0;

  if(nalpha > 0 && nbeta > 0)
    EupEdn = JP.getElectronUpElectronDownParameters();
  
  QMCCorrelationFunctionParameters * EupEup = 0;
  
  if(nalpha > 1)
    EupEup = JP.getElectronUpElectronUpParameters();
  
  QMCCorrelationFunctionParameters * EdnEdn = 0;
  
  if(nbeta > 1 )
    EdnEdn = JP.getElectronDownElectronDownParameters();

  QMCCorrelationFunction *U_Function = 0;

  int index = 0;
  int numP  = 0;

  if(EupEdn != 0)
    {
      U_Function = EupEdn->getCorrelationFunction();
      numP = EupEdn->getTotalNumberOfParameters();

      for(int E1=0; E1<nalpha; E1++)
	for(int E2=nalpha; E2<X.dim1(); E2++)
	  collectForPair(E2,E1,U_Function,X,index,numP);
    }
  
  index += numP;
  
  if(EupEup != 0)
    {
      numP  = EupEup->getTotalNumberOfParameters();
      U_Function = EupEup->getCorrelationFunction();

      for(int E1=0; E1<nalpha; E1++)
	for(int E2=0; E2<E1; E2++)
	  collectForPair(E1,E2,U_Function,X,index,numP);
    } else {
      numP = 0;
    }
  
  if(globalInput.flags.link_Jastrow_parameters == 0)
    index += numP;

  if(EdnEdn != 0)
    {
      numP  = EdnEdn->getTotalNumberOfParameters();
      U_Function = EdnEdn->getCorrelationFunction();

      for(int E1=nalpha; E1<X.dim1(); E1++)
	for(int E2=nalpha; E2<E1; E2++)
	  collectForPair(E1,E2,U_Function,X,index,numP);
    } else {
      numP  = 0;
    }
}

inline void QMCJastrowElectronElectron::collectForPair(int E1, 
						       int E2,
						       QMCCorrelationFunction *U_Function,
						       Array2D<double> & X,
						       int index, int numP)
{
  double temp;
  double Uij, Uij_x, Uij_xx;
  double r = wd->rij(E1, E2);

  U_Function->evaluate(r);      
  Uij    = U_Function->getFunctionValue();
  Uij_x  = U_Function->getFirstDerivativeValue();
  Uij_xx = 2.0*(2.0/r * Uij_x + U_Function->getSecondDerivativeValue());
  
  //Subtract the old value
  wd->U -= wd->Uij(E1, E2);
  //Add the new value
  wd->U += Uij;
  //Save the new value
  wd->Uij(E1, E2)    = Uij;
  
  for(int i=0; i<3; i++)
    {
      temp = wd->Uij_x(E1, E2, i);
      wd->U_x(E1,i) -= temp;
      wd->U_x(E2,i) += temp;
      
      temp = Uij_x * wd->rij_uvec(E1, E2,i);
      wd->Uij_x(E1, E2, i) = temp;
      
      wd->U_x(E1,i) += temp;
      wd->U_x(E2,i) -= temp;
    }
  
  wd->U_xx -= wd->Uij_xx(E1, E2);
  wd->U_xx += Uij_xx;
  wd->Uij_xx(E1, E2) = Uij_xx;

  /*
    These are for calculating the derivative of the energy
    with respect to parameter ai.
  */
  if(globalInput.flags.calculate_Derivatives == 1 &&
     globalInput.flags.optimize_EE_Jastrows == 1)
    {
      double firstDeriv;

      for(int ai=0; ai<numP; ai++)
	{	  
	  p_a(ai+index)    += U_Function->get_p_a(ai);
	  firstDeriv        = U_Function->get_p2_xa(ai);
	  p3_xxa(ai+index) += 2.0*(2.0/r * firstDeriv +
				   U_Function->get_p3_xxa(ai));
	  
	  for(int i=0; i<3; i++)
	    {
	      temp = firstDeriv * wd->rij_uvec(E1, E2,i);
	      
	      (p2_xa(ai+index))(E1,i) += temp;
	      (p2_xa(ai+index))(E2,i) -= temp;
	    }
	}
    }
}

double QMCJastrowElectronElectron::jastrowOnGrid(QMCJastrowParameters & JP,
						 int E,
						 Array2D<double> & R,
						 Array2D<double> & grid,
						 Array1D<double> & integrand)
{
  int nalpha = globalInput.WF.getNumberElectrons(true);
  int nbeta  = globalInput.WF.getNumberElectrons(false);
  double denom = 0.0;

  bool isAlpha = E < nalpha ? true : false;

  Array1D<double> sumU(integrand.dim1());
  sumU = 0.0;

  if(integrand.dim1() != grid.dim1()){
    cout << "Integrand dimensions don't match grid!\n";
  }

  // Get values from JP that will be needed during the calc

  QMCCorrelationFunctionParameters * EupEdn = 0;

  if(nalpha > 0 && nbeta > 0)
    EupEdn = JP.getElectronUpElectronDownParameters();
  
  QMCCorrelationFunctionParameters * EupEup = 0;
  
  if(nalpha > 1)
    EupEup = JP.getElectronUpElectronUpParameters();
  
  QMCCorrelationFunctionParameters * EdnEdn = 0;
  
  if(nbeta > 1 )
    EdnEdn = JP.getElectronDownElectronDownParameters();

  QMCCorrelationFunction *U_Function = 0;

  if(EupEdn != 0)
    {
      U_Function = EupEdn->getCorrelationFunction();

      if(isAlpha)
	{
	  for(int EB=nalpha; EB<R.dim1(); EB++)
	    {
	      denom += U_Function->getFunctionValue(MathFunctions::rij(R,EB,E));
	      for(int gr=0; gr<grid.dim1(); gr++)		
		sumU(gr) += U_Function->getFunctionValue(MathFunctions::rij(grid,R,gr,EB));
	    }
	} else {
	  for(int EA=0; EA<nalpha; EA++)
	    {
	      denom += U_Function->getFunctionValue(MathFunctions::rij(R,E,EA));
	      for(int gr=0; gr<grid.dim1(); gr++)
		sumU(gr) += U_Function->getFunctionValue(MathFunctions::rij(grid,R,gr,EA));
	    }
	}
    }  
  
  if(EupEup != 0)
    {
      U_Function = EupEup->getCorrelationFunction();

      if(isAlpha)
	{
	  for(int EA=0; EA<nalpha; EA++)
	    {
	      if(EA == E) continue;
	      denom += U_Function->getFunctionValue(MathFunctions::rij(R,E,EA));	  
	      for(int gr=0; gr<grid.dim1(); gr++)
		sumU(gr) += U_Function->getFunctionValue(MathFunctions::rij(grid,R,gr,EA));
	    }
	}
    }

  if(EdnEdn != 0)
    {
      U_Function = EdnEdn->getCorrelationFunction();

      if(!isAlpha)
	{
	  for(int EB=nalpha; EB<R.dim1(); EB++)
	    {
	      if(EB == E) continue;
	      denom += U_Function->getFunctionValue(MathFunctions::rij(R,E,EB));
	      for(int gr=0; gr<grid.dim1(); gr++)
		sumU(gr) += U_Function->getFunctionValue(MathFunctions::rij(grid,R,gr,EB));	    
	    }            
	}
    }

  for(int gr=0; gr<grid.dim1(); gr++)
    integrand(gr) *= exp(sumU(gr));
  return exp(denom);
}
