#ifndef Cambridge2CorrelationFunction_H
#define Cambridge2CorrelationFunction_H

#include "QMCCorrelationFunction.h"
#include "QMCPolynomial.h"


/** 
    This correlation function is from the paper:
    Jastrow correlation factor for atoms, molecules, and solids
    Drummond, Towler, Needs Phys Rev B 70, 235119 (2004)
*/

class Cambridge2CorrelationFunction : public QMCCorrelationFunction
{
 protected:
  double FunctionValue;
  double dFunctionValue;
  double d2FunctionValue;
  
 private:    
  /*
    if this is false, then we will return zeros for all values.
    i think this is more useful than quitting the program since this
    jastrow might just represent a correlated sampling run, and bad
    jastrows would be ignored elsewhere.

    the only real problem would be a bad jastrow used to guide
    the calculation.
  */
  bool active;

  //The gamma term for establishing the cusp conditions
  // 1/2 for opposite spin electrons
  // 1/4 for parallel spin electrons
  //-Z_I for electron-nuclui cusp condition
  double g;

  // The "C" parameter. Should be 2 or 3
  int C;

  // The length cutoff
  double L;
  // A scaling factor that might help the L parameter to optimize
  double f;
  double fL;

  /**
     The parameter L is very highly coupled to the other
     parameters. If we allow it to be optimized, it will
     probably mess up the optimization procedure.
  */
  bool optimizeL;

  // alpha_0 from the paper
  double alpha_0;
  double alpha_1;

  // All the alpha parameters in one polynomial
  // alpha_1 is not optimized directly, it depends upon
  // the other parameters.
  Polynomial alpha;

  //the position we evaluated at
  double r;
  
  //Intermediate variables
  double d, d2, dpc;
  double dpc_r, dpc_rr;
  double dpc_L, dpc_Lr, dpc_Lrr;
  double d_a1_dL;

  double P, dP_r, dP_rr;
  double    dP_L, dP_Lr, dP_Lrr;
 public:

  void initializeParameters(Array1D<int> & BeginningIndexOfParameterType, 
			    Array1D<double> &Parameters,
			    Array1D<int> & BeginningIndexOfConstantType, 
			    Array1D<double> & Constants);

  void evaluate( double r );

  bool isSingular();

  Array1D<Complex> getPoles();

  double getFunctionValue();
  double getFunctionValue(double r);

  double get_p_a(int ai);
  
  double getFirstDerivativeValue();
  double get_p2_xa(int ai);
  
  double getSecondDerivativeValue();
  double get_p3_xxa(int ai);

  Array1D<double> getNumeratorCoeffs();

  Array1D<double> getDenominatorCoeffs();

  friend ostream& operator <<(ostream& strm, Cambridge2CorrelationFunction &rhs);

  void print(ostream& strm);
};

#endif
