// -*- C++ -*-

#ifndef THEP8I_ParameterHandler_H
#define THEP8I_ParameterHandler_H

#include "TheP8I/Config/TheP8I.h"


namespace TheP8I {

using namespace ThePEG;

class ParameterHandler {


public:
  typedef map<string,double> PytPars;

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ParameterHandler();

  /**
   * The destructor.
   */
  virtual ~ParameterHandler();
  //@}

  void init(double m2, double bsparameter, PytPars pars);
  
  PytPars GetEffectiveParameters(double h);

protected:

private:

  bool InsertEffectiveParameters(double h);

  bool CalculateEffectiveParameters(double h);

  double IFragmentationF(double a, double b);

  // sets of parameters, ordered in h
  map<double, PytPars> _parameters;

  double a, b, rho, x, y, xi, sigma, m2;
  double a_eff, b_eff, rho_eff, x_eff, y_eff, xi_eff, sigma_eff;
  double  _m2, _bsparameter;
};

}  


#endif