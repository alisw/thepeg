// -*- C++ -*-
#ifndef THEP8I_PythiaPtr_H
#define THEP8I_PythiaPtr_H

#include "TheP8I/Config/Pythia8Interface.h"

namespace TheP8I {

using namespace ThePEG;

struct PythiaPtr{
    Pythia8Interface * pPtr;
    Pythia8Interface * getpPtr(){
      return pPtr;
    }
    double h;
    void used(){
      ++reuse;
    }
    int reuse;
    void nullify(){
      delete pPtr;
      pPtr = NULL;
    }
    bool operator == (PythiaPtr& other){
      return (this->h==other.h);
    }

  };
}
#endif

#ifndef THEP8I_ComparePythias_H
#define THEP8I_ComparePythias_H

namespace TheP8I {

using namespace ThePEG;

struct ComparePythias{
   bool operator() (PythiaPtr a, PythiaPtr b){
      return (a.h < b.h);
   }
};
}
#endif

#ifndef THEP8I_OverlapPythiaHandler_H
#define THEP8I_OverlapPythiaHandler_H

// This is the container class for creating and handling the Overlap Pythia objects used
// for hadronisation in all overlap string models in TheP8I

#include "TheP8I/Config/Pythia8Interface.h"
#include "ThePEG/Handlers/HadronizationHandler.h"

namespace TheP8I {

using namespace ThePEG;

class OverlapPythiaHandler {

public:
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  OverlapPythiaHandler();

  /**
   * The destructor.
   */
  virtual ~OverlapPythiaHandler();
  //@}

  OverlapPythiaHandler(HadronizationHandler * sf, vector<string>arguments);

  /* Returns a pointer to the appropriate Pythia object, given
     an effective string tension
  */

  Pythia8Interface * GetPythiaPtr(double h, bool recycle=false);

  vector<double> GetPythiaParameters(double h);

protected:

private:

  // Create just a single overlap Pythia with enhancement factor h. This is the memory consuming part, and
  // should thus be handled with care
  void CreateSinglePythia(double h);

  // Tentative -- these three methods should be cleaned up a bit

  bool CalculateEffectiveParameters(double h);

  // TODO: Can this be done analytically? Possibly: int_0^1 f(z) dz = (-1)^a Gamma[1 + a] MeijerG[{{}, {1 + a}}, {{0, 0}, {}}, b m^2]
  double IFragmentationF(double a, double b);

  double PythiaParameter(string par);

  int sign(double num);

  string convert(double d);
  
  // The vector of pointers to the Pythia objects
  vector<PythiaPtr> _overlapPythias;

  // The list of additional settings for Pythia
  vector<string> _pythiaSettings;

  HadronizationHandler * _sf;
  double a, b, rho, x, y, xi, sigma, m2;
  double a_eff, b_eff, rho_eff, x_eff, y_eff, xi_eff, sigma_eff;
  double _bsparameter;
};

}  


#endif