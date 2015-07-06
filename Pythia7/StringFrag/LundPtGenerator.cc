// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LundPtGenerator class.
//

#include "LundPtGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"


#ifdef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "LundPtGenerator.tcc"
#endif

using namespace Pythia7;

LundPtGenerator::~LundPtGenerator() {}

TransverseMomentum LundPtGenerator::generate() const {

  TransverseMomentum pT;

  Energy pTmag = Sigma()*sqrt(-log(max(1.e-10, rnd())) );
  
  if ( rndbool(nGaussfraction) ) pTmag *= nGaussfactor; 

  double phi = Constants::twopi*rnd();

  pT.first  =  pTmag*cos(phi);  
  pT.second =  pTmag*sin(phi); 

  return pT;
}


// *** Standard Interfaced functions
void LundPtGenerator::persistentOutput(PersistentOStream & os) const {
  os << ounit(sigma, GeV) << nGaussfraction << nGaussfactor;
}

void LundPtGenerator::persistentInput(PersistentIStream & is, int) {
  is >> iunit(sigma, GeV) >> nGaussfraction >> nGaussfactor;
}


ClassDescription<LundPtGenerator> LundPtGenerator::initLundPtGenerator;

void LundPtGenerator::Init(){

  static ClassDocumentation<LundPtGenerator> documentation
    ("The transverse momentum generator used in the Lund string "
     "fragmentation scheme.");

  static Parameter<LundPtGenerator, Energy> interfaceSigma
    ("sigma",
     "Corresponds to the width in the Gaussian \\f$p_x\\f$ and \\f$p_y\\f$ "
     "transverse momentum distributions for primary hadrons [GeV].",
     &LundPtGenerator::sigma, GeV, 0.36*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<LundPtGenerator, double> interfacenonGaussFraction
    ("nonGaussFraction",
     "The non-Gaussian fraction of the Gaussian transverse momentum "
     "distribution to be enhanced by the factor nGaussfactor [noUnit].",
     &LundPtGenerator::nGaussfraction, 0.01, 0.0, 1.0, false, false, true);

  static Parameter<LundPtGenerator, double> interfacenonGaussFactor
    ("nonGaussFactor",
     "The non-Gaussian tails enhancement factor [noUnit].",
     &LundPtGenerator::nGaussfactor, 2.0, 0.0, 1.0, false, false, true);

  interfaceSigma.rank(10);

}






