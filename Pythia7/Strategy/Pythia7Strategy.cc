// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Pythia7Strategy class.
//

#include "Pythia7Strategy.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "Pythia7Strategy.tcc"
#endif

using namespace Pythia7;

Pythia7Strategy::~Pythia7Strategy() {}

NoPIOClassDescription<Pythia7Strategy> Pythia7Strategy::initPythia7Strategy;

void Pythia7Strategy::Init() {
  static ClassDocumentation<Pythia7Strategy> interfaceDescription
    ("This class represents the default Pythia7 strategy",
     "the default Pythia\\cite{TS94Pythia} strategy",
     "\\bibitem{TS94Pythia} T.~Sj\\\"ostrand, "
     "Comput.~Phys.~Commun.\\ {\\bf 82} (1994) 74.");


}

