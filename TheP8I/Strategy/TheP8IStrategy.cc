// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TheP8IStrategy class.
//

#include "TheP8IStrategy.h"
#include "ThePEG/Interface/ClassDocumentation.h"


using namespace TheP8I;

TheP8IStrategy::~TheP8IStrategy() {}

IBPtr TheP8IStrategy::clone() const {
  return new_ptr(*this);
}

IBPtr TheP8IStrategy::fullclone() const {
  return new_ptr(*this);
}

NoPIOClassDescription<TheP8IStrategy> TheP8IStrategy::initTheP8IStrategy;

void TheP8IStrategy::Init() {
  static ClassDocumentation<TheP8IStrategy> interfaceDescription
    ("This class represents the default TheP8I strategy",
     "the default Pythia8\\cite{Sjostrand:2007gs} strategy",
     "\\bibitem{Sjostrand:2007gs} T.~Sj\\\"ostrand, S.~Mrenna, P.~Skands, "
     "Comput.~Phys.~Commun.\\ {\\bf 178} (2008) 852.");


}

