// -*- C++ -*-
//
// SimpleAlphaEM.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SimpleAlphaEM class.
//

#include "SimpleAlphaEM.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace ThePEG;

IBPtr SimpleAlphaEM::clone() const {
  return new_ptr(*this);
}

IBPtr SimpleAlphaEM::fullclone() const {
  return new_ptr(*this);
}

double SimpleAlphaEM::value(Energy2 scale, const StandardModelBase & sm) const {
  useMe();
  double alem = sm.alphaEM();
  double alem3pi = alem/(3.0*Constants::pi);
  double rpigg;
  double Q2 = scaleFactor()*scale/GeV2;
  if(Q2<=2e-6)
    rpigg=0.;
  else if ( Q2 < 0.09 )
    rpigg = alem3pi*(13.4916 + log(Q2)) + 0.00835*log(1.0 + Q2);
  else if ( Q2 < 9.0 )
    rpigg=alem3pi*(16.3200 + 2.0*log(Q2)) + 0.00238*log(1.0 + 3.927*Q2);
  else if ( Q2 < 10000.0 )
    rpigg=alem3pi*(13.4955 + 3.0*log(Q2)) + 0.00165 + 0.00299*log(1.0 + Q2);
  else
    rpigg=alem3pi*(13.4955 + 3.0*log(Q2)) + 0.00221 + 0.00293*log(1.0 + Q2);
  return alem/(1.0-rpigg);
}

NoPIOClassDescription<SimpleAlphaEM> SimpleAlphaEM::initSimpleAlphaEM;

void SimpleAlphaEM::Init() {
  static ClassDocumentation<SimpleAlphaEM> interfaceDocumentation
    ("This class implements a running \\f$\\alpha_{\\mbox{EM}}\\f$ according "
     "to Buckhardt et al.",
     "In the running of $\\alpha_{EM}$, the parametrization of "
     "H.~Buckhardt et al. was used. See \\cite{KLEISSCERN9808v3pp129}.",
     "\\bibitem{KLEISSCERN9808v3pp129} R.~Kleiss et al, "
     "CERN 89-08, vol.~3, pp 129-131.");
}

