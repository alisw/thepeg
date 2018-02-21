// -*- C++ -*-
//
// O1AlphaS.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the O1AlphaS class.
//

#include "O1AlphaS.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;

IBPtr O1AlphaS::clone() const {
  return new_ptr(*this);
}

IBPtr O1AlphaS::fullclone() const {
  return new_ptr(*this);
}

double O1AlphaS::value(Energy2 scale, const StandardModelBase &) const {
  Energy2 theScale = scaleFactor()*scale;
  return 12.0*Constants::pi/((33.0-2.0*Nf(theScale))*
		  log(max(theScale, sqr(Q0))/sqr(LambdaQCD(Nf(theScale)))));
}

vector<Energy2> O1AlphaS::flavourThresholds() const {
  if ( !quarkMasses().empty() ) {
    int nMasses = quarkMasses().size();
    if ( nMasses != theMaxFlav )
      throw InitException() << "External masses set in O1AlphaS but not for all flavours.";
  }
  vector<Energy2> thresholds;
  for ( long f = 1; f <= theMaxFlav; ++f ) {
    if ( quarkMasses().empty() ) {
      PDPtr p = getParticleData(f);
      if ( p ) thresholds.push_back(sqr((p->mass() + p->CC()->mass())/2.));
    } else {
      thresholds.push_back(sqr(quarkMasses()[f-1]));
    }
  }
  std::sort(thresholds.begin(), thresholds.end());
  return thresholds;
}

vector<Energy> O1AlphaS::LambdaQCDs() const {
  vector<Energy> lambdas(theMaxFlav + 1);
  vector<Energy2> thresholds = flavourThresholds();
  lambdas[theLambdaFlavour] = theLambdaQCD;
  for ( int f = theLambdaFlavour - 1; f >= 0; --f ) {
    if ( thresholds[f] > ZERO ) {
      lambdas[f] =
	sqrt(thresholds[f]*
	     exp(-log(thresholds[f]/sqr(lambdas[f + 1]))*
		 (33.0-2.0*(f+1))/(33.0-2.0*f)));
    } else {
      lambdas[f] = lambdas[f + 1];
    }
  }
  for ( int f = theLambdaFlavour; f < theMaxFlav; ++f ) {
    if ( thresholds[f] > ZERO ) {
      lambdas[f + 1] =
	sqrt(thresholds[f]*
	     exp(-log(thresholds[f]/sqr(lambdas[f]))*
		 (33.0-2.0*f)/(33.0-2.0*(f+1))));
    } else {
      lambdas[f + 1] = lambdas[f];
    }
  }
  return lambdas;
}  

void O1AlphaS::persistentOutput(PersistentOStream & os) const {
  os << ounit(theLambdaQCD, GeV) << theLambdaFlavour << theMaxFlav
     << ounit(Q0, GeV);
}

void O1AlphaS::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theLambdaQCD, GeV) >> theLambdaFlavour >> theMaxFlav
     >> iunit(Q0, GeV);
}

ClassDescription<O1AlphaS> O1AlphaS::initO1AlphaS;

void O1AlphaS::Init() {

   static ClassDocumentation<O1AlphaS> documentation
    ("O1AlphaS inherits from AlphaSBase and implements the leading order "
     "running QCD coupling. The value is determined by the "
     "<interface>LambdaQCD</interface> parameter at a given number of "
     "flavours, <interface>LambdaFlav</interface>. Optionally the coupling "
     "can be frozen under some minimum scale, "
     "<interface>FreezeScale</interface> to avoid divergencies or negative "
     "couplings.");

   static Parameter<O1AlphaS,Energy> interfaceLambdaQCD
     ("LambdaQCD",
      "The \\f$\\Lambda_{QCD}\\f$ in GeV for "
      "<interface>LambdaFlav</interface> active "
      "flavours. The value for other numbers of active flavours is derived by "
      "assuming that \\f$\\alpha_S\\f$ is continuous.",
      &O1AlphaS::theLambdaQCD, GeV, 0.25*GeV, ZERO, ZERO,
      false, false, Interface::lowerlim);

   static Parameter<O1AlphaS,int> interfaceMaxFlav
     ("MaxFlav",
      "The maximum number of flavours used to calculate \\f$\\alpha_S\\f$.",
      &O1AlphaS::theMaxFlav, 6, 3, 8,
      false, false, true);

  typedef void (ThePEG::O1AlphaS::*IFN)(int);

  static Parameter<O1AlphaS,int> interfaceLambdaFlavour
    ("LambdaFlav",
     "The number of active flavours for which LambdaQCD is specified.",
     &O1AlphaS::theLambdaFlavour, 4, 3, 8,
     false, false, true, (IFN)0, 0, 0,
     &O1AlphaS::getMaxFlav, 0);

  static Parameter<O1AlphaS,Energy> interfaceFreezeScale
    ("FreezeScale",
     "The scale in units of GeV below which \\f$\\alpha_S\\f$ is frozen.",
     &O1AlphaS::Q0, GeV, ZERO, ZERO, ZERO,
     true, false, Interface::lowerlim);

  interfaceLambdaQCD.rank(10);
  interfaceLambdaFlavour.rank(9);
  interfaceFreezeScale.rank(8);

}

