// -*- C++ -*-
//
// ME2to2QCD.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ME2to2QCD class.
//

#include "ME2to2QCD.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

ME2to2QCD::~ME2to2QCD() {}

unsigned int ME2to2QCD::orderInAlphaS() const {
  return 2;
}

unsigned int ME2to2QCD::orderInAlphaEW() const {
  return 0;
}

double ME2to2QCD::comfac() const {
  return 32.0*sqr(Constants::pi*SM().alphaS(scale()));
}

void ME2to2QCD::persistentOutput(PersistentOStream & os) const {
  os << theMaxFlavour << theKfac << theKfacA << useInterference;
}

void ME2to2QCD::persistentInput(PersistentIStream & is, int) {
  is >> theMaxFlavour >> theKfac >> theKfacA >> useInterference;
}

AbstractClassDescription<ME2to2QCD> ME2to2QCD::initME2to2QCD;
// Definition of the static class description member.

void ME2to2QCD::Init() {

  static ClassDocumentation<ME2to2QCD> documentation
    ("There is no documentation for the ThePEG::ME2to2QCD class");

  static Parameter<ME2to2QCD,int> interfaceMaxFlavour
    ( "MaxFlavour",
      "The heaviest quark flavour this matrix element is allowed to handle "
      "(if applicable). This applies fo both incoming an outgoing partons.",
      &ME2to2QCD::theMaxFlavour, 5, 0, 8, false, false, true);

  static Parameter<ME2to2QCD,double> interfaceKfac
    ("K-factor",
     "A K-factor for artificially boosting this cross-section.",
     &ME2to2QCD::theKfac, 1.0, 0.0, Constants::MaxFloat, true, false, true);

  static Parameter<ME2to2QCD,double> interfaceKfacA
    ("K-factor-A",
     "A K-factor for artificially boosting the colour annihilation terms "
     "in this cross-section. If less than one, the "
     "<interface>K-factor</interface> will "
     "be used instead.",
     &ME2to2QCD::theKfacA, 1.0, -1.0, Constants::MaxFloat, true, false, true);

  static Switch<ME2to2QCD,bool> interfaceUseInterfecence
    ("Interference",
     "Use of interference terms in the matrix elements. If included the terms "
     "are divided between the differenct possible colour configurations "
     "according to the pole structure of the (string-inspired) matrix "
     "elements for the different colour configurations.",
     &ME2to2QCD::useInterference, true, false, false);

  static SwitchOption interfaceInterferenceOff
    (interfaceUseInterfecence, "Excluded", 
     "Only the parts of the matrix elements with a well-defined colour "
     "structure is used.", false);

  static SwitchOption interfaceInterferenceOn
    (interfaceUseInterfecence, "Included",
     "Use the full matrix element.", true);

  interfaceMaxFlavour.rank(10);

}

