// -*- C++ -*-
//
// ReweightMinPT.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ReweightMinPT class.
//

#include "ReweightMinPT.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

IBPtr ReweightMinPT::clone() const {
  return new_ptr(*this);
}

IBPtr ReweightMinPT::fullclone() const {
  return new_ptr(*this);
}

double ReweightMinPT::weight() const {
  Energy minPt = Constants::MaxEnergy;
  for ( int i = 0, N = subProcess()->outgoing().size(); i < N; ++i )
    if ( !onlyColoured || subProcess()->outgoing()[i]->coloured() )
      minPt = min(minPt, subProcess()->outgoing()[i]->momentum().perp());
  return pow(minPt/scale, power);
}

void ReweightMinPT::persistentOutput(PersistentOStream & os) const {
  os << power << ounit(scale,GeV) << onlyColoured;
}

void ReweightMinPT::persistentInput(PersistentIStream & is, int) {
  is >> power >> iunit(scale,GeV) >> onlyColoured;
}

ClassDescription<ReweightMinPT> ReweightMinPT::initReweightMinPT;
// Definition of the static class description member.

void ReweightMinPT::Init() {

  static ClassDocumentation<ReweightMinPT> documentation
    ("There is no documentation for the ThePEG::ReweightMinPT class");

  static Parameter<ReweightMinPT,double> interfacePower
    ("Power",
     "The power to which the minimum tranverse momentum (divided by a "
     "<interface>Scale</interface>) is raised to give the weight.",
     &ReweightMinPT::power, 4.0, -10.0, 10.0, false, false, true);

  static Parameter<ReweightMinPT,Energy> interfaceScale
    ("Scale",
     "The scale with which the minimum transverse momentum is divided "
     "befor it is raised to a <interface>Power</interface> to give the "
     "weight..",
     &ReweightMinPT::scale, GeV, 50.0*GeV, ZERO, ZERO,
     false, false, Interface::lowerlim);


  static Switch<ReweightMinPT,bool> interfaceOnlyColoured
    ("OnlyColoured",
     "Only consider coloured particles in the SubProcess when finding the minimum transverse momentum.",
     &ReweightMinPT::onlyColoured, false, true, false);
  static SwitchOption interfaceOnlyColouredTrue
    (interfaceOnlyColoured,
     "True",
     "Use only coloured particles.",
     true);
  static SwitchOption interfaceOnlyColouredFalse
    (interfaceOnlyColoured,
     "False",
     "Use all particles.",
     false);
  static SwitchOption interfaceOnlyColouredYes
    (interfaceOnlyColoured,
     "Yes",
     "Use only coloured particles.",
     true);
  static SwitchOption interfaceOnlyColouredNo
    (interfaceOnlyColoured,
     "No",
     "Use all particles.",
     false);


  interfacePower.rank(10);
  interfaceScale.rank(9);
  interfaceOnlyColoured.rank(8);

}

