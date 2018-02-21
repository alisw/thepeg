// -*- C++ -*-
//
// WeakToHadronsDecayer.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WeakToHadronsDecayer class.
//

#include "WeakToHadronsDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

IBPtr WeakToHadronsDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr WeakToHadronsDecayer::fullclone() const {
  return new_ptr(*this);
}

bool WeakToHadronsDecayer::accept(const DecayMode & dm) const {
  if ( dm.orderedProducts().size() < 3 || !dm.productMatchers().empty() ||
       !dm.cascadeProducts().empty() || dm.wildProductMatcher() ) return false;
  int Nqq = 0;
  if ( LeptonMatcher::Check(*dm.orderedProducts()[0]) ) {
    if ( !LeptonMatcher::Check(*dm.orderedProducts()[1]) ) return false;
    if ( dm.orderedProducts()[0]->iCharge() &&
	 dm.orderedProducts()[1]->iCharge() ) return false;
    if ( !dm.orderedProducts()[0]->iCharge() &&
	 !dm.orderedProducts()[1]->iCharge() ) return false;
  }

  for ( int i = 0, N = dm.orderedProducts().size(); i < N; ++i ) {
    if ( !dm.orderedProducts()[i]->coloured() ) continue;
    ++Nqq;
    if ( i == N - 1 ) return false;
    if ( dm.orderedProducts()[i]->iColour() +
	 dm.orderedProducts()[i + 1]->iColour() ) return false;
    if ( i == N - 2 && DiquarkMatcher::Check(*dm.orderedProducts()[i]) &&
	 DiquarkMatcher::Check(*dm.orderedProducts()[i + 1]) ) return false;
    ++i;
  }
  if ( Nqq > 2 ) return false;
  return true;
}

PVector WeakToHadronsDecayer::
getHadrons(int Nh, tcPDVector quarks) const {
  PVector hadrons;
  int Nq = quarks.size();
  tcPDPtr h = flavourGenerator()->
    alwaysGetHadron(quarks[Nq - 2], quarks[Nq - 1]);
  hadrons.push_back(h->produceParticle());
  if ( Nq == 2 ) return hadrons;
  while ( Nh-- > 0 ) {
    int i = irnd(Nq - 2);
    tcPDPair hq = flavourGenerator()->alwaysGenerateHadron(quarks[i]);
    hadrons.push_back(hq.first->produceParticle());
    quarks[i] = hq.second;
  }
  if ( DiquarkMatcher::Check(*quarks[0]) && DiquarkMatcher::Check(*quarks[1]) )
    return PVector();
  h = flavourGenerator()->alwaysGetHadron(quarks[0], quarks[1]);
  hadrons.push_back(h->produceParticle());
  return hadrons;
}

double WeakToHadronsDecayer::
reweight(const Particle & parent, const PVector & children) const {
  if ( LeptonMatcher::Check(children[0]->data()) ) {
    LorentzMomentum ph;
    int il = 1;
    int in = 0;
    if ( children[0]->data().iCharge() ) {
      il = 0;
      in = 1;
    }
    for ( int i = 2, N = children.size(); i < N; ++i )
      ph += children[i]->momentum();
    double w = (parent.mass()*children[il]->momentum().e())*
      (children[in]->momentum()*ph)*16.0/sqr(sqr(parent.mass()));
    return w;
  }
  return 1.0;
}

void WeakToHadronsDecayer::persistentOutput(PersistentOStream &) const {}

void WeakToHadronsDecayer::persistentInput(PersistentIStream &, int) {}

ClassDescription<WeakToHadronsDecayer>
WeakToHadronsDecayer::initWeakToHadronsDecayer;
// Definition of the static class description member.

void WeakToHadronsDecayer::Init() {

  static ClassDocumentation<WeakToHadronsDecayer> documentation
    ("This class performs weak decays of taus and charmed and bottom "
     "hadrons. The intermediate W can either decay leptonically in which "
     "case standard V-A matrix element is used, or it can decay into "
     "quarks in which case the conversion into quarks is performed as for "
     "the QuarkToHadronsDecayer base class. In both cases the W decay "
     "products should be specified first. The spectator system can either "
     "be specified in terms of hadrons or in terms of quarks which will "
     "be collapsed into a single hadron.");

}
