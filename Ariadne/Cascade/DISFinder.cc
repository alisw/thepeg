// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DISFinder class.
//

#include "DISFinder.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/MaxCmp.h"
#include "DipoleState.h"
#include "RemnantParton.h"

using namespace Ariadne5;

DISFinder::DISFinder() {}

DISFinder::~DISFinder() {}

IBPtr DISFinder::clone() const {
  return new_ptr(*this);
}

IBPtr DISFinder::fullclone() const {
  return new_ptr(*this);
}

bool DISFinder::virtualBoson(tPPtr p) const {
  if ( !p->decayed() ) return false;
  if ( p->mass() > ZERO ) return false;
  long id = p->id();
  if ( id == ParticleID::gamma || id == ParticleID::Z0 ||
       id == ParticleID::Wplus || id == ParticleID::Wminus ) return true;
  return false;
}

pair<tPPtr,tPPtr> DISFinder::getScattered(tPPtr inc, SubProcess & sub) const {
  pair<tPPtr,tPPtr> res;
  if ( !LeptonMatcher::Check(inc->data()) ) return res;
  if ( inc->children().size() == 2 ) {
    if ( LeptonMatcher::Check(inc->children()[0]->data()) &&
	 virtualBoson(inc->children()[1]) )
      return make_pair(inc->children()[0], inc->children()[1]);
    if ( LeptonMatcher::Check(inc->children()[1]->data()) &&
	 virtualBoson(inc->children()[0]) )
      return make_pair(inc->children()[1], inc->children()[0]);
  }
  if ( inc->children().size() == 1 &&
       LeptonMatcher::Check(inc->children()[0]->data()) ) {
    tPPtr sc = inc->children()[0];
    if ( sc->parents().size() == 2 ) {
      if ( sc->parents()[0] == inc && virtualBoson(sc->parents()[1]) )
	return make_pair(sc, sc->parents()[1]);
      if ( sc->parents()[1] == inc && virtualBoson(sc->parents()[0]) )
	return make_pair(sc, sc->parents()[0]);
    }
  }
  // Simply go through the final state and find a reasonable
  // lepton. If several choose the one closest in rapidity.
  MaxCmp<double, tPPtr> maxrap;
  int dir = 1;
  if ( inc->momentum().z() < ZERO ) dir = -1;
  for ( int i = 0, N = sub.outgoing().size(); i < N; ++i ) {
    tPPtr p = sub.outgoing()[i];
    if ( p->id()*inc->id() <= 0 ) continue;
    if ( p->id() == inc->id() || 
	 ( inc->id()%2 && abs(inc->id()) + 1 == abs(p->id()) ) ||
	 ( !inc->id()%2 && abs(inc->id()) - 1 == abs(p->id()) ) )
      maxrap(rapidity(p->momentum())*dir, p);
  }
  res.first = maxrap.index();

  return res;
}

tRemParPtr DISFinder::
createRemnant(tPPtr inc, tPPtr lep, tPPtr bos, DipoleState & state) const {
  if ( !lep ) return tRemParPtr();
  tRemParPtr rem = state.create<RemnantParton>();
  rem->orig(lep);
  rem->setupParent(inc);
  LorentzMomentum p = inc->momentum() - lep->momentum();
  tcPDPtr pdb;
  if ( bos ) pdb = bos->dataPtr();
  else if ( inc->id() == lep->id() ) pdb = getParticleData(ParticleID::gamma);
  else if ( inc->id() < lep->id() ) pdb = getParticleData(ParticleID::gamma);
  else pdb = getParticleData(ParticleID::Wplus);
  rem->setExtracted(pdb, p);

  return rem;
}

pair<tRemParPtr,tRemParPtr> DISFinder::
findDISLeptons(SubProcess & sub, DipoleState & state) const {
  pair<tRemParPtr,tRemParPtr> res;
  pair<tPPtr,tPPtr> A = getScattered(sub.incoming().first, sub);
  pair<tPPtr,tPPtr> B = getScattered(sub.incoming().second, sub);

  // Something went wrong silently ignore...
  if ( A.first && A.first == B.first ) return res;
  if ( A.second && A.second == B.second ) return res;

  return make_pair(createRemnant(sub.incoming().first,
				 A.first, A.second, state),
		   createRemnant(sub.incoming().second,
				 B.first, B.second, state));
}

tPPtr DISFinder::getQuark(tRemParPtr lepton, SubProcess & sub) const {
  MaxCmp<double,tPPtr> maxrap;
  int dir = 1;
  if ( lepton->parentMomentum().z() < ZERO ) dir = -1;
  int idb = 0;
  if ( lepton->extractedData().id() == ParticleID::Wplus ) idb = 1;
  if ( lepton->extractedData().id() == ParticleID::Wminus ) idb = -1;
  for ( int i = 0, N = sub.outgoing().size(); i < N; ++i ) {
    tPPtr q = sub.outgoing()[i];
    if ( !QuarkMatcher::Check(q->data()) ) continue;
    if ( lepton->originalExtracted() &&
	 member(q->parents(), lepton->originalExtracted()) )
      maxrap(rapidity(q->momentum())*dir, q);
    else if ( !lepton->extracted() &&
	      ( idb == 0 ||
		( q->id()%2 && q->id()*idb > 0 ) ||
		( !q->id()%2 && q->id()*idb < 0 ) ) )
      maxrap(rapidity(q->momentum())*dir, q);
  }
  return maxrap.index();
}

tRemParPtr DISFinder::
createQuarkRemnant(tRemParPtr l, tPPtr q, DipoleState & state) const {
  if ( !q ) return tRemParPtr();
  tRemParPtr rem = state.create<RemnantParton>();
  rem->orig(q);
  if ( l->originalExtracted() ) rem->setupParent(l->originalExtracted());
  else rem->setupParent(&l->extractedData(), l->extractedMomentum());
  long idq = -q->id();
  if ( abs(l->extractedData().id()) == abs(long(ParticleID::Wplus)) ) {
    if ( idq%2 && idq > 0 ) ++idq;
    if ( idq%2 && idq < 0 ) --idq;
    if ( !idq%2 && idq > 0 ) --idq;
    if ( !idq%2 && idq < 0 ) ++idq;
  }
  rem->setExtracted(getParticleData(idq),
		    l->extractedMomentum() - q->momentum());
  return rem;

}

pair<tRemParPtr,tRemParPtr> DISFinder::
findDISQuarks(pair<tRemParPtr,tRemParPtr> leptons,
	      SubProcess & sub, DipoleState & state) const {
  pair<tRemParPtr,tRemParPtr> res;
  pair<tPPtr,tPPtr> quarks;
  if ( leptons.first ) quarks.first = getQuark(leptons.first, sub);
  if ( leptons.second ) quarks.second = getQuark(leptons.second, sub);

  // Something went wrong or no leptons were present.
  if ( quarks.first == quarks.second ) return res;

  return make_pair(createQuarkRemnant(leptons.first, quarks.first, state),
		   createQuarkRemnant(leptons.second, quarks.second, state));

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).



// The following static variable is needed for the type description
// system in ThePEG.
DescribeNoPIOClass<DISFinder,HandlerBase>
describeAriadne5DISFinder("Ariadne5::DISFinder", "libAriadne5.so");

void DISFinder::Init() {

  static ClassDocumentation<DISFinder> documentation
    ("The DISFinder class is responsible for finding DIS-type scattered "
     "leptons in a SubProcess. This base class only looks at the "
     "mother-daughter relationships in the SubProcess. The case where "
     "these relationships are not available is treated very crudely and "
     "more sophisticated methods may implemented in sub-classes.");

}

