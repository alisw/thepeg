// -*- C++ -*-
//
// RemnantDecayer.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RemnantDecayer class.
//

#include "RemnantDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/PDT/RemnantData.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

RemnantDecayer::~RemnantDecayer() {}

bool RemnantDecayer::accept(const DecayMode &) const {
  return false;
}

bool RemnantDecayer::needsFullStep() const {
  return true;
}

bool RemnantDecayer::
canHandle(tcPDPtr, tcPDPtr) const {
  return false;
}

bool RemnantDecayer::
checkExtract(tcPPtr, tcPPtr, const LorentzMomentum & pnew) const {
  return pnew.e() > ZERO;
}

bool RemnantDecayer::multiCapable() const {
  return false;
}

ParticleVector RemnantDecayer::
decay(const DecayMode &, const Particle &) const {
  return ParticleVector();
}

ParticleVector RemnantDecayer::
decay(const DecayMode & dm, const Particle & p, Step &) const {
  ParticleVector children;
  tcRemPPtr remnant = dynamic_ptr_cast<tcRemPPtr>(&p);
  if ( !remnant ) return children;
  tRemPDPtr rpd = data(remnant);
  PVector ex = extracted(remnant);
  tcPPtr par = parent(remnant);
  if ( !par || ex.empty() || !rpd ) return children;
  children= dm.produceProducts();
  return children;
}

void RemnantDecayer::
fillSubSystem(tPPtr p, set<tPPtr> & sub) const {

  // If this particle has already been added we're done.
  if ( member(sub, p) ) return;

  if ( respectDISKinematics() && !LeptonMatcher::Check(p->data()) &&
       p->momentum().m2() < ZERO ) {
    // If this particle belongs to an electro-weak scattering vertex
    // it should be excluded. (more specifically part of a vertex
    // where the other particles are incoming and outgoing leptons
    // while this is not a lepton and is space-like)
    if ( p->children().size() == 1 && p->children()[0]->parents().size() == 2 &&
	 LeptonMatcher::Check(p->children()[0]->data()) &&
	 ( p->children()[0]->momentum().m2() >= ZERO ||
	   p->children()[0]->children().empty() ) &&
	 ( LeptonMatcher::Check(p->children()[0]->parents()[0]->data()) ||
	   LeptonMatcher::Check(p->children()[0]->parents()[1]->data()) ) )
      return;
    if ( p->parents().size() == 1 && p->parents()[0]->children().size() == 2 &&
	 LeptonMatcher::Check(p->parents()[0]->data()) &&
	 ( ( LeptonMatcher::Check(p->parents()[0]->children()[0]->data()) &&
	     ( p->parents()[0]->children()[0]->momentum().m2() >= ZERO ||
	       p->parents()[0]->children()[0]->children().empty() ) ) ||
	   ( LeptonMatcher::Check(p->parents()[0]->children()[1]->data())  &&
	     ( p->parents()[0]->children()[1]->momentum().m2() >= ZERO ||
	       p->parents()[0]->children()[1]->children().empty() ) ) ) )
      return;
  }

  sub.insert(p);

  // Fill in children.
  if ( p->next() ) fillSubSystem(p->next(), sub);
  for ( int i = 0, N = p->children().size(); i < N; ++i )
    fillSubSystem(p->children()[i], sub);

  // Fill also parents, but only if this is not an incoming parton
  // coming from a particle with a soft remnant. Also if this is an
  // incoming particle (with no parents) is is ignored.
  if ( p->previous() ) fillSubSystem(p->previous(), sub);
  for ( int i = 0, N = p->parents().size(); i < N; ++i ) {
    tPPtr parent = p->parents()[i];
    if ( member(sub, parent) ) continue;
    if ( parent->parents().empty() ) continue;
    for ( int j = 0, M = parent->children().size(); j < M; ++j )
      if ( dynamic_ptr_cast<tcRemPPtr>(parent->children()[j]) ) {
	parent = tPPtr();
	break;
      }
    if ( parent ) fillSubSystem(parent, sub);
  }
}

tPVector RemnantDecayer::getSubSystem(tcPPtr parent, tPPtr parton) const {
  tPVector ret;
  Axis dir = parent->momentum()/GeV;
  set<tPPtr> sub;
  fillSubSystem(parton, sub);
  multimap<double,tPPtr> ordsub;
  for ( set<tPPtr>::iterator it = sub.begin(); it != sub.end(); ++it ) {
    if ( (**it).children().size() || (**it).next() ) continue;
    try {
      ordsub.insert(make_pair(-(**it).momentum().rapidity(dir), *it));
    }
    catch ( ... ) {
      Throw<SubSystemFail>()
	<< "Could not find subsystem for shuffling momenta in RemnantDecayer."
	<< Exception::eventerror;
    }
  }
  ret.reserve(ordsub.size());
  for ( multimap<double,tPPtr>::iterator it = ordsub.begin();
	it != ordsub.end(); ++it )  ret.push_back(it->second);

  return ret;
}

LorentzRotation RemnantDecayer::
getZBoost(const LorentzMomentum & p0, const LorentzMomentum & p) {
  LorentzRotation R;
  if ( p.z() > ZERO )
    R.setBoostZ((sqr(p.plus()) - sqr(p0.plus()))/
		(sqr(p.plus()) + sqr(p0.plus())));
  else
    R.setBoostZ((sqr(p0.minus()) - sqr(p.minus()))/
		(sqr(p0.minus()) + sqr(p.minus())));
  return R;
}

tPVector RemnantDecayer::
decayRemnants(const tPVector & particles, Step & step) {
  tPVector final;
  tPVector remnants;
  for ( int i = 0, N= particles.size(); i < N; ++i )
    if ( dynamic_ptr_cast<tRemPPtr>(particles[i]) )
      remnants.push_back(particles[i]);
    else
      final.push_back(particles[i]);

  while ( !remnants.empty() ) {
    int i = UseRandom::irnd(remnants.size());
    ParticleVector children = Decayer::DecayParticle(remnants[i], step);
    final.insert(final.end(), children.begin(), children.end());
    remnants.erase(remnants.begin() + i);
  }

  // The final particles may have received a pt-kick and would then be
  // copied. Go through and find the final particles.
  for ( int i = 0, N = final.size(); i < N; ++ i) final[i] = final[i]->final();

  return final;
}

bool RemnantDecayer::preInitialize() const {
  return Decayer::preInitialize() || !pTGenerator();
}

void RemnantDecayer::doinit() {
  Decayer::doinit();
  if ( pTGenerator() ) return;
  thePTGenerator = dynamic_ptr_cast<PtGPtr>
    (generator()->preinitCreate("ThePEG::GaussianPtGenerator",
				fullName() + "/PtGen",
				"GaussianPtGenerator.so"));
}

void RemnantDecayer::persistentOutput(PersistentOStream & os) const {
  os << oenum(theRecoilOption) << respectDIS << thePTGenerator;
}

void RemnantDecayer::persistentInput(PersistentIStream & is, int) {
  is >> ienum(theRecoilOption) >> respectDIS >> thePTGenerator;
}

AbstractClassDescription<RemnantDecayer> RemnantDecayer::initRemnantDecayer;
// Definition of the static class description member.

void RemnantDecayer::Init() {

  static ClassDocumentation<RemnantDecayer> documentation
    ("The RemnantDecayer class is the base class to be used for all "
     "decayers capable of decaying a RemnantParticle object produced by a "
     "SoftRemnantHandler object.");

  static Switch<RemnantDecayer,RecoilOption> interfaceRecoilOption
    ("RecoilOption",
     "Different options for how to distribute recoils in the hard subsystem "
     "when taking energy to produce remnants.",
     &RemnantDecayer::theRecoilOption, copyFinal, true, false);
  static SwitchOption interfaceRecoilOptionCopyFinal
    (interfaceRecoilOption,
     "CopyFinal",
     "Boost copies of final state particles in hard subsystem.",
     copyFinal);
  static SwitchOption interfaceRecoilOptionBoostFinal
    (interfaceRecoilOption,
     "BoostFinal",
     "Boost only final state particles in hard subsystem.",
     boostFinal);
  static SwitchOption interfaceRecoilOptionBoostAll
    (interfaceRecoilOption,
     "BoostAll",
     "Boost all particles in the hard subsystem.",
     boostAll);

  static Switch<RemnantDecayer,int> interfaceRespectDISKinematics
    ("RespectDISKinematics",
     "If true, do not boost a scattered lepton (and possible radiated "
     "photons) in a DIS event, to ensure that \f$x\f$ and \f$Q^2\f$ is "
     "unmodified.",
     &RemnantDecayer::respectDIS, 2, true, false);
  static SwitchOption interfaceRespectDISKinematicsYes
    (interfaceRespectDISKinematics,
     "Yes",
     "Do not boost scattered lepton. If that doesn't work, include "
     "the scattered lepton and emit a warning message.",
     1);
  static SwitchOption interfaceRespectDISKinematicsNo
    (interfaceRespectDISKinematics,
     "No",
     "Boost scattered lepton together with the rest of the hard subsystem.",
     0);
  static SwitchOption interfaceRespectDISKinematicsTry
    (interfaceRespectDISKinematics,
     "Try",
     "Do not boost scattered lepton. If that doesn't work, silently "
     "include the scattered lepton.",
     2);

  static Reference<RemnantDecayer,PtGenerator> interfacePTGenerator
    ("PTGenerator",
     "An object capable of generating an intrinsic transverse momentum of "
     "the created remnants. If not set and the controlling EventGenerator "
     "has a default PtGenerator object, this will be used. Otherwise a "
     "GaussianPtGenerator object created with default settings in the "
     "initialization will be used instead.",
     &RemnantDecayer::thePTGenerator, true, false, true, true, true);

}

