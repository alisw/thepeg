// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KTRapidityCut class.
//

#include "KTRapidityCut.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/MatcherBase.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace ThePEG;

KTRapidityCut::~KTRapidityCut() {}

void KTRapidityCut::describe() const {
  CurrentGenerator::log() 
    << fullName() << ":\n"
    << "KT       = " << theMinKT/GeV << " .. " << theMaxKT/GeV << " GeV\n"
    << "Rapidity = " << theMinRapidity << " .. " << theMaxRapidity << "\n\n";
}

IBPtr KTRapidityCut::clone() const {
  return new_ptr(*this);
}

IBPtr KTRapidityCut::fullclone() const {
  return new_ptr(*this);
}

void KTRapidityCut::persistentOutput(PersistentOStream & os) const {
  os << ounit(theMinKT, GeV) << ounit(theMaxKT, GeV)
     << theMinRapidity << theMaxRapidity << theMatcher;
}

void KTRapidityCut::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theMinKT, GeV) >> iunit(theMaxKT, GeV)
     >> theMinRapidity >> theMaxRapidity >> theMatcher;
}

ClassDescription<KTRapidityCut> KTRapidityCut::initKTRapidityCut;
// Definition of the static class description member.

void KTRapidityCut::Init() {

  static ClassDocumentation<KTRapidityCut> documentation
    ("This is a very simple concrete sub-class of OneCutbase simply "
     "requiring a minimum transverse momentum of any outgoing particle. "
     "It is also possible to require a minimum and maximum rapidity. "
     "Optionally the restrictions only apply to particles matching a "
     "specific matcher object.");

  typedef double (ThePEG::KTRapidityCut::*IGFN)() const;
  typedef void (ThePEG::KTRapidityCut::*ISFN)(double);
  typedef Energy (ThePEG::KTRapidityCut::*IGFNK)() const;
  typedef void (ThePEG::KTRapidityCut::*ISFNK)(Energy);

  static Parameter<KTRapidityCut,Energy> interfaceMinKT
    ("MinKT",
     "The minimum allowed value of the transverse momentum of an outgoing "
     "parton.",
     &KTRapidityCut::theMinKT, GeV, 10.0*GeV, ZERO, Constants::MaxEnergy,
     true, false, Interface::limited,
     (ISFNK)0, (IGFNK)0, (IGFNK)0, &KTRapidityCut::maxKTMin, (IGFNK)0);
  interfaceMinKT.setHasDefault(false);

  static Parameter<KTRapidityCut,Energy> interfaceMaxKT
    ("MaxKT",
     "The maximum allowed value of the transverse momentum of an outgoing "
     "parton. Note that this cut does not increase the efficiency of the phase "
     "space generation, but is only applied as a post-cut.",
     &KTRapidityCut::theMaxKT, GeV, Constants::MaxEnergy, ZERO, ZERO,
     true, false, Interface::lowerlim,
     (ISFNK)0, (IGFNK)0,  &KTRapidityCut::minKTMax, (IGFNK)0, (IGFNK)0);
  interfaceMaxKT.setHasDefault(false);

  static Parameter<KTRapidityCut,double> interfaceMinRapidity
    ("MinRapidity",
     "The minimum allowed rapidity of an outgoing parton. "
     "The rapidity is measured in the lab system.",
     &KTRapidityCut::theMinRapidity,
     -Constants::MaxRapidity, 0, Constants::MaxRapidity,
     true, false, Interface::upperlim,
     (ISFN)0, (IGFN)0, (IGFN)0, &KTRapidityCut::maxRapidityMin, (IGFN)0);
  interfaceMinRapidity.setHasDefault(false);

  static Parameter<KTRapidityCut,double> interfaceMaxRapidity
    ("MaxRapidity",
     "The maximum allowed rapidity of an outgoing parton. "
     "The rapidity is measured in the lab system.",
     &KTRapidityCut::theMaxRapidity,
     Constants::MaxRapidity, -Constants::MaxRapidity, 0,
     true, false, Interface::lowerlim,
     (ISFN)0, (IGFN)0,  &KTRapidityCut::minRapidityMax, (IGFN)0, (IGFN)0);
  interfaceMaxRapidity.setHasDefault(false);

  static Reference<KTRapidityCut,MatcherBase> interfaceMatcher
    ("Matcher",
     "If non-null only particles matching this object will be affected "
     "by the cut.",
     &KTRapidityCut::theMatcher, true, false, true, true, false);

  interfaceMinKT.rank(10);
  interfaceMaxKT.rank(6);
  interfaceMinRapidity.rank(9);
  interfaceMaxRapidity.rank(8);
  interfaceMatcher.rank(7);
}

Energy KTRapidityCut::maxKTMin() const {
  return theMaxKT;
}

Energy KTRapidityCut::minKTMax() const {
  return theMinKT;
}

double KTRapidityCut::maxRapidityMin() const {
  return theMaxRapidity;
}

double KTRapidityCut::minRapidityMax() const {
  return theMinRapidity;
}

double KTRapidityCut::maxRapidityMin(tcPDPtr p) const {
  if ( theMatcher )
    if ( !theMatcher->matches(*p) ) 
      return Constants::MaxRapidity;
  return theMaxRapidity;
}

double KTRapidityCut::minRapidityMax(tcPDPtr p) const {
  if ( theMatcher )
    if ( !theMatcher->matches(*p) ) 
      return -Constants::MaxRapidity;
  return theMinRapidity;
}

Energy KTRapidityCut::minKT(tcPDPtr p) const {
  if ( theMatcher && !theMatcher->matches(*p) ) return ZERO;
  return theMinKT;
}

double KTRapidityCut::minEta(tcPDPtr) const {
  return -Constants::MaxRapidity;
}

double KTRapidityCut::maxEta(tcPDPtr) const {
  return Constants::MaxRapidity;
}

bool KTRapidityCut::passCuts(tcCutsPtr parent,
			     tcPDPtr ptype, LorentzMomentum p) const {
  if ( theMatcher && !theMatcher->matches(*ptype) ) return true;
  if ( p.perp() < theMinKT ) return false;
  if ( p.perp() > theMaxKT ) return false;
  double y = p.rapidity() + parent->Y() + parent->currentYHat();
  if ( y > theMaxRapidity ) return false;
  if ( y < theMinRapidity ) return false;
  return true;
}
