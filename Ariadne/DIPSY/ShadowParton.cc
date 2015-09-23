// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShadowParton class.
//

#include "ShadowParton.h"
#include "Parton.h"
#include "DipoleState.h"
#include "DipoleEventHandler.h"
#include "Dipole.h"
#include "ImpactParameters.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Utilities/Current.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/DebugItem.h"

using namespace DIPSY;

ShadowParton::ShadowParton(Parton & p)
  : theOriginal(&p), thePlus(p.plus()),
    thePT(p.pT()), thePT0(p.pT()), theMass(p.mass()), theMinus(p.minus()),
    theY(p.y()), theY0(p.y()), theFlavour(p.flavour()),
    theEmissionFactor(ZERO), theRes(ZERO), isLocked(false),
    hasInteracted(false), isOnShell(false),
    isValence(p.valence()), hasColourSibling(false), memememe(false) {}

ShadowParton::~ShadowParton() {}

void ShadowParton::setupParton() {
  if ( tSPartonPtr sp = original()->shadow() ) {
    sp->theNext = this;
    thePrevious = sp;
  }
  original()->shadow(this);
}

SPartonPtr ShadowParton::createValence(Parton & p) {
  SPartonPtr sp = new_ptr(ShadowParton(p));
  sp->setupParton();
  sp->isValence = true;
  return sp;
}

void ShadowParton::setupEmission(Parton & emitter,
				 Parton & produced, Parton & recoiler) {
  SPartonPtr se = new_ptr(ShadowParton(emitter));
  se->setupParton();
  SPartonPtr sp = new_ptr(ShadowParton(produced));
  sp->setupParton();

  se->theSibling = sp;
  se->hasColourSibling =
    ( emitter.dipoles().first &&
      emitter.dipoles().first == produced.dipoles().second );

  sp->theParent = this;
  sp->theSibling = se;
  sp->hasColourSibling = !se->hasColourSibling;

  theChild = sp;

  InvEnergy2 d2 = se->theRes = sp->theRes = se->dist2(*sp);
  se->theEmissionFactor =
    sp->theEmissionFactor = alphaS(d2)*d2/UseRandom::rnd();
}

void ShadowParton::lock() {
  if ( mother() ) mother()->lock();
  isLocked = true;
}

void ShadowParton::unlock() {
  if ( mother() && locked() ) mother()->unlock();
  isLocked = false;
}

tSPartonPtr ShadowParton::last() {
  if ( next() ) return next();
  return this;
}

tSPartonPtr ShadowParton::initial() {
  if ( previous() ) return previous();
  return this;
}

tSPartonPtr ShadowParton::valenceMother() {
  if ( mother() ) return mother()->valenceMother();
  return this;
} 

tcSPartonPtr ShadowParton::last() const {
  if ( next() ) return next();
  return this;
}

double ShadowParton::alphaS(InvEnergy2 r2) {
  return Current<DipoleEventHandler>()->alphaS(sqrt(r2));
}

// tcSPartonPtr ShadowParton::resolve(InvEnergy2 r2, tPartonPtr stopp) const {
//   if ( resolved(r2, stopp) || !mother() ) return this;
//   return mother()->resolve(r2, stopp);
// }

tSPartonPtr ShadowParton::resolve(InvEnergy2 r2)  {
  if ( resolved(r2) || !mother() ) return this;
  return mother()->resolve(r2);
}
       
tSPartonPtr ShadowParton::resolveInteraction(InvEnergy2 r2)  {
  if ( resolved(r2)
       || !mother()
       || forceEmission()
       || onShell()
       || ( sibling() && sibling()->interacted() ) ) return this;
  return mother()->resolveInteraction(r2);
}
       
void ShadowParton::pTplus(const TransverseMomentum & qt, Energy qp) {
  pT(qt);
  plus(qp);
  minus(mt2()/plus());
  y(log(mt()/plus()));
}

void ShadowParton::setOnShell(int mode) {
  if ( mode < 0 ) return;
  onShell(true);
  if ( mode <= 0 ) return;
  original()->plus(plus());
  original()->pT(pT());
  original()->minus(original()->mt2()/plus());
  original()->y(log(original()->mt()/original()->plus()));
  original()->onShell(true);
}

void ShadowParton::unsetOnShell() {
  if ( onShell() ) {
    original()->onShell(false);
    onShell(false);
  }
}

void ShadowParton::interact() {
  if ( interacted() ) return;
  hasInteracted = true;
  if ( mother() ) mother()->interact();
}

void ShadowParton::resetInteracted() {
  if ( !mother() ) {
    pT(pT0());
    y(y0());
    plus(mt()*exp(-y()));
    minus(mt()*exp(y()));
  }
  hasInteracted = false;
  theInteractions.clear();
  onShell(false);
  if ( next() ) next()->resetInteracted();
  if ( child() ) child()->resetInteracted();
}

void ShadowParton::reset() {
  if ( !mother() ) {
    pT(pT0());
    y(y0());
    plus(mt()*exp(-y()));
    minus(mt()*exp(y()));
  }
  onShell(false);
  theInteractions.clear();
  hasInteracted = false;
  theInteractionRoots.clear();
  if ( next() ) next()->reset();
  if ( child() ) child()->reset();
}

void ShadowParton::insertInteraction(int i) {
  // If this parton has been set on-shell by a previous interaction we
  // do not have to look further. Similarly if the sibling has
  // interacted.
  if ( onShell()
       || ( sibling() && sibling()->interacted() )
       || !mother() ) {
    interactionRoot(i);
  }
  // Otherwise we try with the mother instead.
  else mother()->insertInteraction(i);
}

void ShadowParton::rejectInteraction(int i) {
  if ( hasInteractionRoot(i) ) theInteractionRoots.erase(i);
  else mother()->rejectInteraction(i);
}

void ShadowParton::acceptInteraction(int i) {
  hasInteracted = true;
  if ( !hasInteractionRoot(i) ) mother()->acceptInteraction(i);
}

void ShadowParton::makeIncoming() {
  if ( interacted() && onShell() ) onShell(false);
  else if ( mother() ) mother()->makeIncoming();
}

bool ShadowParton::setEmissionMomenta(const LorentzMomentum & p,
				      bool forced) {
  // *** TODO *** Think this through!  If the emission unresolved but
  // anyway should be done because a valence need to be put on shell,
  // or a subsequent emission needs it, the generated pt is ignored
  // and the emitted parton and the continuing propagator will share
  // the incoming pt, while the generated rapidity of the emission
  // remains untouched.
  TransverseMomentum qT = parent()? pT0(): sibling()->pT0();
  if ( forced ) qT = TransverseMomentum(-p.x()/2.0, -p.y()/2.0);

  if ( parent() ) {
    // If this parton was emitted, pretend it was actually he emitter
    // and continues as a propagator. It will retain its original
    // rapidity, but will get the transverse momentum of the incoming
    // propagator.
    pT(TransverseMomentum(p.x(), p.y()) + qT);
    y(y0());
    plus(mt()*exp(-y()));
    // The sibling will get minus the original transverse momentum and
    // what ever is left of the positive light-cone momenta. (return
    // false if not enough to put on shell).
    sibling()->pT(-qT);
    sibling()->plus(p.plus() - plus());
    if ( sibling()->plus() <= ZERO ) return false;
    sibling()->minus(sibling()->mt2()/sibling()->plus());
    sibling()->y(log(sibling()->mt()/sibling()->plus()));
    // The propagators negative light-cone momentum is trivial,
    minus(p.minus() - sibling()->minus());
  } else {
    // If this parton was assumed emitter in the evolution, life is simpler.
    sibling()->pT(-qT);
    sibling()->y(sibling()->y0());
    sibling()->plus(sibling()->mt()*exp(-sibling()->y()));
    sibling()->minus(sibling()->mt2()/sibling()->plus());
    pT(TransverseMomentum(p.x(), p.y()) + qT);
    plus(p.plus() - sibling()->plus());
    if ( plus() <= ZERO ) return false;
    minus(mt2()/plus());
    y(log(mt()/plus()));
  }
  return true;
}

bool ShadowParton::orderfail(const Propagator & prop) const {
  Energy pos = colourSibling()? prop.acopos: prop.colpos;
  Energy neg = colourSibling()? prop.aconeg: prop.colneg;
  Energy2 ptp = colourSibling()? prop.acoptp: prop.colptp;
  Energy2 ptnow = prop.p.perp2();
  Energy2 ptnxt = (prop.p - sibling()->momentum()).perp2();
  if ( ptnow >= max(ptp, ptnxt) || ptnow <= min(ptp, ptnxt) )
    return false;
  return sibling()->plus() > pos || sibling()->minus() < neg;
}

ShadowParton::Propagator ShadowParton::
propagator(InvEnergy2 r2, int mode) {

  static DebugItem unorderlock("DIPSY::UnorderLock", 6);

  Propagator prop;

  // If we reach a parton that has been set on shell, we will simply
  // return its momentum and set it off-shell.
  if ( onShell() ) {
    if ( mode >= 0 ) unsetOnShell();
    return momentum();
  }
  // If this was a valence, ask the DipoleState about its momentum.
  if ( !mother() ) return dipoleState().incomingMomentum(this, mode);

  // Now, if this emission was not resolved then simply return the
  // propagator of the mother. However, if the emission is needed for
  // subsequent interactions or if the emitted parton was a valence,
  // we need to force it.  with a special procedure.
  bool forced = forceEmission();
  bool unresolved = !resolved(r2);


  // First we check if a resolved emission is possible.
  if ( !unresolved ) {

    // Get the propagator of the mother with the new resolution
    // scale. But don't put anything on-shell yet in case the emission
    // is rejected later.
    prop = mother()->propagator(res(), -1);

    // First check that it was at all kinematically possible perform
    // the emission.
    if ( !prop.fail && setEmissionMomenta(prop.p, false) ) {

      // OK. That seemed to work, but we must also check the ordering.
      // Now depending on the colour line of the emitted parton, it
      // must be ordered in light-cone momenta with previous partons
      // on the same side. If that is not the case, the emissionis
      // marked unresolved.
      if ( !( unorderlock && sibling() && sibling()->locked() ) )
	unresolved = orderfail(prop);
    } else {
      // Since the emission could not be performed, we flag it
      // unresolved.
      unresolved = true;
    }
  }

  // If the emission was really resolved, we get the incoming
  // propagator again, but this time put stuff on-shell. Then we
  // return the outgoing propagator, after setting the ordering for
  // the subsequent emissions, and we're done.
  if ( !unresolved ) {
    if ( mode >= 0 ) prop = mother()->propagator(res(), mode);
    sibling()->setOnShell(mode);
    return setup(prop);
  }

  // OK, so it was unresolved. Get the incoming propagator with the
  // previous resolution scale. If it was not forced, then let it go.
  prop = mother()->propagator(r2, mode);
  if ( !forced ) return prop;

  // However, if we have to force it we need to check that it really works.
  if ( !setEmissionMomenta(prop.p, true) ) {
    // If it doesn't work, just give up the whole thing and fail the
    // whole propagator chain.
    prop.fail = true;
    return prop;
  }

  // If it does work we set things on-shell, return the propagator and
  // we're done.
  sibling()->setOnShell(mode);
  prop.p -= sibling()->momentum();
  return prop;

}

ShadowParton::Propagator ShadowParton::setup(Propagator & prop) const {
  if ( colourSibling() ) {
    prop.acopos = sibling()->plus();
    prop.aconeg = sibling()->minus();
    prop.acoptp = prop.p.perp2();
  } else {    
    prop.colpos = sibling()->plus();
    prop.colneg = sibling()->minus();
    prop.colptp = prop.p.perp2();
  }
  prop.p -= sibling()->momentum();
  return prop;
}

tSPartonPtr ShadowParton::findFirstOnShell() {
  if ( original()->onShell() ) return this;
  if ( !next() ) return tSPartonPtr();
  if ( next()->colourSibling() ) {
    tSPartonPtr ch = child()->findFirstOnShell();
    if ( ch ) return ch;
  }
  return next()->findFirstOnShell();
}
  
tSPartonPtr ShadowParton::findSecondOnShell() {
  if ( original()->onShell() ) return this;
  if ( !next() ) return tSPartonPtr();
  if ( !next()->colourSibling() ) {
    tSPartonPtr ch = child()->findSecondOnShell();
    if ( ch ) return ch;
  }
  return next()->findSecondOnShell();
}

Ariadne5::ClonePtr ShadowParton::clone() const {
  return new_ptr(*this);
}

void ShadowParton::mirror(double yf) {
  y(2.0*yf - y());
  swap(thePlus, theMinus);
  if ( next() ) next()->mirror(yf);
  if ( child() ) child()->mirror(yf);
}

void ShadowParton::translate(const ImpactParameters & b) {
  thePT = b.rotatePT(thePT);
  if ( next() ) next()->translate(b);
  if ( child() ) child()->translate(b);
}

void ShadowParton::rebind(const TranslationMap & trans) {
  thePrevious = trans.translate(thePrevious);
  theParent = trans.translate(theParent);
  theSibling = trans.translate(theSibling);
  theNext = trans.translate(theNext);
  theChild = trans.translate(theChild);
}


DipoleState& ShadowParton::dipoleState() const {
  return original()->dipoleState();
}


double ShadowParton::pTScale() const {
  return Current<DipoleEventHandler>()->emitter().pTScale();
}


Energy ShadowParton::mass() const {
  if ( theMass < ZERO )
    theMass = CurrentGenerator::current().getParticleData(flavour())->mass();
  return theMass;
}

void ShadowParton::debugme() {
  // cout << "data for ShadowParton " << this << endl;
  // cout << "thePosition1: " << thePosition.first*GeV
  //      << ", thePosition2: " << thePosition.second*GeV
  //      << ", thePlus: " << thePlus/GeV
  //      << ", thePT: " << thePT.pt()/GeV
  //      << ", theMinus: " << theMinus/GeV
  //      << ", theY: " << theY
  //      << ", theFlavour: " << theFlavour
  //      << ", hasInteracted: " << hasInteracted
  //      << ", isOnShell: " << isOnShell
  //      << ", isValence: " << isValence
  //      << ", theMass: " << theMass/GeV << endl;
  memememe = true;
  valenceMother()->debugTree("");
  memememe = false;
  
}

void ShadowParton::checkMomentum(Sum20Momentum & sum20,
				 const ImpactParameters * b) const {
  if ( onShell() ) {
    if ( b ) sum20 -= lightCone(minus(), plus(), b->rotatePT(pT()));
    else  sum20 -= momentum();
  }
  if ( next() ) next()->checkMomentum(sum20, b);
  if ( child() ) child()->checkMomentum(sum20, b);
}

void ShadowParton::debugTree(string indent) {
  cerr << indent << (memememe? "==": "--")
       << (valence()? "v": "-") << (locked()? "%": "-");
  if ( interactionRoot() ) {
    set<int>::iterator i = theInteractionRoots.begin();
    cerr << "{" << *i;
    while ( ++i != theInteractionRoots.end() )
      cerr << "," << *i;
    cerr << "}";
  }
  cerr << (interacted()? "+": "-");
  if ( interacting() ) {
    set<int>::iterator i = theInteractions.begin();
    cerr << "[" << *i;
    while ( ++i != theInteractions.end() )
      cerr << "," << *i;
    cerr << "]";
  }
  cerr << (onShell()? "*": "-") << (memememe? "=>": "->");
  if ( onShell() )
    cerr << "[" << pT().x()/GeV << ", "
	 << pT().y()/GeV << ", "
	 << plus()/GeV << "]";
  else if ( interacted() )
    cerr << "(" << pT().x()/GeV << ", "
	 << pT().y()/GeV << ", "
	 << plus()/GeV << ")";
  cerr << endl;
  if ( indent.size() > 1 && indent[indent.size() - 1] == '\\' )
    indent[indent.size() - 1] = ' ';
  if ( next() && child() ) {
    if ( child()->interacted() ) next()->debugTree(indent + " +");
    else next()->debugTree(indent + " |");
  }
  if ( next() && !child() ) next()->debugTree(indent + "  ");
  if ( child() ) child()->debugTree(indent + " \\");
}
/*

----->
 |---->
 | |---->
 | \---->
 |   |---->
 |   \---->
 \---->
 */


void ShadowParton::persistentOutput(PersistentOStream & os) const {
  os << theOriginal << ounit(thePlus, GeV)
     << ounit(thePT, GeV) << ounit(thePT0, GeV) << ounit(theMass, GeV)
     << ounit(theMinus, GeV) << theY << theY0 << theFlavour
     << thePrevious << theParent << theSibling
     << theNext << theChild
     << ounit(theEmissionFactor, 1.0/GeV2) << ounit(theRes, 1.0/GeV2)
     << hasInteracted << isOnShell << isValence << hasColourSibling;
}

void ShadowParton::persistentInput(PersistentIStream & is, int) {
  is >> theOriginal >> iunit(thePlus, GeV)
     >> iunit(thePT, GeV) >> iunit(thePT0, GeV) >> iunit(theMass, GeV)
     >> iunit(theMinus, GeV) >> theY >> theY0 >> theFlavour
     >> thePrevious >> theParent >> theSibling
     >> theNext >> theChild
     >> iunit(theEmissionFactor, 1.0/GeV2) >> iunit(theRes, 1.0/GeV2)
     >> hasInteracted >> isOnShell >> isValence >> hasColourSibling;
}

DescribeClass<ShadowParton,Ariadne5::CloneBase>
describeDIPSYShadowParton("DIPSY::ShadowParton", "libAriadne5.so libDIPSY.so");
// Definition of the static class description member.

void ShadowParton::Init() {}

