// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Parton class.
//

#include "Parton.h"
#include "DipoleState.h"
#include "DipoleEventHandler.h"
#include "Dipole.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Utilities/Throw.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace DIPSY;

Parton::~Parton() {}

Ariadne5::ClonePtr Parton::clone() const {
  return new_ptr(*this);
}

void Parton::rebind(const TranslationMap & trans) {
  theParents.first = trans.translate(theParents.first);
  theParents.second = trans.translate(theParents.second);
  theDipoles.first = trans.translate(theDipoles.first);
  theDipoles.second = trans.translate(theDipoles.second);
  set<tPartonPtr> kids;
  kids.swap(theChildren);
  trans.translate(inserter(theChildren), kids.begin(), kids.end());
}

void Parton::interact() {
  if ( interacted() ) return;
  hasInteracted = true;
  if ( parents().first ) parents().first->interact();
  if ( parents().second ) parents().second->interact();
}

bool Parton::inInteractingLoop() const {
  //search in one direction
  DipolePtr d = dipoles().first;
  while ( d && d != dipoles().second && !(d->interacted()) ) {
    if ( d->neighbors().first )
      d = d->neighbors().first;
    else
      break;
  }
  if ( !d ) {    //if stopped at a quark, look also in other direction
    d = dipoles().second;
    while ( d && d != dipoles().second && !(d->interacted()) ) {
      if ( d->neighbors().first )
	d = d->neighbors().first;
      else
	break;
    }
  }
  if ( d && d->interacted() )
    return true;
  else
    return false;
}

int Parton::nOnShellInChain() const {
  int ret = 0;
  DipolePtr d = dipoles().first;
  if ( onShell() ) ret++;
  while ( d && d != dipoles().second ) {
    if ( d->partons().first->onShell() ) ret++;
    d = d->neighbors().first;
  }
  if ( !d ) {
    d = dipoles().second;
    while ( d ) {
      if ( d->partons().second->valence() ) ret++;
      d = d->neighbors().second;
    }
  }
  return ret;
}

DipoleState& Parton::dipoleState() const {
  if ( dipoles().first ) return dipoles().first->dipoleState();
  else return dipoles().second->dipoleState();
}

bool Parton::valenceParton() const {
  DipoleState & state = dipoleState();
  for ( int i = 0, N = state.initialDipoles().size(); i < N; i++ ) {
    if ( state.initialDipoles()[i]->partons().first == this || 
         state.initialDipoles()[i]->partons().second == this )
      return true;
  }
  return false;
}

bool Parton::inValenceChain() const {
  //search in one direction
  DipolePtr d = dipoles().first;
  if ( valence() ) {
    if ( !d )
      cout << "found valence quark" << endl;
    return true;
  }
  while ( d && d != dipoles().second ) {
    if ( d->partons().first->valence() ) {
      return true;
    }
    d = d->neighbors().first;
  }
  if ( !d ) {    //if stopped at a quark, look also in other direction
    cout << "found non-valence quark at x =" << position().x()*GeV 
	 << ", y = " << position().y()*GeV << endl;
    DipoleState & state = dipoleState();
    for ( int i = 0, N = state.initialDipoles().size(); i < N; i++ ) {
      cout << "valence partons at x = "
	   << state.initialDipoles()[i]->partons().first->position().x()*GeV
	   << ", y = " << state.initialDipoles()[i]->partons().first->position().y()*GeV << endl;
      cout << "valence partons at x = "
	   << state.initialDipoles()[i]->partons().second->position().x()*GeV
	   << ", y = "
	   << state.initialDipoles()[i]->partons().second->position().y()*GeV << endl;
    }
    d = dipoles().second;
    while ( d && d != dipoles().second ) {
      if ( d->partons().second->valence() ) {
      return true;
      }
    d = d->neighbors().second;
    }
  }
  return false;
}

bool Parton::swingedEmission() const {
  if ( !parents().first || !parents().second ) return false;
  else if ( parents().first->parents().second == parents().second ) return false;
  else if ( parents().second->parents().first == parents().first ) return false;
  if ( Debug::level > 5 ) cout << "parton at " << oY() << " is a swinged emission" << endl;
  return true;
}

bool Parton::absorbed() const {
  return dipoles().first->children().first || dipoles().first->children().second ||
    dipoles().second->children().first || dipoles().second->children().second;
}

double Parton::pTScale() const {
  return dipoleState().handler().emitter().pTScale();
}

void Parton::updateMomentum() {
  if ( dipoles().first && dipoles().second ) {
    tPartonPtr p1 = dipoles().first->partons().first;
    tPartonPtr p2 = dipoles().second->partons().second;
    pT( pTScale()*((position() - p1->position())/(position() - p1->position()).pt2() +
	     (position() - p2->position())/(position() - p2->position()).pt2()) );
  }
  else if ( dipoles().first ) {
    tPartonPtr p1 = dipoles().first->partons().first;
    pT( pTScale()*(position() - p1->position())/(position() - p1->position()).pt2());
  }
  else if ( dipoles().second ) {
    tPartonPtr p2 = dipoles().second->partons().second;
    pT( pTScale()*(position() - p2->position())/(position() - p2->position()).pt2() );
  }
  else
    cout << "Parton without any dipoles?? :o" << endl;
  if ( isRightMoving ) {
    y( log(pT().pt()/plus()) );
    minus( pT().pt()*exp(y()) );
  }
  else {
    y( log(minus()/pT().pt()) );
    plus( pT().pt()*exp(-y()) );
  }
}

TransverseMomentum Parton::recoil(tPartonPtr p) const {
  return pTScale()*(position() - p->position())/dist2(*p);
}

Energy Parton::mass() const {
  if ( theMass < ZERO )
    theMass = CurrentGenerator::current().getParticleData(flavour())->mass();
  return theMass;
}

PPtr Parton::produceParticle() const {
  return dipoleState().getParticle(this);
}

void Parton::coutData() {
  cout << "data for Parton " << this << endl;
  cout << "thePosition1: " << thePosition.first*GeV
       << ", thePosition2: " << thePosition.second*GeV
       << ", thePlus: " << thePlus/GeV
       << ", theOriginalY: " << theOriginalY
       << ", thePT: " << thePT.pt()/GeV
       << ", theValencePT: " << theValencePT.pt()/GeV
       << ", theValencePlus: " << theValencePlus/GeV
       << ", theMinus: " << theMinus/GeV
       << ", theY: " << theY
       << ", isRightMoving: " << isRightMoving
       << ", theFlavour: " << theFlavour
       << ", theParents1: " << theParents.first
       << ", theParents2: " << theParents.second
       << ", theChildren: ";
  for( set<tPartonPtr>::iterator it = theChildren.begin();
       it != theChildren.end();it++)
    cout << *it << ", ";
  cout << "theDipoles1: " << theDipoles.first
       << ", theDipoles2: " << theDipoles.second
       << ", hasInteracted: " << hasInteracted
       << ", isOrdered: " << isOrdered
       << ", isOnShell: " << isOnShell
       << ", isValence: " << isValence
       << ", theNumber: " << theNumber
       << ", theMass: " << theMass/GeV << endl;
}

void Parton::persistentOutput(PersistentOStream & os) const {
  os << ounit(thePosition, InvGeV) << ounit(thePlus, GeV) << theOriginalY
     << ounit(thePT, GeV) << ounit(theValencePT, GeV)
     << ounit(theValencePlus, GeV) << ounit(theMinus, GeV) << theY
     << isRightMoving << theFlavour << theParents << theDipoles
     << hasInteracted << isOrdered << isOnShell << isValence
     << ounit(theMass, GeV);
}

void Parton::persistentInput(PersistentIStream & is, int) {
  is >> iunit(thePosition, InvGeV) >> iunit(thePlus, GeV) >> theOriginalY
     >> iunit(thePT, GeV) >> iunit(theValencePT, GeV)
     >> iunit(theValencePlus, GeV) >> iunit(theMinus, GeV) >> theY
     >> isRightMoving >> theFlavour >> theParents >> theDipoles
     >> hasInteracted >> isOrdered >> isOnShell >> isValence
     >> iunit(theMass, GeV);
}

DescribeClass<Parton,Ariadne5::CloneBase>
describeDIPSYParton("DIPSY::Parton", "libAriadne5.so libDIPSY.so");
// Definition of the static class description member.

void Parton::Init() {}

