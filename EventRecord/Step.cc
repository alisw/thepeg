// -*- C++ -*-
//
// Step.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Step class.
//

#include "Step.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/algorithm.h"
#include <iostream>
#include <iomanip>

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
#include "Step.tcc"
#endif

using namespace ThePEG;

Step::Step(const Step & s)
  : Base(s), 
    theParticles(s.theParticles), theIntermediates(s.theIntermediates),
    allParticles(s.allParticles), theCollision(s.theCollision),
    theHandler(s.theHandler) {}

Step::~Step() {
  for ( ParticleSet::iterator it = allParticles.begin();
	it != allParticles.end(); ++it )
    if ( (**it).hasRep() && (**it).birthStep() == this )
      (**it).rep().theBirthStep = tStepPtr();
  theParticles.clear();
  theIntermediates.clear();
  theSubProcesses.clear();
  allParticles.clear();
  theCollision = tCollPtr();
  theHandler = tcEventBasePtr();
}

StepPtr Step::clone() const {
  return ptr_new<StepPtr>(*this);
}

void Step::rebind(const EventTranslationMap & trans) {
  ParticleSet::const_iterator pit;
  allParticles.clear();
  for ( pit = theParticles.begin(); pit != theParticles.end(); ++pit )
    allParticles.insert(trans.translate(*pit));
  theParticles.swap(allParticles);
  allParticles.clear();
  for ( pit = theIntermediates.begin();	pit != theIntermediates.end(); ++pit )
    allParticles.insert(trans.translate(*pit));
  theIntermediates = allParticles;
  allParticles.insert(theParticles.begin(), theParticles.end());
  for ( pit = allParticles.begin(); pit != allParticles.end(); ++pit )
    (**pit).rebind(trans);
}

void Step::addParticle(tPPtr p) {
  if ( !p->birthStep() ) p->rep().theBirthStep = this;
  theParticles.insert(p);
  allParticles.insert(p);
  if ( collision() ) collision()->addParticle(p);
}

void Step::addSubProcess(tSubProPtr sp) {
  if (( !member(allParticles, sp->incoming().first) && collision() &&
	sp->incoming().first != collision()->incoming().first) &&
      std::find(sp->incoming().first->parents().begin(),
		sp->incoming().first->parents().end(),
		collision()->incoming().first)==
      sp->incoming().first->parents().end()) {
      collision()->incoming().first->
	rep().theChildren.push_back(sp->incoming().first);
      sp->incoming().first->rep().theParents.push_back(collision()->incoming().first);
  }
  if (( !member(allParticles, sp->incoming().second) && collision() &&
	sp->incoming().second != collision()->incoming().second) &&
      std::find(sp->incoming().second->parents().begin(),
		sp->incoming().second->parents().end(),collision()->incoming().second)==
      sp->incoming().second->parents().end()) {
    collision()->incoming().second->
      rep().theChildren.push_back(sp->incoming().second);
    sp->incoming().second->rep().theParents.push_back(collision()->incoming().second);
  }
  if ( !sp->incoming().first->birthStep() )
    sp->incoming().first->rep().theBirthStep = this;
  if ( !sp->incoming().second->birthStep() )
    sp->incoming().second->rep().theBirthStep = this;
  addIntermediate(sp->incoming().first);
  addIntermediate(sp->incoming().second);
  addIntermediates(sp->intermediates().begin(), sp->intermediates().end());
  addParticles(sp->outgoing().begin(), sp->outgoing().end());
  theSubProcesses.push_back(sp);
  if ( collision() ) collision()->addSubProcess(sp);
}

void Step::removeSubProcess(tSubProPtr sp) {
  SubProcessVector::iterator sit = ThePEG::find(theSubProcesses, sp);
  if ( sit == theSubProcesses.end() ) return;
  for ( int i = 0, N = sp->outgoing().size(); i < N; ++i )
    removeParticle(sp->outgoing()[i]);
  for ( int i = 0, N = sp->intermediates().size(); i < N; ++i )
    removeParticle(sp->intermediates()[i]);
  removeParticle(sp->incoming().first);
  removeParticle(sp->incoming().second);
  theSubProcesses.erase(sit);
  if ( collision() ) collision()->removeSubProcess(sp);
}

void Step::addIntermediate(tPPtr p) {
  theIntermediates.insert(p);
  ParticleSet::iterator pit = theParticles.find(p);
  if ( pit != theParticles.end() ) theParticles.erase(pit);
  else {
    if ( !p->birthStep() ) p->rep().theBirthStep = this;
    allParticles.insert(p);
    if ( collision() ) collision()->addParticle(p);
  }
}

void Step::
insertIntermediate(tPPtr p, tPPtr parent, tPPtr child) {
  if ( !p->birthStep() ) p->rep().theBirthStep = this;
  addIntermediate(p);
  parent->removeChild(child);
  child->removeParent(parent);
  if ( parent ) {
    parent->rep().theChildren.push_back(p);
    p->rep().theParents.push_back(parent);
  }
  if ( child ) {
    p->rep().theChildren.push_back(child);
    child->rep().theParents.push_back(p);
  }
}

void Step::removeEntry(tPPtr p) {
  ParticleSet::iterator it = allParticles.find(p);
  if ( it == allParticles.end() ) return;
  allParticles.erase(it);
  it = theParticles.find(p);
  if ( it != theParticles.end() ) theParticles.erase(it);

  if ( p->previous() ) {
    it = theIntermediates.find(p->previous());
    if ( it != theIntermediates.end() ) theIntermediates.erase(it);
    theParticles.insert(p->previous());
    allParticles.insert(p->previous());
  }
  while ( !p->parents().empty() ) {
    PPtr parent = p->parents().back();
    p->removeParent(parent);
    parent->removeChild(p);
    if ( !parent->children().empty() ) continue;
    it = theIntermediates.find(parent);
    if ( it != theIntermediates.end() ) theIntermediates.erase(it);
    theParticles.insert(parent);
    allParticles.insert(parent);
  }
  if ( p->hasColourInfo() ) {
    if ( colourNeighbour(p) )
      colourNeighbour(p)->antiColourNeighbour(antiColourNeighbour(p));
    if ( antiColourNeighbour(p) )
      antiColourNeighbour(p)->colourNeighbour(colourNeighbour(p));
    if ( p->incomingColour() ) p->outgoingColour(tPPtr());
    if ( p->incomingAntiColour() ) p->outgoingAntiColour(tPPtr());
  }

  it = theIntermediates.find(p);
  if ( it != theIntermediates.end() ) theIntermediates.erase(it);

}

void Step::removeParticle(tPPtr p) {
  if ( p->next() ) removeParticle(p->next());
  while ( !p->children().empty() ) removeParticle(p->children().back());
  removeEntry(p);
}

bool Step::nullStep() const {
  for ( ParticleSet::const_iterator it = allParticles.begin();
	it != allParticles.end(); ++it )
    if ( (**it).birthStep() == this ) return false;
  return true;
}

tPPtr Step::copyParticle(tcPPtr pin) {
  PPtr cp;
  tPPtr p = const_ptr_cast<tPPtr>(pin);
  if ( !collision() ) return cp;
  ParticleSet::iterator pit = theParticles.find(p);
  if ( collision()->finalStep() != this || p->next() ||
       ! p->children().empty() ||  pit == theParticles.end() ) return cp;
  cp = p->clone();
  cp->rep().thePrevious = p;
  p->rep().theNext = cp;
  if ( p->hasColour() ) p->colourFlow(cp);
  if ( p->hasAntiColour() ) p->antiColourFlow(cp);
  cp->rep().theBirthStep = this;
  theParticles.erase(pit);
  if ( p->birthStep() == this ) theIntermediates.insert(p);
  addParticle(cp);
  return cp;
}

bool Step::setCopy(tcPPtr poldin, tPPtr pnew) {
  if ( poldin->id() != pnew->id() ) return false;
  tPPtr pold = const_ptr_cast<tPPtr>(poldin);
  pold->rep().theNext = pnew;
  pnew->rep().thePrevious = pold;
  theParticles.erase(pold);
  if ( pold->birthStep() == this ) theIntermediates.insert(pold);
  pnew->rep().theBirthStep = this;
  addParticle(pnew);
  return true;
}

tPPtr Step::insertCopy(tcPPtr pin) {
  PPtr cp;
  tPPtr p = const_ptr_cast<tPPtr>(pin);
  if ( !collision() ) return cp;
  if ( collision()->all().find(p) == collision()->all().end() ) return cp;
  cp = p->clone();
  cp->rep().theNext = p;
  cp->rep().theChildren.clear();
  if ( p->previous() ) {
    p->previous()->rep().theNext = cp;
    cp->rep().thePrevious = p->previous();
  } else {
    for ( int i = 0, N = p->parents().size(); i < N; ++i ) {
      tPPtr parent = p->parents()[i];
      for ( int j = 0, M = parent->children().size(); j < M; ++j )
	if ( parent->children()[j] == p ) parent->rep().theChildren[j] = cp;
    }
  }
  p->rep().theParents.clear();
  p->rep().thePrevious = cp;
  if ( p->hasColour() ) cp->colourFlow(p);
  if ( p->hasAntiColour() ) cp->antiColourFlow(p);
  cp->rep().theBirthStep = this;
  theIntermediates.insert(cp);
  return cp;
}


bool Step::removeDecayProduct(tcPPtr par, tPPtr child) {
  if ( !collision() ) return false;
  tPPtr parent = const_ptr_cast<tPPtr>(par->final());
  if ( collision()->all().find(parent) == collision()->all().end() )
    return false;
  if ( !par->hasRep() ) return false;
  PVector::iterator it = ThePEG::find(parent->rep().theChildren, child);
  if ( it == parent->rep().theChildren.end() ) return false;
  parent->rep().theChildren.erase(it);
  ParticleSet::iterator cit = theParticles.find(child);
  if ( cit != theParticles.end() ) {
    theParticles.erase(cit);
    if ( child->birthStep() == this ) theIntermediates.insert(child);
  }
  return true;
}
  
bool Step::addDecayProduct(tcPPtr par, tPPtr child, bool fixColour) {
  if ( !collision() ) return false;
  tPPtr parent = const_ptr_cast<tPPtr>(par->final());
  if ( collision()->finalStep() != this || parent->next() ) return false;
  ParticleSet::iterator pit = theParticles.find(parent);
  if ( pit != theParticles.end() ) {
    theParticles.erase(pit);
    if ( parent->birthStep() == this ) theIntermediates.insert(parent);
  } else {
    if ( parent != collision()->incoming().first &&
	 parent != collision()->incoming().second &&
	 parent->children().empty() ) return false;
  }
  parent->rep().theChildren.push_back(child);
  child->rep().theParents.push_back(parent);
  child->rep().theBirthStep = this;
  addParticle(child);
  if ( !fixColour || !parent->hasColourInfo() || !parent->coloured() ||
       !child->coloured() ) return true;
  if ( parent->hasColour() && child->hasColour() &&
       !parent->outgoingColour() && !child->colourLine() )
    parent->colourFlow(child);
  if ( parent->hasAntiColour() && child->hasAntiColour() &&
       !child->antiColourLine() ) {
    if ( parent->outgoingAntiColour() )
      parent->antiColourLine()->
	removeAntiColoured(parent->outgoingAntiColour());
    parent->antiColourFlow(child);
  }
  return true;
}

void Step::addDecayNoCheck(tPPtr parent, tPPtr child) {
  ParticleSet::iterator pit = theParticles.find(parent);
  if ( pit != theParticles.end() ) {
    theParticles.erase(pit);
    if ( parent->birthStep() == this ) theIntermediates.insert(parent);
  }
  child->rep().theBirthStep = this;
  addParticle(child);
}

void Step::addDecayProduct(tPPtr child) {
  for ( int i = 0, N = child->parents().size(); i < N; ++i ) {
    ParticleSet::iterator pit = theParticles.find(child->parents()[i]);
    if ( pit != theParticles.end() ) {
      theParticles.erase(pit);
      if ( child->parents()[i]->birthStep() == this )
	theIntermediates.insert(child->parents()[i]);
    }
  }
  child->rep().theBirthStep = this;
  addParticle(child);
}

void Step::fixColourFlow() {
  tParticleVector news;
  for ( ParticleSet::iterator pi = theParticles.begin();
	pi != theParticles.end(); ++pi )
    if ( (**pi).birthStep() == this ) news.push_back(*pi);
  for ( int i = 0, N = news.size(); i < N; ++i ) {
    tPPtr p = news[i];
    if ( p->hasColour() && !antiColourNeighbour(p) ) {
      tPPtr ng = p;
      while ( ( ng = ng->incomingColour() ) && !antiColourNeighbour(ng) ) {}
      if ( ng ) {
	ng = antiColourNeighbour(ng);
	if ( !ng->outgoingColour() ) ng = copyParticle(ng);
	while ( ng->outgoingColour() ) ng = ng->outgoingColour();
	p->antiColourConnect(ng);
      }
    }
    if ( p->hasAntiColour() && !colourNeighbour(p) ) {
      tPPtr ng = p;
      while ( ( ng = ng->incomingAntiColour() ) && !colourNeighbour(ng) ) {}
      if ( ng ) {
	ng = colourNeighbour(ng);
	if ( !ng->outgoingAntiColour() ) ng = copyParticle(ng);
	while ( ng->outgoingAntiColour() ) ng = ng->outgoingAntiColour();
	p->colourConnect(ng);
      }
    }
  }
}

tPPtr Step::antiColourNeighbour(tcPPtr p) const {
  return colourNeighbour(p, true);
}

tPPtr Step::colourNeighbour(tcPPtr p, bool anti) const {
  if ( !member(particles(), const_ptr_cast<tPPtr>(p)) ) return tPPtr();
  tColinePtr line = p->colourLine(!anti);
  if ( !line ) return tPPtr();
  for ( ParticleSet::const_iterator it = particles().begin();
	it != particles().end(); ++it )
    if ( (**it).hasColourLine(line, anti) ) return *it;
  return tPPtr();
}

vector<tPVector> Step::getSinglets(tParticleSet & left) {
  vector<tPVector> ret;
  while ( !left.empty() ) {
    tPPtr first = *left.begin();
    left.erase(first);
    if ( !first->hasColourInfo() || !first->coloured() ) continue;
    tPPtr last = first;
    tPPtr test;
    while ( ( test = last->antiColourNeighbour(left.begin(), left.end()) ) &&
	    test != first )
      last = test;
    while ( ( test = first->colourNeighbour(left.begin(), left.end()) ) &&
	    test != last )
      first = test;
    ret.push_back(tPVector());
    for ( ; first != last;
	  first = first->antiColourNeighbour(left.begin(), left.end()) ) {
      left.erase(first);
      ret.back().push_back(first);
    }
    left.erase(first);
    ret.back().push_back(first);
  }
  return ret;
}

ostream & ThePEG::operator<<(ostream & os, const Step & s) {
  if ( !s.intermediates().empty() ) os << "--- intermediates:" << endl;
  Particle::PrintParticles(os, s.intermediates(), &s);
  os << "--- final:" << endl;
  LorentzMomentum sum;
  Energy2 sumx = Energy2();
  Energy2 sumy = Energy2();
  Energy2 sumz = Energy2();
  Particle::PrintParticles(os, s.particles(), &s);
  for ( ParticleSet::const_iterator it = s.particles().begin();
	it != s.particles().end(); ++it ) {
    sum += (**it).momentum();
    sumx += sqr((**it).momentum().x());
    sumy += sqr((**it).momentum().y());
    sumz += sqr((**it).momentum().z());
  }
  os << string(78, '-') << endl  << "     Sum of momenta:        ";
  int oldprecision = os.precision();
  Energy sumx1
    = ( sqr(sum.x()) > Constants::epsilon*sumx ? 
	sum.x(): ZERO );
  Energy sumy1
    = ( sqr(sum.y()) > Constants::epsilon*sumy ? 
	sum.y(): ZERO );
  Energy sumz1
    = ( sqr(sum.z()) > Constants::epsilon*sumz ? 
	sum.z(): ZERO );
  os << setprecision(3) << setw(10) << sumx1/GeV << setw(10) << sumy1/GeV
     << setw(10) << sumz1/GeV << setw(10) << sum.e()/GeV
     << setw(10) << sum.m()/GeV << endl << setprecision(oldprecision);
  return os;
}

void Step::debugme() const {
  cerr << *this;
  EventRecordBase::debugme();
}

void Step::persistentOutput(PersistentOStream & os) const {
  os << theParticles << theIntermediates << theSubProcesses << allParticles
     << theCollision;
  EventConfig::putHandler(os, theHandler);
}

void Step::persistentInput(PersistentIStream & is, int) {
  is >> theParticles >> theIntermediates >> theSubProcesses >> allParticles
     >> theCollision;
  EventConfig::getHandler(is, theHandler);
}

ClassDescription<Step> Step::initStep;

void Step::Init() {}
