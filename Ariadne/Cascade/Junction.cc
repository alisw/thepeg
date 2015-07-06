// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Junction class.
//

#include "Junction.h"
#include "QCDDipole.h"
#include "DipoleState.h"

#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

Junction::Junction(): Parton(true) {}

Junction::~Junction() {}

ClonePtr Junction::clone() const {
  return new_ptr(*this);
}

tPPtr Junction::produceParticle(const LorentzRotation & r) {
  return tPPtr();
}

void Junction::fillReferences(CloneSet & cset) const {
  Parton::fillReferences(cset);
}

void Junction::rebind(const TranslationMap & trans) {
  Parton::rebind(trans);
  p1 = trans.translate(p1);
  p2 = trans.translate(p2);
  p3 = trans.translate(p3);
  d1 = trans.translate(d1);
  d2 = trans.translate(d2);
  d3 = trans.translate(d3);
}

bool Junction::replace(tQCDPtr oldd, tQCDPtr newd) {
  if ( oldd == d1 ) {
    d1 = newd;
    p1 = newd->iPart() == this? newd->oPart(): newd->iPart();
  }
  else if ( oldd == d2 ) {
    d2 = newd;
    p2 = newd->iPart() == this? newd->oPart(): newd->iPart();
  }
  else if ( oldd == d3 ) {
    d3 = newd;
    p3 = newd->iPart() == this? newd->oPart(): newd->iPart();
  }
  else return false;
  return true;
}

pair<tcQCDPtr,tcQCDPtr> Junction::getNeighbours(tcQCDPtr d) const {
  pair<tQCDPtr,tQCDPtr> ret(d1, d2);
  if ( d == d1 ) ret.first = d3;
  if ( d == d2 ) ret.second = d3;
  return ret;
}


tParPtr Junction::getNeighbour(tcParPtr p) const {
  tParPtr ret;
  if ( p == p1 ) ret = p2;
  if ( p == p2 ) ret = p3;
  if ( p == p3 ) ret = p1;
  if ( tJunctionPtr j = dynamic_ptr_cast<tJunctionPtr>(ret) )
    return j->getNeighbour(this);
  return ret;
}

bool Junction::source() const {
  return d3 && d3->oPart() == this;
}

bool Junction::sink() const {
  return d3 && d3->iPart() == this;
}

tParPtr Junction::getRandomRecoiler(tcQCDPtr d) const {
  pair<tcQCDPtr,tcQCDPtr> nbs = getNeighbours(d);
  tcQCDPtr sel = UseRandom::rndbool()? nbs.first: nbs.second;
  tParPtr p = ( sel->iPart() == this? sel->oPart(): sel->iPart() );
  if (  tJunctionPtr j = dynamic_ptr_cast<tJunctionPtr>(p) )
    return j->getRandomRecoiler(sel);
  return p;
}

bool Junction::touchedNeighbours(tcQCDPtr d) const {
  if ( d1 != d && touchedBranch(d1) ) return true;
  if ( d2 != d && touchedBranch(d2) ) return true;
  if ( d3 != d && touchedBranch(d3) ) return true;
  return false;
}

bool Junction::touchedBranch(tcQCDPtr d) const {
  if ( d->touched() ) return true;
  tcParPtr p = (d->iPart() != this? d->iPart(): d->oPart());
  if ( p->touched() || p->touchedNeighbours(d) ) return true;
  return false;
}

void Junction::persistentOutput(PersistentOStream & os) const {
  os << d1 << d2 << d3 << p1 << p2 << p3;
}

void Junction::persistentInput(PersistentIStream & is, int) {
  is >> d1 >> d2 >> d3 >> p1 >> p2 >> p3;
}

void Junction::debugme() const {
  CascadeBase::debugme();
  cerr << "J" << setw(3) << state()->index(this);
  if ( source() ) cerr << " (source)";
  else if ( sink() ) cerr << " (sink)";
  else cerr << " (incomplete)";
  if ( d1 ) cerr << setw(3) << state()->index(d1);
  if ( d2 ) cerr << setw(3) << state()->index(d2);
  if ( d3 ) cerr << setw(3) << state()->index(d3);
}

// The following static variable is needed for the type description
// system in ThePEG.
DescribeClass<Junction,Parton>
  describeAriadne5Junction("Ariadne5::Junction", "libAriadne5.so");

void Junction::Init() {}

