// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleAbsorber class.
//

#include "DipoleAbsorber.h"
#include "Ariadne/DIPSY/DipoleState.h"
#include "Ariadne/DIPSY/Parton.h"
#include "Ariadne/DIPSY/Dipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/Current.h"
#include "DipoleEventHandler.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

DipoleAbsorber::DipoleAbsorber() {}

DipoleAbsorber::~DipoleAbsorber() {}

bool rapSorty(tcPartonPtr part1, tcPartonPtr part2) {
  return (part1->y() >= part2->y());
}

void DipoleAbsorber::removeParton(tPartonPtr p) const {
  DipolePtr d1 = p->dipoles().first;
  DipolePtr d2 = p->dipoles().second;

  //if a quark (ie only one colour line out), then just remove
  //the parton and the connected dipole.
  if ( !d1 ) {
    PartonPtr g = d2->partons().second;
    if ( !g->dipoles().first ) return;
    d2->neighbors().second->firstNeighbor(d1);
    g->dipoles(make_pair(d1, d2->neighbors().second));
    d2->children().second = d2->neighbors().second;
    g->flavour(p->flavour());
    return;
  }
  if ( !d2 ) {
    PartonPtr g = d1->partons().first;
    if ( !g->dipoles().second ) return;
    d1->neighbors().first->secondNeighbor(d2);
    g->dipoles(make_pair(d1->neighbors().first, d2));
    d1->children().second = d1->neighbors().first;
    g->flavour(p->flavour());
    return;
  }

  PartonPtr p1 = d1->partons().first;
  PartonPtr p2 = d2->partons().second;

  if ( p1->dipoles().second->partons().second == p2 ) return;

  d2->children().second = d1;
  p2->dipoles(make_pair(d1,p2->dipoles().second));
  d1->partons(make_pair(p1,p2));
  d1->secondNeighbor(d2->neighbors().second);
  if ( d1->neighbors().second ) d1->neighbors().second->firstNeighbor(d1);

  if ( d1->size()/(d1->size()+d2->size()) < UseRandom::rnd() )
    d1->colour(d2->colour());
}

void DipoleAbsorber::absorbVirtualPartons( DipoleState & state ) const {
  list<PartonPtr> virtualPartons = state.virtualPartons();
  virtualPartons.sort(rapSorty);
  while ( !(virtualPartons.empty()) ) {
    virtualPartons.sort(rapSorty);
    absorbParton( virtualPartons.front(), state );
    virtualPartons.pop_front();
  }
}

bool DipoleAbsorber::swing( DipolePtr d1, DipolePtr d2 ) const {
  if ( d1->neighbors().second == d2 ) return false;
  if ( d2->neighbors().second == d1 ) return false;

  if (d1 == d2) Throw<Exception>()
		  << "lol, swinging with itself?? >_>" << Exception::abortnow;
  
  if (d2->children().second || d2->children().first)
    Throw<Exception>()
      << "swinging with someone thats emitted already" << Exception::abortnow;

  PartonPtr p11 = d1->partons().first;
  PartonPtr p12 = d1->partons().second;
  PartonPtr p21 = d2->partons().first;
  PartonPtr p22 = d2->partons().second;
  
  d1->partons(make_pair(p11, p22));
  d1->neighbors(make_pair(p11->dipoles().first, p22->dipoles().second));
  if ( d1->neighbors().second ) d1->neighbors().second->firstNeighbor(d1);
  
  d2->partons(make_pair(p21, p12));
  d2->neighbors(make_pair(p21->dipoles().first, p12->dipoles().second));
  if ( d2->neighbors().second ) d2->neighbors().second->firstNeighbor(d2);

  p12->dipoles(make_pair(d2, p12->dipoles().second));
  p22->dipoles(make_pair(d1, p22->dipoles().second));

  return true;

}

void DipoleAbsorber::isolateNonParticipating( DipoleState & state ) const {
  for( int i = 0; i < int(state.initialDipoles().size()); i++ ) {
    DipolePtr d = state.initialDipoles()[i];
    if ( !(d->participating()) ) {
      if ( d->partons().first->dipoles().second != d->partons().second->dipoles().first )
	swing( d->partons().first->dipoles().second, d->partons().second->dipoles().first );
      d->partons().first->dipoles().second->participating( false );
    }
  }
}

void DipoleAbsorber::swingLoop( DipolePtr d, DipoleState & state ) const {
  //make a list of the partons and dipoles in the loop.
  list<PartonPtr> partons;
  list<DipolePtr> dipoles;
  PartonPtr p = d->partons().second;
  partons.push_back(p);
  dipoles.push_back(d);
  bool forward = true;
  while ( p != d->partons().first && (forward || p->dipoles().first) ) {
    if ( !(p->dipoles().second) ) {
      cout << "found a quark in SmallDipoleAbsorber::absorbLoop!!!!!!!" << endl;
      state.diagnosis(true);
      p = d->partons().first;
      partons.push_back(p);
      forward = false;
      continue;
    }
    if ( forward ) {
      p = p->dipoles().second->partons().second;
      dipoles.push_back(p->dipoles().first);
    }
    else {
      p = p->dipoles().first->partons().first;
      dipoles.push_back(p->dipoles().second);
    }
    partons.push_back(p);
  }

  //find and perform a probable swing. Don't bother with colours...
  list<DipolePtr> otherDips = d->dipoleState().getDipoles();
  for ( list<DipolePtr>::iterator i = dipoles.begin(); i != dipoles.end(); i++ ) {
    for ( list<DipolePtr>::iterator j = otherDips.begin(); j != otherDips.end(); j++ ) {
      if ( *i == *j )
  	j = otherDips.erase(j)--;
    }
  }

  DipolePtr d1, d2;
  double maxy = 0.0;
  for ( list<DipolePtr>::iterator i = dipoles.begin(); i != dipoles.end(); i++ ) {
    for ( list<DipolePtr>::iterator j = otherDips.begin(); j != otherDips.end(); j++ ) {
      double R = log( UseRandom::rnd() );
      InvEnergy2 a = (*i)->partons().first->dist2(*(*i)->partons().second);
      InvEnergy2 b = (*j)->partons().first->dist2(*(*j)->partons().second);
      InvEnergy2 c = (*i)->partons().second->dist2(*(*j)->partons().first);
      InvEnergy2 d = (*i)->partons().first->dist2(*(*j)->partons().second);
      a = sqr(Current<DipoleEventHandler>()->rMax())/
  	(Current<DipoleEventHandler>()->alphaS(sqrt(a)))*
  	sqr(exp(sqrt(a)/Current<DipoleEventHandler>()->rMax()) - 1.0);
      b = sqr(Current<DipoleEventHandler>()->rMax())/
  	(Current<DipoleEventHandler>()->alphaS(sqrt(b)))*
  	sqr(exp(sqrt(b)/Current<DipoleEventHandler>()->rMax()) - 1.0);
      c = sqr(Current<DipoleEventHandler>()->rMax())/
  	(Current<DipoleEventHandler>()->alphaS(sqrt(c)))*
  	sqr(exp(sqrt(c)/Current<DipoleEventHandler>()->rMax()) - 1.0);
      d = sqr(Current<DipoleEventHandler>()->rMax())/
  	(Current<DipoleEventHandler>()->alphaS(sqrt(d)))*
  	sqr(exp(sqrt(d)/Current<DipoleEventHandler>()->rMax()) - 1.0);
      double A = c*d/(a*b);
      double y = R*A;
      if ( y < maxy || maxy == 0.0 ) {
  	maxy = y;
  	d1 = *i;
  	d2 = *j;
      }
    }
  }
  swing(d1, d2);
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).



// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeAbstractNoPIOClass<DipoleAbsorber,HandlerBase>
  describeDIPSYDipoleAbsorber("DIPSY::DipoleAbsorber", "libAriadne5.so libDIPSY.so");

void DipoleAbsorber::Init() {

  static ClassDocumentation<DipoleAbsorber> documentation
    ("DipoleAbsorber is a base class to be used for models describing how "
     "non-interacted dipoles will be reabsorbed.");

}

