// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LargePTDipoleAbsorber class.
//
// WARNING!
// This class is not in use any longer, and is outdated.
// Can be removed, or need to be completely revised if
// taken back in use.
//

#include <algorithm>

#include "LargePTDipoleAbsorber.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Parton.h"
#include "Dipole.h"
#include "DipoleState.h"
#include "ThePEG/Repository/UseRandom.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

LargePTDipoleAbsorber::LargePTDipoleAbsorber() {}

LargePTDipoleAbsorber::~LargePTDipoleAbsorber() {}

IBPtr LargePTDipoleAbsorber::clone() const {
  return new_ptr(*this);
}

IBPtr LargePTDipoleAbsorber::fullclone() const {
  return new_ptr(*this);
}

//checks if the parton is part of an interacting dipole or not
bool interacted ( tPartonPtr p ) {
  return p->dipoles().first->interacted() || p->dipoles().second->interacted();
}

//checks that the parton hasnt been absorbed yet.
bool absorbed ( tPartonPtr p ) {
  return p->dipoles().first->partons().second != p || 
         p->dipoles().second->partons().first != p;
}

void LargePTDipoleAbsorber::hiddenRecoil(DipoleState & state) const {
  cout << "hidden recoil not implemented for large PT absorber" << endl;
}

void LargePTDipoleAbsorber::
absorbParton( tPartonPtr p, DipoleState & state) const {

  if ( p == p->dipoles().first->neighbors().first->partons().first ) {
    if ( interacted(p) ) {
      return;
    }
    state.sortDipoles();
    if( !(state.forceSwing(p, 0.0, 0.0)) ) {
      cout << "failed forceswing, will not absorb" << endl;
      return;   //TODO: move to this class. think through.
    }
  }

  //Find the best parton to merge with.
  list<PartonPtr> partons = state.getPartons();
  PartonPtr absP;
  double bestPTLeft = 1.0;
  for ( list<PartonPtr>::const_iterator it = partons.begin() ; it != partons.end(); it++ ) {
    if ( *it == p ) continue;
    // Here is physics! May be done differently!!
    double pTLeft = ((*it)->pT() + p->pT()).pt()/((*it)->pT().pt() + p->pT().pt());
    if ( pTLeft < bestPTLeft ) {
      if (  !(interacted(p)) || !(interacted(*it)) ) {
        bestPTLeft = pTLeft;
        absP = *it;
      }
    }
  }
  if ( !absP ) {
//     cout << "couldnt find a parton to merge with, will not absorb" << endl;
    return;
  }

  //if p is part of an interacting dipole, then let absP be absorbed instead.
  if ( interacted(p) ) {
    PartonPtr dummy = p;
    p = absP;
    absP = dummy;
  }

  //Move pT to the absorbing parton, conserving pT and p+
  absP->pT( absP->pT() + p->pT() );
  absP->plus( absP->plus() + p->plus() );
  absP->y( log( absP->pT().pt()/absP->plus() ) );
  absP->minus( absP->pT().pt()*exp(absP->y()) );

  //recconect colour flow
  DipolePtr d1 = p->dipoles().first;
  DipolePtr d2 = p->dipoles().second;
  DipolePtr d3 = d2->neighbors().second;
  PartonPtr p2 = d2->partons().second;

  d2->children().second = d1;
  d1->neighbors( make_pair(d1->neighbors().first, d3) );
  d3->neighbors( make_pair(d1, d3->neighbors().second) );
  d1->partons( make_pair(d1->partons().first, p2) );
  p2->dipoles( make_pair(d1, p2->dipoles().second) );

  //chose colour of the new dipole from the bigger of the two old ones.
  if ( d1->size()/(d1->size()+d2->size()) < UseRandom::rnd() )
    d1->colour(d2->colour());

}

void LargePTDipoleAbsorber::
absorbDipole( tDipolePtr d, DipoleState & state ) const{
  // IMPLEMENT!!!
}

bool pTSort( tcPartonPtr p1, tcPartonPtr p2 ) {
  return p1->pT().pt() > p2->pT().pt();
}

//This algorithm often pairs up pT incorectly, leaving pT too large soemtimes.
void LargePTDipoleAbsorber::reabsorb(DipoleState & state) const {

  absorbVirtualPartons(state);

  //Find the largest pT of an interacting dipole. Maybe this scale should be
  //calculated in another way.
  Energy pTScale = 0.0*GeV;
  list<DipolePtr> dipoles = state.getDipoles();
  for ( list<DipolePtr>::const_iterator it = dipoles.begin() ; it != dipoles.end(); it++ ) {
    if ( (*it)->interacted() && 1.0/(*it)->size() > pTScale )
      pTScale = 1.0/(*it)->size(); //or should it be a factor 2 here?? consistency...
  }

  //Go through the partons and see if they should be absorbed or not.
  list<PartonPtr> partons = state.getPartons();
  while ( !(partons.empty()) ) {
    list<PartonPtr>::iterator it = max_element( partons.begin(), 
						partons.end(), pTSort );
    if ( sqr(pTScale/(*it)->pT().pt()) < UseRandom::rnd() && !(absorbed(*it)) )
      absorbParton( *it, state);
    partons.erase( it );
  }
}

void LargePTDipoleAbsorber::absorbVirtualPartons( DipoleState & state ) const {
  list<PartonPtr> virtualPartons = state.virtualPartons();
  while ( !(virtualPartons.empty()) ) {
    list<PartonPtr>::iterator it = max_element( virtualPartons.begin(), 
						virtualPartons.end(), pTSort );
    absorbParton( *it, state);
    virtualPartons.erase( it );
  }
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void LargePTDipoleAbsorber::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void LargePTDipoleAbsorber::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<LargePTDipoleAbsorber,DIPSY::DipoleAbsorber>
  describeDIPSYLargePTDipoleAbsorber("DIPSY::LargePTDipoleAbsorber",
				     "LargePTDipoleAbsorber.so");


void LargePTDipoleAbsorber::Init() {

  static ClassDocumentation<LargePTDipoleAbsorber> documentation
    ("The LargePTDipoleAbsorber class implements a model describing how "
     "non-interacted dipoles will be reabsorbed. The model is based on "
     "the removal of dipoles with partons with large and opposite "
     "transverse momenta.");

}

