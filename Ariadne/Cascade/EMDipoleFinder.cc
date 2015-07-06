// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EMDipoleFinder class.
//

#include "EMDipoleFinder.h"
#include "EMDipole.h"
#include "Resonance.h"
#include "ResonanceParton.h"
#include "DipoleState.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/MaxCmp.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

EMDipoleFinder::EMDipoleFinder() {}

EMDipoleFinder::~EMDipoleFinder() {}

IBPtr EMDipoleFinder::clone() const {
  return new_ptr(*this);
}

IBPtr EMDipoleFinder::fullclone() const {
  return new_ptr(*this);
}

vector<tEMDipPtr> EMDipoleFinder::findDipoles(DipoleState & state) const {
  // *** ATTENTION *** Is this really reasonable?
  vector<tEMDipPtr> ret;
  list< pair<tParPtr,int> > pos;
  list< pair<tParPtr,int> > neg;
  MaxCmp<int> subsysmax;
  typedef list< pair<tParPtr,int> >::iterator ListIter;
  for ( set<tParPtr>::const_iterator it = state.finalState().begin();
	it != state.finalState().end(); ++it ) {
    tParPtr p = *it;
    int ss = 0;
    if ( tResParPtr rp = dynamic_ptr_cast<tResParPtr>(p) )
      ss = rp->resonance()->decaySystem();
    subsysmax(ss);
    for ( int i = 0; i < p->data().iCharge(); ++i )
      pos.insert(pos.end(), make_pair(p, ss));
    for ( int i = 0; i > p->data().iCharge(); --i )
      neg.insert(neg.end(), make_pair(p, ss));
  }

  while ( pos.size() && neg.size() ) {
    // While there is still charged partons left first find the pair
    // with the smallest invariant mass.
    MinCmp<Energy2,pair<ListIter,ListIter> > smin;
    // First consider only pairs within the same system.
    for ( ListIter itp = pos.begin(); itp != pos.end(); ++itp )
      for ( ListIter itn = neg.begin(); itn != neg.end(); ++itn )
	if ( itp->second == itn->second )
	  smin((itp->first->momentum() + itn->first->momentum()).m2(),
	       make_pair(itp, itn));
    // If no pair in the same system is found consider any pair.
    if ( !smin.index().first->first )
      for ( ListIter itp = pos.begin(); itp != pos.end(); ++itp )
	for ( ListIter itn = neg.begin(); itn != neg.end(); ++itn )
	  smin((itp->first->momentum() + itn->first->momentum()).m2(),
	       make_pair(itp, itn));
    tEMDipPtr dip = state.create<EMDipole>();
    dip->iPart(smin.index().first->first);
    dip->oPart(smin.index().second->first);
    dip->system(min(smin.index().first->second, smin.index().second->second));
    pos.erase(smin.index().first);
    pos.erase(smin.index().second);
    ret.push_back(dip);
  }

  // Now continue with all left-over charged partons
  for ( ListIter itp = pos.begin(); itp != pos.end(); ++itp ) {
    tEMDipPtr dip = state.create<EMDipole>();
    dip->iPart(itp->first);
    ret.push_back(dip);
  }
  for ( ListIter itn = neg.begin(); itn != neg.end(); ++itn ) {
    tEMDipPtr dip = state.create<EMDipole>();
    dip->oPart(itn->first);
    ret.push_back(dip);
  }

  return ret;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void EMDipoleFinder::persistentOutput(PersistentOStream & os) const {}

void EMDipoleFinder::persistentInput(PersistentIStream & is, int) {}


// The following static variable is needed for the type description
// system in ThePEG.
DescribeClass<EMDipoleFinder,HandlerBase>
  describeAriadne5EMDipoleFinder("Ariadne5::EMDipoleFinder", "libAriadne5.so");

void EMDipoleFinder::Init() {

  static ClassDocumentation<EMDipoleFinder> documentation
    ("The EMDipoleFinder class and its sub-classes are responsible for "
     "identifying and introducing of electro-magnetic dipoles in the setup "
     "phase of a given DipoleState.");

}

