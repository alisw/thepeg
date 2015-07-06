// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QCDDipoleFinder class.
//

#include "QCDDipoleFinder.h"
#include "QCDDipole.h"
#include "ResonanceParton.h"
#include "Junction.h"
#include "Resonance.h"
#include "DipoleState.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/Throw.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

QCDDipoleFinder::QCDDipoleFinder() {}

QCDDipoleFinder::~QCDDipoleFinder() {}

IBPtr QCDDipoleFinder::clone() const {
  return new_ptr(*this);
}

IBPtr QCDDipoleFinder::fullclone() const {
  return new_ptr(*this);
}

vector<tQCDPtr> QCDDipoleFinder::findDipoles(DipoleState & state) const {
  typedef multimap<tColinePtr,tParPtr>::iterator MapIt;
  vector<tQCDPtr> ret;
  multimap<tColinePtr,tParPtr> icol;
  multimap<tColinePtr,tParPtr> ocol;
  multimap<tColinePtr,tParPtr> iorphans;
  for ( set<tParPtr>::const_iterator it = state.finalState().begin();
	it != state.finalState().end(); ++it ) {
    tColinePtr ic = (**it).origICol();
    tColinePtr oc = (**it).origOCol();
    if ( oc ) {
      ocol.insert(make_pair(oc, *it));
      tColinePair jsource = oc->sourceNeighbours();
      if ( jsource.first && jsource.second ) {
	tJunctionPtr j = state.create<Junction>();
	icol.insert(make_pair(oc, j));
	icol.insert(make_pair(jsource.first, j));
	icol.insert(make_pair(jsource.second, j));
      }
    }
    if ( ic ) {icol.insert(make_pair(ic, *it));
      tColinePair jsink = ic->sinkNeighbours();
      if ( jsink.first && jsink.second ) {
	tJunctionPtr j = state.create<Junction>();
	ocol.insert(make_pair(ic, j));
	ocol.insert(make_pair(jsink.first, j));
	ocol.insert(make_pair(jsink.second, j));
      }
    }
  }
  while ( icol.size() ) {
    MapIt it1 = icol.begin();
    MapIt it2 = ocol.find(it1->first);
    if ( it2 != ocol.end() ) {
      ret.push_back(createDipole(it1->second, it2->second, state));
      ocol.erase(it2);
    } else {
      Throw<QCDFinderException>()
	<< "The QCDDipoleFinder '" << name()
	<< "' was not able to construct dipoles for all colour lines "
	<< "in the final state. This is a serious problem. "
	<< "Please contact the author." << Exception::runerror;
    }
    icol.erase(it1);
  }
  if ( ocol.size() )
    Throw<QCDFinderException>()
      << "The QCDDipoleFinder '" << name()
      << "' was not able to construct dipoles for all colour lines "
      << "in the final state. This is a serious problem. "
      << "Please contact the author." << Exception::runerror;

  // Now we need to connect the dipoles with eachother.
  for ( int i = 0, N = ret.size(); i < N; ++i )
    if ( ret[i]->oPart()->isG() ) {
      for ( int j = 0; j < N; ++j )
	if ( ret[i]->oPart() == ret[j]->iPart() ) {
	  ret[i]->next(ret[j]);
	  ret[j]->prev(ret[i]);
	  break;
	}
      if ( !ret[i]->next() )
      Throw<QCDFinderException>()
	<< "The QCDDipoleFinder '" << name()
	<< "' was not able to construct dipoles for all colour lines "
	<< "in the final state. This is a serious problem. "
	<< "Please contact the author." << Exception::runerror;
    }

  // Finally we need to generate colour indices.
  for ( int i = 0, N = ret.size(); i < N; ++i )
    ret[i]->generateColourIndex();

  return ret;
}

tQCDPtr QCDDipoleFinder::
createDipole(tParPtr pi, tParPtr po, DipoleState & state) const {
  tQCDPtr dip = state.create<QCDDipole>();
  dip->oPart(pi);
  dip->iPart(po);
  // If either parton is a junction, we cannot connect dipoles in the
  // normal way. Instead we setup the junction.
  if ( tJunctionPtr j = dynamic_ptr_cast<tJunctionPtr>(pi) )
    j->replace(tQCDPtr(), dip);
  if ( tJunctionPtr j = dynamic_ptr_cast<tJunctionPtr>(po) )
    j->replace(tQCDPtr(), dip);

  // If both partons come from a resonance we must put them in another
  // colour system.
  tResParPtr rpi = dynamic_ptr_cast<tResParPtr>(pi);
  tResParPtr rpo = dynamic_ptr_cast<tResParPtr>(po);
  if ( rpi && rpo )
    dip->colourIndex(ColourIndex(min(rpi->resonance()->decaySystem(),
				     rpo->resonance()->decaySystem())));
  return dip;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void QCDDipoleFinder::persistentOutput(PersistentOStream & os) const {}

void QCDDipoleFinder::persistentInput(PersistentIStream & is, int) {}


// The following static variable is needed for the type description
// system in ThePEG.
DescribeClass<QCDDipoleFinder,HandlerBase>
  describeAriadne5QCDDipoleFinder("Ariadne5::QCDDipoleFinder",
				  "libAriadne5.so");

void QCDDipoleFinder::Init() {

  static ClassDocumentation<QCDDipoleFinder> documentation
    ("The QCDDipoleFinder class and its sub-classes are responsible for "
     "identifying and introducing of QCD dipoles in the setup phase of a "
     "given DipoleState.");

}

