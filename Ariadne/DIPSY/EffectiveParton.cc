// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EffectiveParton class.
//

#include "EffectiveParton.h"
#include "Dipole.h"
#include "DipoleState.h"
#include "DipoleEventHandler.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/MaxCmp.h"
#include "ThePEG/Utilities/UnitIO.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

EffectiveParton::EffectiveParton() {}

template <typename T>
bool sortof(const T & t1, const T & t2) {
  return sqr(t1 - t2 ) <= 1.0e-16*sqr(t1 + t2);
}

template <typename T>
bool sortof(const pair<T,T> & t1, const pair<T,T> & t2) {
  return sortof(t1.first, t2.first) && sortof(t1.second, t2.second);
}

template <typename T>
bool sortof(const Transverse<T> & t1, const Transverse<T> & t2) {
  return sortof(t1.first, t2.first) && sortof(t1.second, t2.second);
}

EffectiveParton::EffectiveParton(Parton & p): Parton(p) {}

void EffectiveParton::uncache() {
    plus(cache[currange].plus);
    pT(cache[currange].pT);
    y(cache[currange].y);
    minus(cache[currange].minus);
    dipoles(cache[currange].dipair);
}

void EffectiveParton::recache() {
  for ( int ir = 0, Nr = cache.size(); ir < Nr; ++ir ) {
    cache[ir].plus = ZERO;
    cache[ir].pT = TransverseMomentum();
    for ( int i = 0, N = cache[ir].partons.size(); i < N; ++i ) {
      cache[ir].plus += cache[ir].partons[i]->plus();
      cache[ir].pT += cache[ir].partons[i]->pT();
    }
    cache[ir].y = log(cache[ir].pT.pt()/cache[ir].plus);
    cache[ir].minus = cache[ir].pT.pt()*exp(cache[ir].y);
  }
  uncache();
}

/*
void EffectiveParton::checkCache() const {
  if ( cache.empty() ) return;
  if ( internalPartons().empty() ) return;
  if ( !sortof(plus(), cache[currange].plus) )
    breakThePEG();
  if ( !sortof(pT(), cache[currange].pT) )
    breakThePEG();
  if (internalPartons().size() != cache[currange].partons.size() )
    breakThePEG();
  set<tPartonPtr> ps = internalPartons();
  for ( int i = 0, N = cache[currange].partons.size(); i < N; ++i )
    if ( !ps.erase(cache[currange].partons[i]) ) breakThePEG();
  if (theInternalDipoles.size() != cache[currange].dipoles.size() )
    breakThePEG();
  else
    for ( int i = 0, N = cache[currange].dipoles.size(); i < N; ++i )
      if ( cache[currange].dipoles[i] != theInternalDipoles[i] )
	breakThePEG();
  if ( dipoles().first != cache[currange].dipair.first ||
       dipoles().second != cache[currange].dipair.second )
    breakThePEG();
}
*/

EffectivePartonPtr EffectiveParton::create(Parton & p, InvEnergy range) {
  EffectivePartonPtr effp = new_ptr(EffectiveParton(p));
  effp->theOriginalParton = & p;
  if ( Current<DipoleEventHandler>()->effectivePartonMode() >= 2 )
    effp->setupCache(range);
  else
    effp->setRange(range);
  while ( effp->pT().pt() <= 10.0*Constants::epsilon*p.pT().pt() ) {
    MaxCmp<InvEnergy2> maxdist;
    // If pt is zero, shrink the range below the parton which is furthest away.
    for ( set<tPartonPtr>::const_iterator it = effp->internalPartons().begin();
	  it != effp->internalPartons().end(); ++it )
      maxdist(effp->dist2(**it));
    effp->setRange((1.0 - 10.0*Constants::epsilon)*sqrt(maxdist.value()));
  }
  effp->y(log(effp->pT().pt()/effp->plus()));
  effp->minus(effp->pT().pt2()/effp->plus());

  return effp;
}

EffectiveParton::~EffectiveParton() {}

Ariadne5::ClonePtr EffectiveParton::clone() const {
  return new_ptr(*this);
}

void EffectiveParton::rebind(const TranslationMap & trans) {
  theOriginalParton = trans.translate(theOriginalParton);
  vector<tDipolePtr> idips;
  idips.swap(theInternalDipoles);
  trans.translate(inserter(theInternalDipoles), idips.begin(), idips.end());
  set<tPartonPtr> ipars;
  ipars.swap(theInternalPartons);
  trans.translate(inserter(theInternalPartons), ipars.begin(), ipars.end());
  for ( int ir = 0, Nr = cache.size(); ir < Nr; ++ir ) {
    vector<tPartonPtr> pv;
    pv.swap(cache[ir].partons);
    trans.translate(inserter(cache[ir].partons), pv.begin(), pv.end());
    vector<tDipolePtr> dv;
    dv.swap(cache[ir].dipoles);
    trans.translate(inserter(cache[ir].dipoles), dv.begin(), dv.end());
    cache[ir].dipair.first = trans.translate(cache[ir].dipair.first);
    cache[ir].dipair.second = trans.translate(cache[ir].dipair.second);
  }
}

void EffectiveParton::printData() {
  const set<tPartonPtr> & partons = internalPartons();
  Energy sum = 0.0*GeV;
  for ( set<tPartonPtr>::const_iterator it = partons.begin();
	it != partons.end(); it++ ) {
    tPartonPtr p = *it;
    sum += p->plus();
  }
  if ( abs(sum/plus() - 1.0) > 0.000001 || true ) {
    for ( set<tPartonPtr>::iterator it = partons.begin(); it != partons.end(); it++ ) {
      tPartonPtr p = *it;
    }
  }
}

void EffectiveParton::setRange( InvEnergy range ) {
  if ( cache.size() ) {
    currange = 1;
    while ( currange < int(cache.size()) &&
	    sqr(range) >= cache[currange].range ) ++currange;
    --currange;
    uncache();
    return;
  } 
  theInternalDipoles.clear();
  theInternalPartons.clear();
  dipoles( theOriginalParton->dipoles() );
  addPartons(range);
  if ( cache.size() ) uncache();
}

void EffectiveParton::setupCache(InvEnergy maxrange) {
  if ( maxrange > Current<DipoleEventHandler>()->coherenceRange() )
    maxrange = Current<DipoleEventHandler>()->coherenceRange();
  InvEnergy2 cut = sqr(maxrange);
  cache.clear();
  cache.resize(1);
  cache.back().range = ZERO;
  cache.back().plus = theOriginalParton->plus();
  cache.back().pT = theOriginalParton->pT();
  cache.back().partons.push_back(theOriginalParton);
  cache.back().dipair = dipoles();

  pair<InvEnergy2,InvEnergy2> ranges(cut*2.0, cut*2.0);
  pair<tDipolePtr, tDipolePtr> dips = dipoles();
  if ( dips.first && dips.first->partons().first != theOriginalParton )
    ranges.first = dips.first->partons().first->dist2(*this);
  if ( dips.second && dips.second->partons().second != theOriginalParton )
    ranges.second = dips.second->partons().second->dist2(*this);

  while ( min(ranges.first, ranges.second) <= cut ) {
    tPartonPtr mergeparton;
    tDipolePtr mergedipole;
    InvEnergy2 range = ZERO;
    if ( ranges.first < ranges.second ) {
      mergeparton = dips.first->partons().first;
      range = ranges.first;
      mergedipole = dips.first;
      dips.first = dips.first->neighbors().first;
      if ( dips.first && dips.first->partons().first != theOriginalParton )
	ranges.first = dips.first->partons().first->dist2(*this);
    } else {
      mergeparton = dips.second->partons().second;
      range = ranges.second;
      mergedipole = dips.second;
      dips.second = dips.second->neighbors().second;
      if ( dips.second && dips.second->partons().second != theOriginalParton )
	ranges.second = dips.second->partons().second->dist2(*this);
    }
    if ( range > cache.back().range ) {
      cache.resize(cache.size() + 1, cache.back());
      cache.back().range = range;
    }
    cache.back().partons.push_back(mergeparton);
    cache.back().dipoles.push_back(mergedipole);
    cache.back().dipair = dips;
  }

  currange = cache.size() - 1;
  recache();

}

void EffectiveParton::addPartons( InvEnergy range ) {
  if ( range > Current<DipoleEventHandler>()->coherenceRange() )
    range = Current<DipoleEventHandler>()->coherenceRange();
  if ( Current<DipoleEventHandler>()->effectivePartonMode() == 1 )
    return addRelatives(range, theOriginalParton);
  PartonPtr mergeParton;
  if ( dipoles().first &&
       dipoles().first->partons().first->dist2( *this ) <= sqr(range) &&
       dipoles().first->partons().first != theOriginalParton ) {
    mergeParton = dipoles().first->partons().first;
    theInternalDipoles.push_back( dipoles().first );
    dipoles( make_pair(dipoles().first->neighbors().first, dipoles().second));
    addPartons( range );
  }
  else if ( dipoles().second &&
            dipoles().second->partons().second->dist2( *this ) <= sqr(range) &&
	    dipoles().second->partons().second != theOriginalParton ) {
    mergeParton = dipoles().second->partons().second;
    theInternalDipoles.push_back( dipoles().second );
    dipoles( make_pair(dipoles().first, dipoles().second->neighbors().second));
    addPartons( range );
  }
  if ( mergeParton ) merge(mergeParton);
  else {
    plus(ZERO);
    pT(TransverseMomentum());
    merge(theOriginalParton);
  }
}

void EffectiveParton::addRelatives( InvEnergy range, tPartonPtr p ) {
  if ( !p ) return;
  //if already added, dont do anything
  if ( theInternalPartons.find(p) != theInternalPartons.end() ) return;
  //we may only have one valence in an effective parton
  //if not in range, dont do anything
  if ( p->dist2(*this) > sqr(range) ) return;
  //otherwise add the new parton
  merge(p);
  //and recur to parents and children in range
  addRelatives(range, p->parents().first);
  addRelatives(range, p->parents().second);
  if ( p->children().empty() ) return;
  const set<tPartonPtr> & kids = p->children();
  for ( set<tPartonPtr>::const_iterator it = kids.begin();
	it != kids.end(); it++) {
    addRelatives(range, *it);
  }
}

void EffectiveParton::merge(tPartonPtr p) {
  //  if ( theInternalPartons.find(p) != theInternalPartons.end() ) return;
  if ( theInternalPartons.empty() ) {
    plus(ZERO);
    pT(TransverseMomentum());
  }
  if ( theInternalPartons.insert(p).second ) {
    plus( plus() + p->plus() );
    pT( pT() + p->pT() );
  }
}

bool EffectiveParton::changed() {
  if ( cache.size() ) {
    for ( int i = 0, N = cache[currange].dipoles.size(); i < N; ++i )
      if ( cache[currange].dipoles[i]->children().first ||
	   cache[currange].dipoles[i]->children().second )
	return true;
  } else {
    for (int i = 0; i < int(theInternalDipoles.size()); i++ ) {
      if ( theInternalDipoles[i]->children().first ||
	   theInternalDipoles[i]->children().second ) {
	return true;
      }
    }
  }
  if ( (dipoles().first &&  dipoles().first->children().first) ||
       (dipoles().first &&  dipoles().first->children().second) ||
       (dipoles().second &&  dipoles().second->children().first) ||
       (dipoles().second &&  dipoles().second->children().second) ) {
    return true;
  }

  TransverseMomentum sum = TransverseMomentum();
  if ( cache.size() ) {
    for ( int i = 0, N = cache[currange].partons.size(); i < N; ++i )
      sum += cache[currange].partons[i]->pT();
  } else {
    const set<tPartonPtr> & partons = internalPartons();
    for ( set<tPartonPtr>::const_iterator it = partons.begin();
	  it != partons.end(); it++ ) {
      tPartonPtr p = *it;
      sum += p->pT();
    }
  }
  if ( (sum - pT()).pt() > 0.000000001*GeV ) return true;

  return false;
}

bool EffectiveParton::
recoilsOverYMax(TransverseMomentum recoilPT,
		Energy recoilPlus, double yMax) const {
  if ( cache.size() ) {
    for ( int i = 0, N = cache[currange].partons.size(); i < N; ++i ) {
      tPartonPtr p = cache[currange].partons[i];
      if ( p != theOriginalParton &&
	   log(p->pT().pt()/(p->plus() - recoilPlus*p->plus()/plus())) > yMax )
	return true;
    }
  } else {
    const set<tPartonPtr> & partons = internalPartons();
    for ( set<tPartonPtr>::const_iterator it = partons.begin();
	  it != partons.end(); it++ ) {
      tPartonPtr p = *it;
      if ( p != theOriginalParton &&
	   log( p->pT().pt()/(p->plus() - recoilPlus*p->plus()/plus()) ) > yMax )
	return true;
    }
  }
  if ( log( (theOriginalParton->pT() - recoilPT).pt()/
	    (theOriginalParton->plus()*(1.0 - recoilPlus/plus())) ) > yMax )
    return true;
  return false;
}

void EffectiveParton::recoil( TransverseMomentum recoilPT, Energy recoilPlus ) {
  int opteffw = Current<DipoleEventHandler>()->eventFiller().effectiveWeights();
  if ( cache.size() ) {
    for ( int i = 0, N = cache[currange].partons.size(); i < N; ++i ) {
      tPartonPtr p = cache[currange].partons[i];
      if ( opteffw == 0 )
	p->pT( p->pT() - recoilPT*p->plus()/plus() );
      else if ( opteffw == 1 )
	p->pT( p->pT() - recoilPT/double(N) );
      else if ( opteffw == 2 ) {
	if ( p == theOriginalParton ) p->pT( p->pT() - recoilPT );
      }
      p->plus( p->plus() - recoilPlus*p->plus()/plus() );
      p->y( log(p->pT().pt()/p->plus()) );
      p->minus( p->pT().pt()*exp(p->y()) );
    }
    recache();
  } else {
    const set<tPartonPtr> & partons = internalPartons();
    for ( set<tPartonPtr>::const_iterator it = partons.begin();
	  it != partons.end(); it++ ) {
      tPartonPtr p = *it;
      if ( opteffw == 0 )
	p->pT( p->pT() - recoilPT*p->plus()/plus() );
      else if ( opteffw == 1 )
	p->pT( p->pT() - recoilPT/double(partons.size()) );
      else if ( opteffw == 2 ) {
	if ( p == theOriginalParton ) p->pT( p->pT() - recoilPT );
      }
      p->plus( p->plus() - recoilPlus*p->plus()/plus() );
      p->y( log(p->pT().pt()/p->plus()) );
      p->minus( p->pT().pt()*exp(p->y()) );
    }
    pT( pT() - recoilPT );
    plus( plus() - recoilPlus );
    y( log(pT().pt()/plus()) );
    minus( pT().pt()*exp(y()) );
  }
}

void EffectiveParton::newRecoil( Energy recoilPlus ) {
  //first update pT of the original parton (the only pt that changes, since
  // its the only dipole that changes.)
  theOriginalParton->updateMomentum();
  //then recoil plus, and update y and minus correspondingly
  if ( cache.size() ) {
    for ( int i = 0, N = cache[currange].partons.size(); i < N; ++i ) {
      tPartonPtr p = cache[currange].partons[i];
      p->plus( p->plus() - recoilPlus*p->plus()/plus() );
      p->y( log(p->pT().pt()/p->plus()) );
      if ( p->y() > 0.0 )
	p->minus( p->pT().pt()*exp(p->y()) );
    } 
  } else {
    const set<tPartonPtr> & partons = internalPartons();
    for ( set<tPartonPtr>::const_iterator it = partons.begin();
	  it != partons.end(); it++ ) {
      tPartonPtr p = *it;
      p->plus( p->plus() - recoilPlus*p->plus()/plus() );
      p->y( log(p->pT().pt()/p->plus()) );
      if ( p->y() > 0.0 )
	p->minus( p->pT().pt()*exp(p->y()) );
    }
  }
  //update total momentum
  plus( plus() - recoilPlus );
  updateMomentum();
}

bool EffectiveParton::checkSums() const {
  TransverseMomentum sumPT = TransverseMomentum();
  Energy sumPlus = ZERO;
  if ( cache.size() ) {
    for ( int i = 0, N = cache[currange].partons.size(); i < N; ++i ) {
      sumPT += cache[currange].partons[i]->pT();
      sumPlus += cache[currange].partons[i]->plus();
    }
  } else {
    const set<tPartonPtr> & partons = internalPartons();
    for ( set<tPartonPtr>::const_iterator it = partons.begin();
	  it != partons.end(); it ++ ) {
      sumPT += (*it)->pT();
      sumPlus += (*it)->plus();
    }
  }
  if ( (sumPT - pT()).pt() > 0.00000001*GeV ) {
    dipoles().first->dipoleState().diagnosis(true);
    return false;
  }
  if ( abs(sumPlus - plus()) > 0.00000001*GeV ) {
    return false;
  }
  return true;
}

void EffectiveParton::persistentOutput(PersistentOStream & os) const {
  os << theOriginalParton << theInternalDipoles << theInternalPartons;
  os << cache.size();
  for ( int ir = 0, Nr = cache.size(); ir < Nr; ++ ir )
    os << ounit(cache[ir].range, sqr(InvGeV)) << ounit(cache[ir].pT, GeV)
       << ounit(cache[ir].plus, GeV) << ounit(cache[ir].minus, GeV)
       << cache[ir].partons << cache[ir].y << cache[ir].dipoles
       << cache[ir].dipair;
}

void EffectiveParton::persistentInput(PersistentIStream & is, int) {
  is >> theOriginalParton >> theInternalDipoles >> theInternalPartons;
  long Nr = 0;
  is >> Nr;
  cache.resize(Nr);
  for ( int ir = 0; ir < Nr; ++ ir )
    is >> iunit(cache[ir].range, sqr(InvGeV)) >> iunit(cache[ir].pT, GeV)
       >> iunit(cache[ir].plus, GeV) >> iunit(cache[ir].minus, GeV)
       >> cache[ir].partons >> cache[ir].y >> cache[ir].dipoles
       >> cache[ir].dipair;  
}

// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<EffectiveParton,DIPSY::Parton>
  describeDIPSYEffectiveParton("DIPSY::EffectiveParton", "libAriadne5.so libDIPSY.so");

void EffectiveParton::Init() {}

