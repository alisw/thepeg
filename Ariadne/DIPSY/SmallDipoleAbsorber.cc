// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SmallDipoleAbsorber class.
//

#include "SmallDipoleAbsorber.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Utilities/Current.h"
#include "DipoleEventHandler.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

SmallDipoleAbsorber::SmallDipoleAbsorber() {}

SmallDipoleAbsorber::~SmallDipoleAbsorber() {}

IBPtr SmallDipoleAbsorber::clone() const {
  return new_ptr(*this);
}

IBPtr SmallDipoleAbsorber::fullclone() const {
  return new_ptr(*this);
}

void printInfo( DipolePtr d ) {
  cout << "dipole " << d << " has size " << d->size()*GeV <<
    ". Has neighbours " << d->neighbors().first << " (" << d->neighbors().first->DGLAPsafe()
       << "), " << d->neighbors().second <<" ("<< d->neighbors().second->DGLAPsafe() << 
"). Has children " << d->children().first << ", " << d->children().second <<
    ". Is DGLAPsafe: " << d->DGLAPsafe() << endl;
}

void SmallDipoleAbsorber::
absorbParton(tPartonPtr p, DipoleState & state ) const {

  DipolePtr d1 = p->dipoles().first;
  DipolePtr d2 = p->dipoles().second;

  if ( !d1 ) {
    if ( !(d2->neighbors().second) )
      cout << "absorbing qqbar pair!!!!!!!!!!!!!!!!!! :o (bad!)" << endl;
    PartonPtr g = d2->partons().second;
    d2->neighbors().second->firstNeighbor(d1);
    g->dipoles(make_pair(d1, d2->neighbors().second));
    d2->children().second = d2->neighbors().second;
    g->plus(g->plus() + p->plus());
    g->updateMomentum();
    g->flavour(p->flavour());
    if ( p->valence() ) {
      g->valencePT(p->valencePT());
      p->parents(make_pair(g,g));
    }
    state.save();
    return;
  }
  if ( !d2 ) {
    if ( !(d1->neighbors().first) )
      cout << "absorbing qqbar pair!!!!!!!!!!!!!!!!!! :o (bad!)" << endl;
    PartonPtr g = d1->partons().first;
    d1->neighbors().first->secondNeighbor(d2);
    g->dipoles(make_pair(d1->neighbors().first, d2));
    d1->children().second = d1->neighbors().first;
    g->plus(g->plus() + p->plus());
    g->updateMomentum();
    g->flavour(p->flavour());
    if ( p->valence() ) {
      g->valencePT(p->valencePT());
      p->parents(make_pair(g,g));
    }
    state.save();
    return;
  }

  PartonPtr p1 = d1->partons().first;
  PartonPtr p2 = d2->partons().second;
  if ( p1 == p2 ) {
    state.sortDipoles();
    if( !(state.forceSwing(p, 5.0, 0.0)) ) {
      cout << "couldnt find swing" << endl;
      return;   //TODO: move to this class. think through.
    }
    d1 = p->dipoles().first;
    d2 = p->dipoles().second;
    p1 = d1->partons().first;
    p2 = d2->partons().second;
  }

  //always save interacting dip, and otherwise the bigger DGLAP safe dip.
  if ( ( (d1->DGLAPsafe() && ( !(d2->DGLAPsafe()) || d1->size() >= d2->size()) ) ||
         ( !(d2->DGLAPsafe()) && d1->size() >= d2->size()) ||
	 ( d1->interacted() )) && !(d2->interacted()) ) {
    d2->children().second = d1;
    p2->dipoles(make_pair(d1,p2->dipoles().second));
    d1->partons(make_pair(p1,p2));
    d1->secondNeighbor(d2->neighbors().second);
    if ( d1->neighbors().second ) d1->neighbors().second->firstNeighbor(d1);
  }
  else {
    d1->children().second = d2;
    p1->dipoles(make_pair(p1->dipoles().first, d2));
    d2->partons(make_pair(p1,p2));
    d2->firstNeighbor(d1->neighbors().first);
    if ( d2->neighbors().first ) d2->neighbors().first->secondNeighbor(d2);
  }

  //Conserve pT and p+ based on dipole sizes. y and p- adapt.
  Parton::Point r1 = p->position() - p1->position();
  Parton::Point r2 = p->position() - p2->position();
  double P1 = sqr(r2.pt())/( sqr(r1.pt()) + sqr(r2.pt()) );

  //dont give p+ to a parton ahead in rapidity if possible.
  if ( p1->y() < p->y() && p->y() < p2->y() ) {
    P1 = 1.0;
  }
  else if ( p2->y() < p->y() && p->y() < p1->y() ) {
    P1 = 0.0;
  }
  else if ( p->y() < p2->y() && p->y() < p1->y() && p1->y() > p2->y() ) {
    P1 = 0.0;
  }
  else if ( p->y() < p2->y() && p->y() < p1->y() && p2->y() > p1->y() ) {
    P1 = 1.0;
  }

  //with valence-status comes always also the p+ though. or?? >_>
  if ( p->valence() ) {
    if ( r1.pt() < r2.pt() || p2->valence() ) {
      p1->valencePT(p->valencePT());
      p->parents(make_pair(p1,p1));
      P1 = 1.0;
    }
    else {
      p2->valencePT(p->valencePT());
      p->parents(make_pair(p2,p2));
      P1 = 0.0;
    }
  }

  double P2 = 1.0 - P1;

  p1->plus( p1->plus() + P1*p->plus() );
  p2->plus( p2->plus() + P2*p->plus() );
  p1->updateMomentum();
  p2->updateMomentum();

  if ( d1->size()/(d1->size()+d2->size()) < UseRandom::rnd() )
    d1->colour(d2->colour());

  state.save();
}

void SmallDipoleAbsorber::setDGLAPsafe ( tDipolePtr dip ) const {
  dip->DGLAPsafe( true );
  InvEnergy firstLength = dip->size();
  InvEnergy secondLength = dip->size();
  if ( dip->interacted() ) {
    firstLength = dip->interactionLengths().first;
    secondLength = dip->interactionLengths().second;
  }
  if ( dip->neighbors().first && !(dip->neighbors().first->interacted()) &&
       dip->neighbors().first->size() > firstLength )
    setDGLAPsafe( dip->neighbors().first );
  if ( dip->neighbors().second && !(dip->neighbors().second->interacted()) &&
       dip->neighbors().second->size() > secondLength )
    setDGLAPsafe( dip->neighbors().second );
}

template < typename Sort >
void SmallDipoleAbsorber::setOrdered ( tPartonPtr p, Sort sorted ) const {
  p->ordered( true );
  if ( p->dipoles().first && sorted( p, p->dipoles().first->partons().first ) )
    setOrdered( p->dipoles().first->partons().first, sorted );
  if ( p->dipoles().second && sorted(p, p->dipoles().second->partons().second) )
    setOrdered( p->dipoles().second->partons().second, sorted );
}

void SmallDipoleAbsorber::DGLAPabsorb(tPartonPtr p, DipoleState & state ) const {
  DipolePtr shortDip, longDip, shortNeighbor;
  PartonPtr shortP, secondShortP;
  InvEnergy firstLength = p->dipoles().first->size();
  InvEnergy secondLength = p->dipoles().second->size();
  if ( p->dipoles().first->interacted() )
    firstLength = p->dipoles().first->interactionLengths().second;
  if ( p->dipoles().second->interacted() )
    secondLength = p->dipoles().second->interactionLengths().first;
  if (firstLength > secondLength ) {
    shortDip = p->dipoles().second;
    longDip = p->dipoles().first;
    shortP = shortDip->partons().second;
    shortNeighbor = shortDip->neighbors().second;
    secondShortP = shortNeighbor->partons().second;
  }
  else {
    shortDip = p->dipoles().first;
    longDip = p->dipoles().second;
    shortP = shortDip->partons().first;
    shortNeighbor = shortDip->neighbors().first;
    secondShortP = shortNeighbor->partons().first;
  }
  if ( p->valence() && shortP->valence() ) {
    setDGLAPsafe(shortDip);
    return;
  }

  if ( shortNeighbor->interacted() ) {
    if ( longDip->interacted() ) {
      //move the larger of the interacting dips
      if ( shortDip->neighbors().first->interactionLengths().second >
	   shortDip->neighbors().second->interactionLengths().first )
	absorbParton( shortDip->partons().first, state );
      else
	absorbParton( shortDip->partons().second, state );
    }
    else
      absorbParton(p, state);
    return; //dont recur
  }

  if ( shortNeighbor->DGLAPsafe() ) {
    if ( longDip->size() > shortNeighbor->size() )
      absorbParton( p, state );
    else
      absorbParton( shortP, state );
    return; //dont recur
  }

  absorbParton(shortP, state);
  if ( p->dist2(*secondShortP) < p->dist2(*shortP) ) {
    DGLAPabsorb(p, state);
  }
}

bool rapSort(PartonPtr part1,PartonPtr part2) {
  return (part1->y() >= part2->y());
}

bool minusSort(PartonPtr p1, PartonPtr p2) {
  return ( p1->minus() >= p2->minus() );
}

bool plusSort(PartonPtr p1, PartonPtr p2) {
  return ( p1->plus() <= p2->plus() );
}

bool SmallDipoleAbsorber::swing( tDipolePtr d1, tDipolePtr d2 ) const {
  if ( DipoleAbsorber::swing(d1, d2) ) return true;
  DipoleState & state = d1->dipoleState();
  state.sortDipoles();
  if( !(state.forceSwing(d1->partons().second, 5.0, 0.0)) ) return false;

  return DipoleAbsorber::swing(d1->partons().first->dipoles().second,
			       d2->partons().second->dipoles().first);

}

// void SmallDipoleAbsorber::isolateNonParticipating( DipoleState & state ) const {
//   for( int i = 0; i < int(state.initialDipoles().size()); i++ ) {
//     DipolePtr d = state.initialDipoles()[i];
//     if ( !(d->participating()) ) {
//       if ( d->partons().first->dipoles().second != d->partons().second->dipoles().first )
// 	swing( d->partons().first->dipoles().second, d->partons().second->dipoles().first );
//       d->partons().first->dipoles().second->participating( false );
//     }
//   }
// }

void SmallDipoleAbsorber::hiddenRecoil(DipoleState & state) const {
  cout << "entering hiddenrecoil--------------------" << endl;
  state.diagnosis(true);
  list<PartonPtr> partons = state.getPartons();
  list<PartonPtr> toBeAbsorbed;
  cout << partons.size() << " number o partons." << endl;
  //reset momenta
  for ( list<PartonPtr>::iterator it = partons.begin();
	it != partons.end(); it++ )
    (*it)->pT(0.0*(*it)->pT());
  for ( list<PartonPtr>::iterator it = partons.begin();
	it != partons.end(); it++ ) {
    tPartonPtr p = *it;
    if ( p->valence() ) {
      p->pT(p->pT() + p->valencePT());
      if ( p->pT().pt() == ZERO )
	cout << "zero momentum valence parton in hiddenrecoil" << endl;
      if ( p->rightMoving() ) {
	p->y( log(p->pT().pt()/p->plus()) );
	p->minus( p->pT().pt()*exp(p->y()) );
      }
      else {
	p->y( log(p->minus()/p->pT().pt()) );
	p->plus( p->pT().pt()*exp(-p->y()) );
      }
      cout << "added valence momentum" << endl;
      state.diagnosis(true);
      continue;
    }
    tPartonPtr p1 = p;
    tPartonPtr p2 = p;
    //find first real ancestor on each side, and recoil with them
    do { p1 = p1->parents().first; } while (p1->absorbed());
    do { p2 = p2->parents().second; } while (p2->absorbed());
    recoil( p, p1, state.ymax() );
    recoil( p, p2, state.ymax() );
    cout << "did normal recoil" << endl;
    state.diagnosis(true);
    if ( p->y() > state.ymax() && !(p->dipoles().first->interacted() &&
			   p->dipoles().second->interacted()) )
      toBeAbsorbed.push_back(p);
  }

  //if any partons needed to be absorbed, absorb them and restart.
  if ( toBeAbsorbed.size() != 0 ) {
    bool absorbed = false;
    for ( list<PartonPtr>::iterator it = toBeAbsorbed.begin();
	  it != toBeAbsorbed.end(); it++ ) {
      cout << "absorbing parton" << endl;
      state.diagnosis(true);
      absorbParton(*it, state);
      if (   (*it)->dipoles().first->children().second  ||
	     (*it)->dipoles().second->children().second  )
	absorbed = true;
    }
    if ( absorbed ) {
      cout << "recur" << endl;
      hiddenRecoil(state);
    }
  }
}

void SmallDipoleAbsorber::
recoil(tPartonPtr p1, tPartonPtr p2, double ymax) const {
  Parton::Point r = p1->position() - p2->position();
  TransverseMomentum recoil = p1->pTScale()*r/sqr(r.pt());
  p1->pT( p1->pT() + recoil );
  p2->pT( p2->pT() - recoil );

  if ( p1->rightMoving() ) {
    p1->y( log(p1->pT().pt()/p1->plus()) );
    p1->minus( p1->pT().pt()*exp(p1->y()) );
  }
  else {
    p1->y( log(p1->minus()/p1->pT().pt()) );
    p1->plus( p1->pT().pt()*exp(-p1->y()) );
  }
  if ( p2->rightMoving() ) {
    p2->y( log(p2->pT().pt()/p2->plus()) );
    p2->minus( p2->pT().pt()*exp(p2->y()) );
  }
  else {
    p2->y( log(p2->minus()/p2->pT().pt()) );
    p2->plus( p2->pT().pt()*exp(-p2->y()) );
  }
  //check if one of the partons is between two interacting dipoles,
  //and got pushed over ymax. Should happen very rarely.
  if ( ( p1->y() > ymax
	 && (p1->dipoles().first->interacted() && p1->dipoles().second->interacted()) )
       ||
       ( p2->y() > ymax
	 && (p2->dipoles().first->interacted() && p2->dipoles().second->interacted()) ) )
    cout << "parton got pushed over ymax!!" << endl;
}

template <typename Sort>
void SmallDipoleAbsorber::absorbNonOrderedPartons( DipoleState & state, Sort ordered ) const {
  list < PartonPtr > partons = state.getPartons();
  list < PartonPtr > marked;

  while ( true ) {
    for ( list<PartonPtr>::iterator p = partons.begin();
	  p != partons.end(); p++ ) {
      (*p)->ordered(false);
    }

    //mark the ones in a sorted chain from an interacting parton
    for ( list<PartonPtr>::iterator p = partons.begin();
	  p != partons.end(); p++ ) {
      if ( (**p).ordered() ) continue;
      if ( ((*p)->dipoles().first && (*p)->dipoles().first->interacted()) ||
	   ((*p)->dipoles().second && (*p)->dipoles().second->interacted()) ) {
	setOrdered( (*p), ordered );
      }
    }

    //manually set the marked partons safe.
    for ( list<PartonPtr>::const_iterator p = marked.begin(); p != marked.end(); p++ ) {
      setOrdered( *p, ordered );
    }

    //find the highest non-ordered parton absP with an ordered neighbour
    PartonPtr absP;
    for ( list<PartonPtr>::iterator p = partons.begin();
	  p != partons.end(); p++ ) {
      if ( !((*p)->ordered()) && 
	   ( ((*p)->dipoles().first &&
	      (*p)->dipoles().first->partons().first->ordered()) ||
	     ((*p)->dipoles().second &&
	      (*p)->dipoles().second->partons().second->ordered()) ) )
	if ( !absP || ordered( *p, absP) ) absP = *p;
    }

    //check that the unordered partons isnt a valence parton.
    if ( absP && absP->valence() ) {
      bool firstCanBeAbsorbed = false;
      bool secondCanBeAbsorbed = false;
      if ( absP->dipoles().first ) {
	PartonPtr p1 = absP->dipoles().first->partons().first;
        if ( !(p1->valence()) &&
	   !( (p1->dipoles().first && p1->dipoles().first->interacted()) ||
	      (p1->dipoles().second && p1->dipoles().second->interacted()) ) &&
	     p1->ordered() )
	  firstCanBeAbsorbed = true;
      }
     if ( absP->dipoles().second ) {
	PartonPtr p2 = absP->dipoles().second->partons().second;
        if ( !(p2->valence()) &&
	   !( (p2->dipoles().first && p2->dipoles().first->interacted()) ||
	      (p2->dipoles().second && p2->dipoles().second->interacted()) ) &&
	   p2->ordered() )
	  secondCanBeAbsorbed = true;
      }
     if ( firstCanBeAbsorbed && !secondCanBeAbsorbed ) {
       absP = absP->dipoles().first->partons().first;
     }
     else if ( !firstCanBeAbsorbed && secondCanBeAbsorbed ) {
       absP = absP->dipoles().second->partons().second;
     }
     else if ( firstCanBeAbsorbed && secondCanBeAbsorbed ) {
       if ( ordered ( absP->dipoles().first->partons().first,
                      absP->dipoles().second->partons().second ) )
         absP = absP->dipoles().first->partons().first;
       else
         absP = absP->dipoles().second->partons().second;
     }
     else if ( !firstCanBeAbsorbed && !secondCanBeAbsorbed ) {
       marked.push_back(absP);
       continue;
     }

    }

    if ( absP ) {
      absorbParton( absP, state);
      partons.remove(absP);
    } else  { //if there is none, look for unordered loops
      PartonPtr highestP;
      for ( list<PartonPtr>::iterator p = partons.begin();
	    p != partons.end(); p++ ) {
	if ( !((*p)->ordered()) )
	  if ( !highestP || ordered((*p), highestP) )
	    highestP = *p; 
      }
      if ( highestP ) { //mark the highest parton in the loop as safe
	marked.push_back(highestP);
      }
      else {//no unordered partons, we are done
	break;
      }
    }
  }
}

bool SmallDipoleAbsorber::checkDGLAPsafe ( InvEnergy r1, InvEnergy r2 ) const {
  return sqr(r1)/Current<DipoleEventHandler>()->alphaS(r1) >
    (sqr(r2)/Current<DipoleEventHandler>()->alphaS(r2))*UseRandom::rnd();
}

void SmallDipoleAbsorber::absorbSmallDipoles( DipoleState & state ) const {
  set<DipolePtr> marked;  //the manually marked safe dipoles.
  //start by adding the interacting dipoles as safe.
  list< DipolePtr > dipoles = state.getDipoles();
  InvEnergy minScale = 0.0/GeV;
  for ( list<DipolePtr>::const_iterator it = dipoles.begin(); 
	it != dipoles.end(); it++ ) { DipolePtr dip = *it;
    if ( dip->interacted() ) {
      InvEnergy thisScale =
	min( sqrt(dip->partons().first->dist2(*(dip->interacted()->partons().second))),
	     sqrt(dip->partons().second->dist2(*(dip->interacted()->partons().first))));
      if ( minScale == 0.0/GeV || thisScale < minScale )
	minScale = thisScale;
      marked.insert(dip);
    }
  }
  //find all loops, and decide if they should be absorbed or not.
  // first mark all partons as nonordered, so it can be used for safe loops.
  list<PartonPtr> partons = state.getPartons();
for ( list<PartonPtr>::const_iterator p = partons.begin(); p != partons.end(); p++ ) {
  (*p)->ordered(false);
 }
 set< list<PartonPtr> > loops = state.loops();
 for ( set< list<PartonPtr> >::iterator it = loops.begin();
       it != loops.end(); it++ ) { list<PartonPtr> loop = *it;
   bool keep = false;
   DipolePtr maxDip;

   for( list<PartonPtr>::const_iterator jt = loop.begin();
	jt != loop.end(); jt++ ) { PartonPtr p = *jt;
     if ( !(p->dipoles().first) ) continue;
     if ( p->dipoles().first->interacted() ) keep = true;
     if ( p->ordered() ) keep = true;
     if ( !maxDip || p->dipoles().first->size() > maxDip->size() )
       maxDip = p->dipoles().first;
   }

   if ( !keep && checkDGLAPsafe( maxDip->size(), minScale ) ) {
     keep = true;
     maxDip->partons().first->ordered(true);
   }
   if ( !keep ) { //otherwise absorb loop, update and restart loop iteration
     absorbLoop( (*(loop.begin()))->dipoles().first, state );
     loops = state.loops();
     it = loops.begin();
   }
 }

 //now mark the largest dips in non-interacting loop
 for ( set< list<PartonPtr> >::iterator it = loops.begin();
       it != loops.end(); it++ ) { list<PartonPtr> loop = *it;
   bool interacting = false;
   DipolePtr maxDip;
   for( list<PartonPtr>::const_iterator jt = loop.begin();
	jt != loop.end(); jt++ ) { PartonPtr p = *jt;
     if ( !(p->dipoles().first) ) continue;
     if ( !maxDip || p->dipoles().first->size() > maxDip->size() )
       maxDip = p->dipoles().first;
     if ( p->dipoles().first->interacted() )
       interacting = true;
   }
   if ( !interacting )
     marked.insert(maxDip);
 }
  //now look through all non-safe dipoles, and mark them as safe, or absorb them.
  while( true ) {
    //first update the lists and safe flags
    for( list<DipolePtr>::iterator it = dipoles.begin();
	 it != dipoles.end(); it++ ) { DipolePtr dip = *it;
      if ( dip->children().first || dip->children().second ) {
	it = dipoles.erase(it);
	it--;
	continue;
      }
      else if ( dip->partons().first->valence() && dip->partons().second->valence() ) {
	marked.insert(dip);
      }
      dip->DGLAPsafe(false);
    }
    for( set<DipolePtr>::iterator it = marked.begin();
	 it != marked.end(); it++ ) { DipolePtr dip = *it;
      setDGLAPsafe(dip);
    }
    //then find a safe - nonsafe neighbouring pair.
    DipolePtr safeDip;
    DipolePtr nonSafeDip;
    bool found = false;
    InvEnergy smallest = 0.0/GeV;
    for( list<DipolePtr>::iterator it = dipoles.begin();
	 it != dipoles.end(); it++ ) { DipolePtr dip = *it;
      if ( !(dip->DGLAPsafe()) ) {
	if ( dip->neighbors().first->DGLAPsafe() &&
	     (!found || dip->size() < smallest) ) {
	  nonSafeDip = dip;
	  safeDip = dip->neighbors().first;
	  found = true;
	  smallest = dip->size();
	}
	if ( dip->neighbors().second->DGLAPsafe() &&
             (!found || dip->size() < smallest || 
	      (dip == nonSafeDip &&
	       dip->neighbors().second->size() < safeDip->size()) ) ) {
	  nonSafeDip = dip;
	  safeDip = dip->neighbors().second;
	  found = true;
	  smallest = dip->size();
	}
      }
    }
    //if a pair is found, decide if the nonsafe one should be absorbed or marked
    if ( found ) {
      InvEnergy safeScale = safeDip->size();
      if ( safeDip->interacted() ) { //use the interacted dipsize for interacting dips
	if ( nonSafeDip == safeDip->neighbors().first )
	  safeScale = safeDip->interactionLengths().first;
	else
	  safeScale = safeDip->interactionLengths().second;
      }
      if ( checkDGLAPsafe( nonSafeDip->size(), safeScale ) ) {
	marked.insert(nonSafeDip);
      }
      else {
	if ( safeDip == nonSafeDip->neighbors().first )
	  DGLAPabsorb(safeDip->partons().second, state);
	else if ( safeDip == nonSafeDip->neighbors().second )
	  DGLAPabsorb(safeDip->partons().first, state);
	else cout << "safedip and nonsafe dip not neighbors! :(" << endl;
      }
    }
    //if no nonsafe dips left, we're done.
    else
      break;
  }
}

void SmallDipoleAbsorber::reabsorb( DipoleState & state) const {
  realValence( state );
  state.save();
  absorbVirtualPartons( state );
  // isolateNonParticipating( state );
  absorbSmallDipoles( state );
  absorbNonOrderedPartons( state, minusSort );
  absorbNonOrderedPartons( state, plusSort );
  hiddenRecoil( state );
}

void SmallDipoleAbsorber::realValence( DipoleState & state ) const {
  vector<DipolePtr> dips = state.initialDipoles();
  for ( int i = 0; i < int(dips.size());i++) {
    dips[i]->partons().first->interact();
    dips[i]->partons().second->interact();
  }
}

DipolePtr SmallDipoleAbsorber::largestNonsafeDipole(list<DipolePtr> dipoles,
						    DipoleState & state ) const {
  DipolePtr largestDip;
  InvEnergy DGLAPScale = 0.0/GeV;
  InvEnergy size = 0.0/GeV;
  //find largest dipole that is not safe, and the smallest safe dipole.
  // ***************** ATTENTION PHYSICS HERE!! WHAT SCALE TO PICK!?!?! ***********
  for ( list<DipolePtr>::const_iterator dPtr = dipoles.begin(); 
	dPtr != dipoles.end(); dPtr++ ) {
    if ( !((*dPtr)->children().first || (*dPtr)->children().second )) {
      if ( (*dPtr)->DGLAPsafe() ) {
	if (DGLAPScale == 0.0/GeV || (*dPtr)->size() < DGLAPScale )
	  DGLAPScale = (*dPtr)->size();
      }
      else if( (*dPtr)->size() > size ) {
	largestDip = *dPtr;
	size = (*dPtr)->size();
      }
    }
  }
  //compare them, and mark as safe or absorb the entire loop.
  //do not automatically save valence chains, as this creates many high pT events
  if ( largestDip ) {
    if ( (sqr(largestDip->size())/Current<DipoleEventHandler>()->alphaS(largestDip->size())) >
	 (sqr(DGLAPScale)/Current<DipoleEventHandler>()->alphaS(DGLAPScale))*
	 UseRandom::rnd() ) {
      return largestDip;
    }
    else {
      absorbLoop( largestDip, state );
      return largestNonsafeDipole(dipoles, state);
    }
  }
  else {
    return largestDip;
  }
}

void SmallDipoleAbsorber::absorbLoop( DipolePtr d, DipoleState & state ) const {
  //make a list of the partons and dipoles in the loop.
  list<PartonPtr> partons;
  list<DipolePtr> dipoles;
  tPartonPtr p = d->partons().second;
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
  if ( dipoles.size() == 0 )
    cout << "wtf, no dipoles in loop!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1" << endl;
  if ( otherDips.size() == 0 ) {
    cout << "wtf, no other dipoles in loop!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1" << endl;
    state.diagnosis(true);
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
  //then sort the partons in rapidity, and absorb them.
  //leave all valence partons, or if none, the lowest rapidity parton.
  partons.sort(rapSort);
  bool valence;
  for ( list<PartonPtr>::const_iterator it = partons.begin() ;
	(!valence && it != partons.end()--) || (valence && it == partons.end()); it++ ) {
    if ( (*it)->valence() )
      valence = true;
    else
      absorbParton((*it), state);
  }
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void SmallDipoleAbsorber::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void SmallDipoleAbsorber::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<SmallDipoleAbsorber,DIPSY::DipoleAbsorber>
  describeDIPSYSmallDipoleAbsorber("DIPSY::SmallDipoleAbsorber", "libAriadne5.so libDIPSY.so");


void SmallDipoleAbsorber::Init() {

  static ClassDocumentation<SmallDipoleAbsorber> documentation
    ("The SmallDipoleAbsorber class implements a model describing how "
     "non-interacted dipoles will be reabsorbed. The model is based on "
     "the removal of small dipoles.");

}

