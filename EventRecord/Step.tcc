// -*- C++ -*-
//
// Step.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the Step class.
//

namespace ThePEG {

template <typename OutputIterator>
void Step::select(OutputIterator r, const SelectorBase & s) const {
  if ( s.finalState() ) copyIfCheck(r, particles(), s);
  if ( s.intermediate() ) copyIfCheck(r, intermediates(), s);
}

template <typename Iterator>
void Step::addParticles(Iterator first, Iterator last) {
  theParticles.insert(first, last);
  allParticles.insert(first, last);
  if ( collision() ) collision()->addParticles(first, last);
  for ( ; first != last; ++first )
    if ( !(**first).birthStep() )
      (**first).rep().theBirthStep = this;
}

template <typename Iterator>
void Step::addIntermediates(Iterator first, Iterator last) {
  theIntermediates.insert(first, last);
  allParticles.insert(first, last);
  if ( collision() ) collision()->addParticles(first, last);
  for ( ; first != last; ++first ) {
    if ( !(**first).birthStep() ) (**first).rep().theBirthStep = this;
    ParticleSet::iterator pit = theParticles.find(*first);
    if ( pit != theParticles.end() ) theParticles.erase(pit);
  }
}

template <typename Iterator>
bool Step::
addDecayProduct(Iterator firstParent, Iterator lastParent, tPPtr child,
		bool checkfinal) {
  if ( !collision() ) return false;
  if ( collision()->finalStep() != this ) return false;
  if ( checkfinal ) {
    for ( Iterator it = firstParent; it != lastParent; ++it ) {
      tPPtr parent = const_ptr_cast<tPPtr>((**it).final());
      if ( member(theParticles, parent) ) continue;
      if ( parent->children().empty() ||
	   !member(theParticles, parent->children()[0]->final()) ) return false;
    }
  }
  for ( Iterator it = firstParent; it != lastParent; ++it ) {
    tPPtr parent = const_ptr_cast<tPPtr>((**it).final());
    ParticleSet::iterator pit = theParticles.find(parent);
    if ( pit != theParticles.end() ) {
      theParticles.erase(pit);
      if ( parent->birthStep() == this ) theIntermediates.insert(parent);
    }
    parent->rep().theChildren.push_back(child);
    child->rep().theParents.push_back(parent);
  }
  child->rep().theBirthStep = this;
  addParticle(child);
  return true;
}

template <typename PIterator, typename CIterator>
bool Step::addDecayProduct(PIterator firstParent, PIterator lastParent,
			   CIterator firstChild, CIterator lastChild) {
  if ( !collision() ) return false;
  if ( collision()->finalStep() != this ) return false;
  for ( PIterator it = firstParent; it != lastParent; ++it ) {
    tPPtr parent = const_ptr_cast<tPPtr>((**it).final());
    if ( member(theParticles, parent) ) continue;
    if ( parent->children().empty() ||
	 !member(theParticles, parent->children()[0]->final()) ) return false;
  }
  for ( PIterator it = firstParent; it != lastParent; ++it ) {
    tPPtr parent = const_ptr_cast<tPPtr>((**it).final());
    ParticleSet::iterator pit = theParticles.find(parent);
    if ( pit != theParticles.end() ) {
      theParticles.erase(pit);
      if ( parent->birthStep() == this ) theIntermediates.insert(parent);
    }
    for ( CIterator cit = firstChild; cit != lastChild; ++cit ) {
      parent->rep().theChildren.push_back(*cit);
      (**cit).rep().theParents.push_back(parent);
    }
  }
  for ( CIterator cit = firstChild; cit != lastChild; ++cit ) {
    (**cit).rep().theBirthStep = this;
    addParticle(*cit);
  }
  return true;
}

template <typename Iterator>
tParticleSet Step::getCurrent(Iterator first, Iterator last) const {
  tParticleSet res;
  for ( Iterator i = first; i != last; ++i )
    addIfFinal(inserter(res), *i);
  return res;
}

template <typename Inserter, typename PPointer>
void Step::addIfFinal(Inserter o, PPointer p) {
  if ( member(theParticles, p) ) *o++ = p;
  else if ( p->next() ) addIfFinal(o, p->next());
  else
    for ( int i = 0, N = p->children().size(); i < N; ++i )
      addIfFinal(o, p->children()[i]);
}

}
