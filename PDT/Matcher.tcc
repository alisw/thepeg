// -*- C++ -*-
//
// Matcher.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the Matcher class.
//

namespace ThePEG {

template <class T>
Matcher<T>::~Matcher() {
  assert ( initMatcher.check() );
} 

template <class T>
NoPIOClassDescription< Matcher<T> > Matcher<T>::initMatcher;

template <class T>
PMPtr Matcher<T>::Create(const string & newName, string antiName) {
  typedef typename Ptr< Matcher<T> >::pointer MatcherPtr;
  typedef typename Ptr< Matcher<typename T::CC> >::pointer AMatcherPtr;
  PMPtr pm = new_ptr<MatcherPtr>();
  registerRepository(pm, newName);
  if ( typeid(T) == typeid(typename T::CC) ) return pm;
  if ( antiName.empty() ) antiName = newName + "~";
  PMPtr apm = new_ptr<AMatcherPtr>();
  setCC(pm, apm);
  registerRepository(apm, antiName);
  return pm;
}

template <class T>
PMPtr Matcher<T>::pmclone() const {
  return new_ptr(*this);
}

template <class T>
IBPtr Matcher<T>::clone() const {
  return pmclone();
}

template <class T>
IBPtr Matcher<T>::fullclone() const {
  PMPtr pm = pmclone();
  registerRepository(pm);
  if ( !CC() ) return pm;
  PMPtr apm = CC()->pmclone();
  setCC(pm, apm);
  registerRepository(apm);
  return pm;
}

}

