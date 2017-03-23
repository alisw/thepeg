// -*- C++ -*-
//
// HandlerGroup.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HandlerGroup class.
//

#include "ThePEG/Handlers/HandlerGroup.h"
#include "ThePEG/Handlers/StepHandler.h"
#include "ThePEG/Handlers/Hint.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/Throw.h"
#include <algorithm>

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
#include "HandlerGroup.tcc"
#endif

using namespace ThePEG;

HandlerGroupBase::HandlerGroupBase()
  : isEmpty(true) {}

HandlerGroupBase::HandlerGroupBase(const HandlerGroupBase & hg)
  : isEmpty(hg.isEmpty), theDefaultPreHandlers(hg.theDefaultPreHandlers),
    theDefaultPostHandlers(hg.theDefaultPostHandlers),
    thePreHandlers(hg.thePreHandlers), theHints(hg.theHints),
    thePostHandlers(hg.thePostHandlers) {}

HandlerGroupBase::~HandlerGroupBase() {}

HandlerGroupBase::StepWithHint HandlerGroupBase::next() {
  StepWithHint sh = make_pair(StepHdlPtr(), HintPtr());
  if ( isEmpty ) return sh;
  if (  !thePreHandlers.empty() ) {
    sh = thePreHandlers.back();
    thePreHandlers.pop_back();
    return sh;
  }
  if ( handler() ) {
    if ( !theHints.empty() ) {
      sh.first = handler();
      sh.second = theHints.back();
      theHints.pop_back();
      return sh;
    } else
      setHandler();
  }
  if (  !thePostHandlers.empty() ) {
    sh = thePostHandlers.back();
    thePostHandlers.pop_back();
    return sh;
  }
  isEmpty = true;
  return sh;
}

void HandlerGroupBase::addPreHandler(tStepHdlPtr s, tHintPtr h,
				     const HandlerGroupBase & ext) {
  if ( !s ) return;
  if ( !handler() ) refillDefaults(ext);
  thePreHandlers.push_back(make_pair(s, h));
  isEmpty = false;
}

void HandlerGroupBase::addPostHandler(tStepHdlPtr s, tHintPtr h,
				      const HandlerGroupBase & ext) {
  if ( !s ) return;
  if ( empty() ) refillDefaults(ext);
  thePostHandlers.push_back(make_pair(s, h));
  isEmpty = false;
}

void HandlerGroupBase::
addHint(tHintPtr h, const HandlerGroupBase & ext) {
  if ( !handler() || theHints.empty() ) refillDefaults(ext);
  if ( count(theHints.begin(), theHints.end(), h) ) return;
  theHints.push_back(h);
  isEmpty = false;
}

void HandlerGroupBase::clear() {
  thePreHandlers.clear();
  theHints.clear();
  thePostHandlers.clear();
  isEmpty = true;
}

void HandlerGroupBase::refillDefaults(const HandlerGroupBase & ext) {
  checkInsert(thePreHandlers, preHandlers());
  checkInsert(thePreHandlers, ext.preHandlers());
  refillDefaultHandler(ext.defaultHandler());
  checkInsert(thePostHandlers, postHandlers());
  checkInsert(thePostHandlers, ext.postHandlers());
}

void HandlerGroupBase::
checkInsert(StepHintVector & handlers, const StepVector & defHandlers) {
  for ( StepVector::const_reverse_iterator r = defHandlers.rbegin();
	r != defHandlers.rend(); ++r ) {
    try {
      for ( StepHintVector::iterator i = handlers.begin();
	    i != handlers.end(); ++i )
	if ( i->first == *r && i->second == Hint::Default() )
	  throw int();
      handlers.push_back(make_pair(*r, Hint::Default()));
      isEmpty = false;
    } catch (int) {}
  }
}

namespace
{
  void warning(tStepHdlPtr p, const HandlerGroupBase::StepVector & v)
  {
    for ( HandlerGroupBase::StepVector::const_iterator i = v.begin();
	  i != v.end(); ++i) {
      if ( p == *i )
	{
	  Throw<InterfaceException>()
	    << "\n\nWarning: Double insertion of "
	    << p->fullName() << ".\n"
	    << "         Do you intend to run the handler more than once?\n\n"
	    << Exception::warning;
	  break;
	}
    }
  }
}

void HandlerGroupBase::interfaceSetPrehandler(StepHdlPtr p, int i) {
  if ( i >= 0 && unsigned(i) < preHandlers().size() && preHandlers()[i] != p ) {
    warning( p, preHandlers() );
    preHandlers()[i] = p;
  }
}

void HandlerGroupBase::interfaceSetPosthandler(StepHdlPtr p, int i) {
  if ( i >= 0 && unsigned(i) < postHandlers().size() && postHandlers()[i] != p ) {
    warning( p, postHandlers() );
    postHandlers()[i] = p;
  }
}

void HandlerGroupBase::interfaceInsertPrehandler(StepHdlPtr p, int i) {
  if ( i >= 0 && unsigned(i) <= preHandlers().size() ) {
    warning( p, preHandlers() );
    preHandlers().insert(preHandlers().begin() + i, p);
  }
}

void HandlerGroupBase::interfaceInsertPosthandler(StepHdlPtr p, int i) {
  if ( i >= 0 && unsigned(i) <= postHandlers().size() ) {
    warning( p, postHandlers() );
    postHandlers().insert(postHandlers().begin() + i, p);
  }
}

void HandlerGroupBase::interfaceErasePrehandler(int i) {
  if ( i >= 0 && unsigned(i) < preHandlers().size() )
    preHandlers().erase(preHandlers().begin() + i);
}

void HandlerGroupBase::interfaceErasePosthandler(int i) {
  if ( i >= 0 && unsigned(i) < postHandlers().size() )
    postHandlers().erase(postHandlers().begin() + i);
}

vector<StepHdlPtr> HandlerGroupBase::interfaceGetPrehandlers() const {
  return preHandlers();
}

vector<StepHdlPtr> HandlerGroupBase::interfaceGetPosthandlers() const {
  return postHandlers();
}

void HandlerGroupBase::write(PersistentOStream & os) const {
  os << isEmpty << theDefaultPreHandlers << theDefaultPostHandlers;
  os << thePreHandlers.size();
  for ( StepHintVector::const_iterator it = thePreHandlers.begin();
	it != thePreHandlers.end(); ++it )
    os << it->first << (it->second == Hint::Default()? HintPtr(): it->second);
  os << theHints.size();
  for ( HintVector::const_iterator it = theHints.begin();
	it != theHints.end(); ++it )
    os << ( *it == Hint::Default()? HintPtr(): *it );
  os << thePostHandlers.size();
  for ( StepHintVector::const_iterator it = thePostHandlers.begin();
	it != thePostHandlers.end(); ++it )
    os << it->first << (it->second == Hint::Default()? HintPtr(): it->second);
}

void HandlerGroupBase::read(PersistentIStream & is) {
  is >> isEmpty >> theDefaultPreHandlers >> theDefaultPostHandlers;
  tStepHdlPtr sh;
  HintPtr h;
  long size;
  thePreHandlers.clear();
  is >> size;
  while ( size-- ) {
    is >> sh >> h;
    if ( !h ) h = Hint::Default();
    thePreHandlers.push_back(make_pair(sh,h));
  }
  is >> size;
  while ( size-- ) {
    is >> h;
    if ( !h ) h = Hint::Default();
    theHints.push_back(h);
  }
  is >> size;
  while ( size-- ) {
    is >> sh >> h;
    if ( !h ) h = Hint::Default();
    thePostHandlers.push_back(make_pair(sh,h));
  }
}
