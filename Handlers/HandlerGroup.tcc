// -*- C++ -*-
//
// HandlerGroup.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the HandlerGroup class.
//

#include "ThePEG/Utilities/ClassTraits.h"
#include "ThePEG/Handlers/Hint.h"

namespace ThePEG {

template <typename HDLR>
HandlerGroup<HDLR>::~HandlerGroup() {}

template <typename HDLR>
bool HandlerGroup<HDLR>::setHandler(tStepHdlPtr hin,
				    const HandlerGroupBase & ext) {
  tHdlPtr h = dynamic_ptr_cast<HdlPtr>(hin);
  if ( !h ) return false;
  if ( !theHandler ) refillDefaults(ext);
  theHandler = h;
  isEmpty = false;
  return true;
}

template <typename HDLR>
void HandlerGroup<HDLR>::refillDefaultHandler(tStepHdlPtr h) {
  tHdlPtr ext = dynamic_ptr_cast<HdlPtr>(h);
  if ( ext ) theHandler = ext;
  else theHandler = theDefaultHandler;
  if ( theHandler ) {
    for ( int i = 0, N = theHints.size(); i < N; ++i )
      if ( theHints[i] == Hint::Default() ) return;
    theHints.push_front(Hint::Default());
    isEmpty = false;
  }
}

template <typename HDLR>
void HandlerGroup<HDLR>::clear() {
  theHandler = HdlPtr();
  HandlerGroupBase::clear();
}

template <typename HDLR>
string HandlerGroup<HDLR>::handlerClass() const {
  return ClassTraits<HDLR>::className();
}

template <typename HDLR>
void HandlerGroup<HDLR>::interfaceSetHandler(HdlPtr p) {
  theDefaultHandler = p;
}

template <typename HDLR>
typename HandlerGroup<HDLR>::HdlPtr
HandlerGroup<HDLR>::interfaceGetHandler() const {
  return theDefaultHandler;
}

}
