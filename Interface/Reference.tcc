// -*- C++ -*-
//
// Reference.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the Reference class.
//

namespace ThePEG {

template <class T, class R>
void Reference<T,R>::set(InterfacedBase & i, IBPtr newRef, bool chk) const
  {
  if ( readOnly() ) throw InterExReadOnly(*this, i);
  T * t = dynamic_cast<T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( noNull() && !newRef ) throw InterExNoNull(*this, i);
  RefPtr r = dynamic_ptr_cast<RefPtr>(newRef);
  if ( !r && newRef) throw RefExSetRefClass(*this, i, newRef);
  RefPtr oldRef = dynamic_ptr_cast<RefPtr>(get(i));
  if ( theSetFn && ( chk || !theMember ) ) {
    try { (t->*theSetFn)(r); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw RefExSetUnknown(*this, i, r); }
  } else {
    if ( theMember ) t->*theMember = r;
    else throw InterExSetup(*this, i);
  }
  if ( !InterfaceBase::dependencySafe() && oldRef != get(i) ) i.touch();  
}

template <class T, class R>
IBPtr Reference<T,R>::get(const InterfacedBase & i) const
  {
  const T * t = dynamic_cast<const T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( theGetFn ) {
    try { return (t->*theGetFn)(); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw RefExGetUnknown(*this, i); }
  }
  if ( theMember ) return t->*theMember;
  throw InterExSetup(*this, i);
}

template <class T, class R>
bool Reference<T,R>::check(const InterfacedBase & i, cIBPtr ir) const
  {
  const T * t = dynamic_cast<const T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( noNull() && !ir ) return false;
  cRefPtr r = dynamic_ptr_cast<cRefPtr>(ir);
  if ( !r && ir ) return false;
  if ( !theCheckFn ) return true;
  return (t->*theCheckFn)(r);
}

}

