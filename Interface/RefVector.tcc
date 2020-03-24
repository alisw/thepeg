// -*- C++ -*-
//
// RefVector.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the RefVector class.
//

namespace ThePEG {

template <class T, class R>
RefVector<T,R>::
RefVector(string newName, string newDescription,
	  Member newMember, int newSize, bool depSafe,
	  bool readonly, bool rebind, bool nullable, SetFn newSetFn,
	  InsFn newInsFn, DelFn newDelFn, GetFn newGetFn, CheckFn newCheckFn)
  : RefVectorBase(newName, newDescription, ClassTraits<T>::className(),
		  typeid(T), ClassTraits<R>::className(), typeid(R), newSize,
		  depSafe, readonly, !rebind, nullable, false),
    theMember(newMember), theSetFn(newSetFn), theInsFn(newInsFn),
    theDelFn(newDelFn), theGetFn(newGetFn), theCheckFn(newCheckFn) {}

template <class T, class R>
RefVector<T,R>::
RefVector(string newName, string newDescription,
	  Member newMember, int newSize, bool depSafe,
	  bool readonly, bool rebind, bool nullable, bool defnull,
	  SetFn newSetFn, InsFn newInsFn, DelFn newDelFn, GetFn newGetFn,
	  CheckFn newCheckFn)
  : RefVectorBase(newName, newDescription, ClassTraits<T>::className(),
		  typeid(T), ClassTraits<R>::className(), typeid(R), newSize,
		  depSafe, readonly, !rebind, nullable, defnull),
    theMember(newMember), theSetFn(newSetFn), theInsFn(newInsFn),
    theDelFn(newDelFn), theGetFn(newGetFn), theCheckFn(newCheckFn) {}

template <class T, class R>
void RefVector<T,R>::
set(InterfacedBase & i, IBPtr newRef, int place, bool chk) const
  {
  if ( readOnly() ) throw InterExReadOnly(*this, i);
  T * t = dynamic_cast<T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( noNull() && !newRef ) throw InterExNoNull(*this, i);
  RefPtr r = dynamic_ptr_cast<RefPtr>(newRef);
  if ( !r && newRef ) throw RefVExRefClass(*this, i, newRef, "set");
  IVector oldVector = get(i);
  if ( theSetFn && ( chk || !theMember ) ) {
    try { (t->*theSetFn)(r, place); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw RefVExSetUnknown(*this, i, r, place, "set"); }
  } else {
    if ( !theMember ) throw RefVExNoSet(*this, i);
    if ( place < 0 || static_cast<unsigned long>(place) >= (t->*theMember).size() )
      throw RefVExIndex(*this, i, place);
    (t->*theMember)[place] = r;
  }
  if ( !InterfaceBase::dependencySafe() && oldVector != get(i) ) i.touch();
}

template <class T, class R>
void RefVector<T,R>::
insert(InterfacedBase & i, IBPtr newRef, int place, bool chk) const
  {
  if ( readOnly() ) throw InterExReadOnly(*this, i);
  if ( size() > 0 ) throw RefVExFixed(*this, i);
  T * t = dynamic_cast<T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( noNull() && !newRef ) throw InterExNoNull(*this, i);
  RefPtr r = dynamic_ptr_cast<RefPtr>(newRef);
  if ( !r && newRef ) throw RefVExRefClass(*this, i, newRef, "insert");
  IVector oldVector = get(i);
  if ( theInsFn && ( chk || !theMember ) ) {
    try { (t->*theInsFn)(r, place); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw RefVExSetUnknown(*this, i, r, place, "insert"); }
  } else {
    if ( !theMember ) throw RefVExNoIns(*this, i);
    if ( place < 0 || static_cast<unsigned long>(place) > (t->*theMember).size() )
      throw RefVExIndex(*this, i, place);
    (t->*theMember).insert((t->*theMember).begin()+place, r);
  }
  if ( !InterfaceBase::dependencySafe() && oldVector != get(i) ) i.touch();
}

template <class T, class R>
void RefVector<T,R>::erase(InterfacedBase & i, int place) const
  {
  if ( readOnly() ) throw InterExReadOnly(*this, i);
  if ( size() > 0 ) throw RefVExFixed(*this, i);
  T * t = dynamic_cast<T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  IVector oldVector = get(i);
  if ( theDelFn ) {
    try { (t->*theDelFn)(place); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw RefVExDelUnknown(*this, i, place); }
  } else {
    if ( !theMember ) throw RefVExNoDel(*this, i);
    if ( place < 0 || static_cast<unsigned long>(place) >= (t->*theMember).size() )
      throw RefVExIndex(*this, i, place);
    (t->*theMember).erase((t->*theMember).begin()+place);
  }
  if (  !InterfaceBase::dependencySafe() && oldVector != get(i) ) i.touch();
}

template <class T, class R>
void RefVector<T,R>::clear(InterfacedBase & i) const
  {
  if ( readOnly() ) throw InterExReadOnly(*this, i);
  if ( size() > 0 ) throw RefVExFixed(*this, i);
  T * t = dynamic_cast<T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( !theMember ) throw RefVExNoDel(*this, i);
  (t->*theMember).clear();
  if (  !InterfaceBase::dependencySafe() ) i.touch();
}

template <class T, class R>
IVector RefVector<T,R>::get(const InterfacedBase & i) const
  {
  const T * t = dynamic_cast<const T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( theGetFn ) {
    try {
      vector<RefPtr> ret = (t->*theGetFn)();
      return IVector(ret.begin(), ret.end());
    }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw RefVExGetUnknown(*this, i); }
  }
  if ( theMember ) return IVector((t->*theMember).begin(),
				  (t->*theMember).end());
  throw InterExSetup(*this, i);
}

template <class T, class R>
bool RefVector<T,R>::check(const InterfacedBase & i, cIBPtr ir, int place) const
  {
  const T * t = dynamic_cast<const T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( noNull() && !ir ) return false;
  cRefPtr r = dynamic_ptr_cast<cRefPtr>(ir);
  if ( !r && ir ) return false;
  if ( theCheckFn ) return (t->*theCheckFn)(r, place);
  if ( !theMember ) return true;
  return place >= 0 && static_cast<unsigned long>(place) <= (t->*theMember).size(); 
}

}

