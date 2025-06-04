// -*- C++ -*-
//
// ParMap.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the ParMap and ParMapTBase classes.
//

#include "ParMap.xh"

namespace ThePEG {

template <typename Type>
string ParMapTBase<Type>::type() const {
  if ( std::numeric_limits<Type>::is_integer ) return "Vi";
  if ( typeid(Type) == typeid(string) ) return "Vs";
  return "Vf";
}

template <typename Type>
string ParMapTBase<Type>::doxygenType() const {
  ostringstream os;
  if ( size() <= 0 ) os << "Varying size ";
  else os << "Fixed size (" << size() << ") ";
  os << "map of ";
  string lim = "";
  if ( !limited() ) lim = " unlimited";
  if ( std::numeric_limits<Type>::is_integer ) os << lim << "integer ";
  else if ( typeid(Type) == typeid(string) ) os << "string ";
  else os << lim;
  os << "parameters";
  return os.str();
}

template <typename Type>
string ParMapTBase<Type>::fullDescription(const InterfacedBase & ib) const {
  return ParMapBase::fullDescription(ib) + def() + "\n";
}

template <typename Type>
void ParMapTBase<Type>::setDef(InterfacedBase & i, int place) const
  {
  if ( place >= 0 ) tset(i, tdef(i, place), place);
  int sz = get(i).size();
  for ( int j = 0; j < sz; ++j ) tset(i, tdef(i, j), j);
}

template <typename Type>
inline void ParMapTBase<Type>::
setImpl(InterfacedBase & i, string newValue, int place, StandardT) 
  const {
  istringstream is(newValue);
  if ( unit() > Type() ) {
    double t;
    is >> t;
    tset(i, Type(t*unit()), place);
  } else {
    Type t = Type();
    is >> t;
    tset(i, t, place);
  }
}

template<>
inline void ParMapTBase<bool>::
setImpl(InterfacedBase & i, string newValue, int place, StandardT) 
  const {
  istringstream is(newValue);
  bool t;
  is >> t;
  tset(i, t, place);
}
  
template <typename Type>
inline void ParMapTBase<Type>::
setImpl(InterfacedBase & i, string newValue, int place, DimensionT) 
  const {
  istringstream is(newValue);
  double t;
  is >> t;
  tset(i, t*unit(), place);
}

template <typename T>
void ParMapTBase<T>::
set(InterfacedBase & i, string newValue, int place) const 
  {
  setImpl(i, newValue, place, typename TypeTraits<T>::DimType());
}

template <typename Type>
inline void ParMapTBase<Type>::
insertImpl(InterfacedBase & i, string newValue, int place, StandardT) 
  const {
  istringstream is(newValue);
  if ( unit() > Type() ) {
    double t;
    is >> t;
    tinsert(i, Type(t*unit()), place);
  } else {
    Type t = Type();
    is >> t;
    tinsert(i, t, place);
  }
}

template <>
inline void ParMapTBase<bool>::
insertImpl(InterfacedBase & i, string newValue, int place, StandardT) 
  const {
  istringstream is(newValue);
  bool t;
  is >> t;
  tinsert(i, t, place);
}

template <typename Type>
inline void ParMapTBase<Type>::
insertImpl(InterfacedBase & i, string newValue, int place, DimensionT) 
  const {
  istringstream is(newValue);
  double t;
  is >> t;
  tinsert(i, t*unit(), place);
}

template <typename T>
void ParMapTBase<T>::
insert(InterfacedBase & i, string newValue, int place) const 
  {
  insertImpl(i, newValue, place, typename TypeTraits<T>::DimType());
}

template <typename Type>
typename ParMapTBase<Type>::StringMap ParMapTBase<Type>::
get(const InterfacedBase & i) const {
  TypeMap tres = tget(i);
  StringMap res;
  for ( typename TypeMap::iterator it = tres.begin();
	it != tres.end(); ++it ) {
    ostringstream os;
    putUnit(os, it->second);
    res.insert(make_pair(it->first,os.str()));
  }
  return res;
}

template <typename Type>
string ParMapTBase<Type>::
minimum(const InterfacedBase & i, int place) const {
  ostringstream os;
  putUnit(os, tminimum(i,place));
  return os.str();
}

template <typename Type>
string ParMapTBase<Type>::
maximum(const InterfacedBase & i, int place) const {
  ostringstream os;
  putUnit(os, tmaximum(i, place));
  return os.str();
}

template <typename Type>
string ParMapTBase<Type>::
def(const InterfacedBase & i, int place) const {
  ostringstream os;
  putUnit(os, tdef(i,place));
  return os.str();
}

template <typename Type>
string ParMapTBase<Type>::def() const {
  ostringstream os;
  putUnit(os, tdef());
  return os.str();
}

template <typename T, typename Type>
Type ParMap<T,Type>::tdef() const {
  return theDef;
}

template <typename T, typename Type>
void ParMap<T,Type>::tset(InterfacedBase & i, Type newValue, int place) const
  {
  if ( InterfaceBase::readOnly() ) throw InterExReadOnly(*this, i);
  T * t = dynamic_cast<T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( ( ParMapBase::lowerLimit() && newValue < tminimum(*t, place) ) ||
       ( ParMapBase::upperLimit() && newValue > tmaximum(*t, place) ) )
    throw ParMExLimit(*this, i, newValue);
  TypeMap oldMap = tget(i);
  if ( theSetFn ) {
    try { (t->*theSetFn)(newValue, place); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw ParMExUnknown(*this, i, newValue, place, "set"); }
  } else {
    if ( !theMember ) throw InterExSetup(*this, i);
    (t->*theMember)[place] = newValue;
  }
  if ( !InterfaceBase::dependencySafe() && oldMap != tget(i) ) i.touch();
}

template <typename T, typename Type>
void ParMap<T,Type>::
tinsert(InterfacedBase & i, Type newValue, int place) const
  {
  if ( InterfaceBase::readOnly() ) throw InterExReadOnly(*this, i);
  if ( ParMapBase::size() > 0 ) throw ParMExFixed(*this, i);
  T * t = dynamic_cast<T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( ( ParMapBase::lowerLimit() && newValue < tminimum(*t, place) ) ||
       ( ParMapBase::upperLimit() && newValue > tmaximum(*t, place) ) )
    throw ParMExLimit(*this, i, newValue);
  TypeMap oldMap = tget(i);
  if ( theInsFn ) {
    try { (t->*theInsFn)(newValue, place); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw ParMExUnknown(*this, i, newValue, place, "insert"); }
  } else {
    if ( !theMember ) throw InterExSetup(*this, i);
    (t->*theMember).insert(make_pair(place, newValue));
  }
  if ( !InterfaceBase::dependencySafe() && oldMap != tget(i) ) i.touch();
}

template <typename T, typename Type>
void ParMap<T,Type>::
erase(InterfacedBase & i, int place) const {
  if ( InterfaceBase::readOnly() ) throw InterExReadOnly(*this, i);
  if ( ParMapBase::size() > 0 ) throw ParMExFixed(*this, i);
  T * t = dynamic_cast<T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  TypeMap oldMap = tget(i);
  if ( theDelFn ) {
    try { (t->*theDelFn)(place); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) {  throw ParMExDelUnknown(*this, i, place); }
  } else {
    if ( !theMember ) throw InterExSetup(*this, i);
    (t->*theMember).erase(place);
  }
  if ( !InterfaceBase::dependencySafe() && oldMap != tget(i) ) i.touch();
}

template <class T, class R>
void ParMap<T,R>::clear(InterfacedBase & i) const
  {
  if ( ParMapBase::readOnly() ) throw InterExReadOnly(*this, i);
  if ( ParMapBase::size() > 0 ) throw ParMExFixed(*this, i);
  T * t = dynamic_cast<T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  (t->*theMember).clear();
  if (  !InterfaceBase::dependencySafe() ) i.touch();
}

template <typename T, typename Type>
typename ParMap<T,Type>::TypeMap ParMap<T,Type>::
tget(const InterfacedBase & i) const {
  const T * t = dynamic_cast<const T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( theGetFn ) {
    try { return (t->*theGetFn)(); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw ParMExGetUnknown(*this, i, "current"); }
  }
  if ( theMember ) return t->*theMember;
  throw InterExSetup(*this, i);
}

template <typename T, typename Type>
typename ParMap<T,Type>::StringMap ParMap<T,Type>::
get(const InterfacedBase & i) const {
  if ( !theStringGetFn ) return ParMapTBase<Type>::get(i);
  const T * t = dynamic_cast<const T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  try { return (t->*theStringGetFn)(); }
  catch (InterfaceException & e) { throw e; }
  catch ( ... ) { throw ParMExGetUnknown(*this, i, "current"); }
}


template <typename T, typename Type>
Type ParMap<T,Type>::tdef(const InterfacedBase & i, int place) const
  {
  if ( !theDefFn ) return theMin;
  const T * t = dynamic_cast<const T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  try { return (t->*theDefFn)(place); }
  catch (InterfaceException & e) { throw e; }
  catch ( ... ) { throw ParMExGetUnknown(*this, i, "default"); }
}

template <typename T, typename Type>
Type ParMap<T,Type>::tminimum(const InterfacedBase & i, int place) const
  {
  if ( !theMinFn ) return theMin;
  const T * t = dynamic_cast<const T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  try { return (t->*theMinFn)(place); }
  catch (InterfaceException & e) { throw e; }
  catch ( ... ) { throw ParMExGetUnknown(*this, i, "minimum"); }
}

template <typename T, typename Type>
Type ParMap<T,Type>::tmaximum(const InterfacedBase & i, int place) const
  {
  if ( !theMaxFn ) return theMax;
  const T * t = dynamic_cast<const T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  try { return (t->*theMaxFn)(place); }
  catch (InterfaceException & e) { throw e; }
  catch ( ... ) { throw ParMExGetUnknown(*this, i, "maximum"); }
}

template <typename T, typename Type>
void ParMap<T,Type>::doxygenDescription(ostream & os) const {
  ParMapTBase<Type>::doxygenDescription(os);
  os << "<b>Default value:</b> ";
  this->putUnit(os, theDef);
  if ( theDefFn ) os << " (May be changed by member function.)";
  if ( ParMapBase::lowerLimit() ) {
    os << "<br>\n<b>Minimum value:</b> ";
    this->putUnit(os, theMin);
    if ( theMinFn ) os << " (May be changed by member function.)";
  }
  if ( ParMapBase::upperLimit() ) {
    os << "<br>\n<b>Maximum value:</b> ";
    this->putUnit(os, theMax);
    if ( theMaxFn ) os << " (May be changed by member function.)";
  }
  os << "<br>\n";
}

namespace {
  template <typename T>
  inline
  void ostreamInsert2(ostream & os, T v, DimensionT) {
    os << ounit(v, T::baseunit());
  }
  
  template <typename T>
  inline
  void ostreamInsert2(ostream & os, T v, StandardT) {
    os << v;
  }
}

template <typename T>
ParMExLimit::ParMExLimit(const InterfaceBase & i,
			 const InterfacedBase & o, T v) {
  theMessage << "Could not set/insert ";
  ostreamInsert2(theMessage,v,typename TypeTraits<T>::DimType() );
  theMessage << " in the parameter map \""
	     << i.name() << "\" for the object \"" << o.name()
	     << "\" because the value is outside the specified limits.";
  severity(setuperror);
}

template <typename T>
ParMExUnknown::ParMExUnknown(const InterfaceBase & i, const InterfacedBase & o,
			     T v, int j, const char * s) {
  theMessage << "Could not " << s << " the value ";
  ostreamInsert2(theMessage,v,typename TypeTraits<T>::DimType() );
  theMessage << " at position "
	     << j << " in the parameter map \"" << i.name()
	     << "\" for the object \"" << o.name() << "\" because the "
	     << s << " function threw an unknown exception.";
  severity(setuperror);
}

}
