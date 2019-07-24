// -*- C++ -*-
//
// ParVector.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the ParVector and ParVectorTBase classes.
//

#include "ParVector.xh"

namespace ThePEG {

template <typename Type>
string ParVectorTBase<Type>::type() const {
  if ( std::numeric_limits<Type>::is_integer ) return "Vi";
  if ( typeid(Type) == typeid(string) ) return "Vs";
  return "Vf";
}

template <typename Type>
string ParVectorTBase<Type>::doxygenType() const {
  ostringstream os;
  if ( size() <= 0 ) os << "Varying size ";
  else os << "Fixed size (" << size() << ") ";
  os << "vector of ";
  string lim = "";
  if ( !limited() ) lim = " unlimited";
  if ( std::numeric_limits<Type>::is_integer ) os << lim << "integer ";
  else if ( typeid(Type) == typeid(string) ) os << "string ";
  else os << lim;
  os << "parameters";
  return os.str();
}

template <typename Type>
string ParVectorTBase<Type>::fullDescription(const InterfacedBase & ib) const {
  return ParVectorBase::fullDescription(ib) + def() + "\n";
}

template <typename Type>
void ParVectorTBase<Type>::setDef(InterfacedBase & i, int place) const
  {
  if ( place >= 0 ) tset(i, tdef(i, place), place);
  int sz = get(i).size();
  for ( int j = 0; j < sz; ++j ) tset(i, tdef(i, j), j);
}

template <typename Type>
inline void ParVectorTBase<Type>::
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
inline void ParVectorTBase<bool>::
setImpl(InterfacedBase & i, string newValue, int place, StandardT) 
  const {
  istringstream is(newValue);
  bool t;
  is >> t;
  tset(i, t, place);
}
  
template <typename Type>
inline void ParVectorTBase<Type>::
setImpl(InterfacedBase & i, string newValue, int place, DimensionT) 
  const {
  istringstream is(newValue);
  double t;
  is >> t;
  tset(i, t*unit(), place);
}

template <typename T>
void ParVectorTBase<T>::
set(InterfacedBase & i, string newValue, int place) const 
  {
  setImpl(i, newValue, place, typename TypeTraits<T>::DimType());
}

template <typename Type>
inline void ParVectorTBase<Type>::
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
inline void ParVectorTBase<bool>::
insertImpl(InterfacedBase & i, string newValue, int place, StandardT) 
  const {
  istringstream is(newValue);
  bool t;
  is >> t;
  tinsert(i, t, place);
}

template <typename Type>
inline void ParVectorTBase<Type>::
insertImpl(InterfacedBase & i, string newValue, int place, DimensionT) 
  const {
  istringstream is(newValue);
  double t;
  is >> t;
  tinsert(i, t*unit(), place);
}

template <typename T>
void ParVectorTBase<T>::
insert(InterfacedBase & i, string newValue, int place) const 
  {
  insertImpl(i, newValue, place, typename TypeTraits<T>::DimType());
}

template <typename Type>
typename ParVectorTBase<Type>::StringVector ParVectorTBase<Type>::
get(const InterfacedBase & i) const {
  TypeVector tres = tget(i);
  StringVector res;
  for ( typename TypeVector::iterator it = tres.begin();
	it != tres.end(); ++it ) {
    ostringstream os;
    putUnit(os, *it);
    res.push_back(os.str());
  }
  return res;
}

template <typename Type>
string ParVectorTBase<Type>::
minimum(const InterfacedBase & i, int place) const {
  ostringstream os;
  putUnit(os, tminimum(i,place));
  return os.str();
}

template <typename Type>
string ParVectorTBase<Type>::
maximum(const InterfacedBase & i, int place) const {
  ostringstream os;
  putUnit(os, tmaximum(i, place));
  return os.str();
}

template <typename Type>
string ParVectorTBase<Type>::
def(const InterfacedBase & i, int place) const {
  ostringstream os;
  putUnit(os, tdef(i,place));
  return os.str();
}

template <typename Type>
string ParVectorTBase<Type>::def() const {
  ostringstream os;
  putUnit(os, tdef());
  return os.str();
}

template <typename T, typename Type>
Type ParVector<T,Type>::tdef() const {
  return theDef;
}

template <typename T, typename Type>
void ParVector<T,Type>::tset(InterfacedBase & i, Type newValue, int place) const
  {
  if ( InterfaceBase::readOnly() ) throw InterExReadOnly(*this, i);
  T * t = dynamic_cast<T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( ( ParVectorBase::lowerLimit() && newValue < tminimum(*t, place) ) ||
       ( ParVectorBase::upperLimit() && newValue > tmaximum(*t, place) ) )
    throw ParVExLimit(*this, i, newValue);
  TypeVector oldVector = tget(i);
  if ( theSetFn ) {
    try { (t->*theSetFn)(newValue, place); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw ParVExUnknown(*this, i, newValue, place, "set"); }
  } else {
    if ( !theMember ) throw InterExSetup(*this, i);
    if ( place < 0 || unsigned(place) >= (t->*theMember).size() )
      throw ParVExIndex(*this, i, place);
    (t->*theMember)[place] = newValue;
  }
  if ( !InterfaceBase::dependencySafe() && oldVector != tget(i) ) i.touch();
}

template <typename T, typename Type>
void ParVector<T,Type>::
tinsert(InterfacedBase & i, Type newValue, int place) const
  {
  if ( InterfaceBase::readOnly() ) throw InterExReadOnly(*this, i);
  if ( ParVectorBase::size() > 0 ) throw ParVExFixed(*this, i);
  T * t = dynamic_cast<T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( ( ParVectorBase::lowerLimit() && newValue < tminimum(*t, place) ) ||
       ( ParVectorBase::upperLimit() && newValue > tmaximum(*t, place) ) )
    throw ParVExLimit(*this, i, newValue);
  TypeVector oldVector = tget(i);
  if ( theInsFn ) {
    try { (t->*theInsFn)(newValue, place); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw ParVExUnknown(*this, i, newValue, place, "insert"); }
  } else {
    if ( !theMember ) throw InterExSetup(*this, i);
    if ( place < 0 || unsigned(place) > (t->*theMember).size() )
      throw ParVExIndex(*this, i, place);
    (t->*theMember).insert((t->*theMember).begin()+place, newValue);
  }
  if ( !InterfaceBase::dependencySafe() && oldVector != tget(i) ) i.touch();
}

template <typename T, typename Type>
void ParVector<T,Type>::
erase(InterfacedBase & i, int place) const {
  if ( InterfaceBase::readOnly() ) throw InterExReadOnly(*this, i);
  if ( ParVectorBase::size() > 0 ) throw ParVExFixed(*this, i);
  T * t = dynamic_cast<T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  TypeVector oldVector = tget(i);
  if ( theDelFn ) {
    try { (t->*theDelFn)(place); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) {  throw ParVExDelUnknown(*this, i, place); }
  } else {
    if ( !theMember ) throw InterExSetup(*this, i);
    if ( place < 0 || unsigned(place) >= (t->*theMember).size() )
      throw ParVExIndex(*this, i, place);
    (t->*theMember).erase((t->*theMember).begin()+place);
  }
  if ( !InterfaceBase::dependencySafe() && oldVector != tget(i) ) i.touch();
}

template <class T, class R>
void ParVector<T,R>::clear(InterfacedBase & i) const
  {
  if ( ParVectorBase::readOnly() ) throw InterExReadOnly(*this, i);
  if ( ParVectorBase::size() > 0 ) throw ParVExFixed(*this, i);
  T * t = dynamic_cast<T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  (t->*theMember).clear();
  if (  !InterfaceBase::dependencySafe() ) i.touch();
}

template <typename T, typename Type>
typename ParVector<T,Type>::TypeVector ParVector<T,Type>::
tget(const InterfacedBase & i) const {
  const T * t = dynamic_cast<const T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( theGetFn ) {
    try { return (t->*theGetFn)(); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw ParVExGetUnknown(*this, i, "current"); }
  }
  if ( theMember ) return t->*theMember;
  throw InterExSetup(*this, i);
}

template <typename T, typename Type>
typename ParVector<T,Type>::StringVector ParVector<T,Type>::
get(const InterfacedBase & i) const {
  if ( !theStringGetFn ) return ParVectorTBase<Type>::get(i);
  const T * t = dynamic_cast<const T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  try { return (t->*theStringGetFn)(); }
  catch (InterfaceException & e) { throw e; }
  catch ( ... ) { throw ParVExGetUnknown(*this, i, "current"); }
}


template <typename T, typename Type>
Type ParVector<T,Type>::tdef(const InterfacedBase & i, int place) const
  {
  if ( place < 0 || !theDefFn ) return theMin;
  const T * t = dynamic_cast<const T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  try { return (t->*theDefFn)(place); }
  catch (InterfaceException & e) { throw e; }
  catch ( ... ) { throw ParVExGetUnknown(*this, i, "default"); }
}

template <typename T, typename Type>
Type ParVector<T,Type>::tminimum(const InterfacedBase & i, int place) const
  {
  if ( place < 0 || !theMinFn ) return theMin;
  const T * t = dynamic_cast<const T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  try { return (t->*theMinFn)(place); }
  catch (InterfaceException & e) { throw e; }
  catch ( ... ) { throw ParVExGetUnknown(*this, i, "minimum"); }
}

template <typename T, typename Type>
Type ParVector<T,Type>::tmaximum(const InterfacedBase & i, int place) const
  {
  if ( place < 0 || !theMaxFn ) return theMax;
  const T * t = dynamic_cast<const T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  try { return (t->*theMaxFn)(place); }
  catch (InterfaceException & e) { throw e; }
  catch ( ... ) { throw ParVExGetUnknown(*this, i, "maximum"); }
}

template <typename T, typename Type>
void ParVector<T,Type>::doxygenDescription(ostream & os) const {
  ParVectorTBase<Type>::doxygenDescription(os);
  os << "<b>Default value:</b> ";
  this->putUnit(os, theDef);
  if ( theDefFn ) os << " (May be changed by member function.)";
  if ( ParVectorBase::lowerLimit() ) {
    os << "<br>\n<b>Minimum value:</b> ";
    this->putUnit(os, theMin);
    if ( theMinFn ) os << " (May be changed by member function.)";
  }
  if ( ParVectorBase::upperLimit() ) {
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
ParVExLimit::ParVExLimit(const InterfaceBase & i,
			 const InterfacedBase & o, T v) {
  theMessage << "Could not set/insert ";
  ostreamInsert2(theMessage,v,typename TypeTraits<T>::DimType() );
  theMessage << " in the parameter vector \""
	     << i.name() << "\" for the object \"" << o.name()
	     << "\" because the value is outside the specified limits.";
  severity(setuperror);
}

template <typename T>
ParVExUnknown::ParVExUnknown(const InterfaceBase & i, const InterfacedBase & o,
			     T v, int j, const char * s) {
  theMessage << "Could not " << s << " the value ";
  ostreamInsert2(theMessage,v,typename TypeTraits<T>::DimType() );
  theMessage << " at position "
	     << j << " in the parameter vector \"" << i.name()
	     << "\" for the object \"" << o.name() << "\" because the "
	     << s << " function threw an unknown exception.";
  severity(setuperror);
}

}
