// -*- C++ -*-
//
// Parameter.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the Parameter and ParameterTBase classes.
//

namespace ThePEG {

template <typename Type>
string ParameterTBase<Type>::type() const {
  if ( std::numeric_limits<Type>::is_integer ) return "Pi";
  if ( typeid(Type) == typeid(string) ) return "Ps";
  return "Pf";
}

template <typename Type>
string ParameterTBase<Type>::doxygenType() const {
  string lim = "";
  if ( !limited() ) lim = "Unlimited ";
  if ( std::numeric_limits<Type>::is_integer ) return lim + "Integer parameter";
  if ( typeid(Type) == typeid(string) ) return "Character string parameter";
  return lim + "Parameter";
}

template <typename Type>
inline void 
ParameterTBase<Type>::setImpl(InterfacedBase & i, 
			      string newValue, StandardT) 
  const {
  istringstream is(newValue);
  if ( unit() > Type() ) {
    double t;
    is >> t;
    tset(i, umult(t, unit()));
  } else {
    Type t = Type();
    is >> t;
    tset(i, t);
  }
}

template <typename Type>
inline void 
ParameterTBase<Type>::setImpl(InterfacedBase & i, 
			      string newValue, DimensionT) 
  const {
  istringstream is(newValue);
  double t;
  is >> t;
  // if 'is' has no more chars, all stream ops below are no-ops
  is.ignore(); // skip the connecting char
  string suffix;
  is >> suffix;
  checkUnitConsistency(suffix);
  tset(i, umult(t, unit()));
}
  
template <typename Type>
inline void 
ParameterTBase<Type>::setImpl(InterfacedBase & i, 
			      string newValue, EnumT) 
  const {
  istringstream is(newValue);
  int t;
  is >> t;
  tset(i, Type(t));
}
  
    // Macs need a visible template specialization.
template <>
void ParameterTBase<Energy>::
checkUnitConsistency(string suffix) const;
  
template <>
void ParameterTBase<Energy2>::
checkUnitConsistency(string suffix) const;
  
template <>
void ParameterTBase<Length>::
checkUnitConsistency(string suffix) const;
  

template <typename T>
void ParameterTBase<T>::
checkUnitConsistency(string suffix) const {
  if ( ! suffix.empty() ) {
     Throw<InterfaceException>()
       << name() 
       << ": unit suffix " << suffix << " will be ignored.\n"
       << "The unit specified in the parameter definition is used instead.\n\n"
       << "To proceed, remove the unit suffix in the input file or \n"
       << "request unit support for " << suffix << " to be added.\n\n"
       << Exception::setuperror;
  }
}
  
  

template <typename T>
void ParameterTBase<T>::
set(InterfacedBase & i, string newValue) const {
  setImpl(i, newValue, typename TypeTraits<T>::DimType());
}

template <typename Type>
string ParameterTBase<Type>::
get(const InterfacedBase & i) const {
  ostringstream os;
  putUnit(os, tget(i));
  return os.str();
}

template <typename Type>
string ParameterTBase<Type>::
minimum(const InterfacedBase & i) const {
  ostringstream os;
  if ( ParameterBase::lowerLimit() ) putUnit(os, tminimum(i));
  return os.str();
}

template <typename Type>
string ParameterTBase<Type>::
maximum(const InterfacedBase & i) const {
  ostringstream os;
  if ( ParameterBase::upperLimit() ) putUnit(os, tmaximum(i));
  return os.str();
}
    
template <typename Type>
string ParameterTBase<Type>::
def(const InterfacedBase & i) const {
  ostringstream os;
  putUnit(os, tdef(i));
  return os.str();
}

template <typename T, typename Type>
void Parameter<T,Type>::tset(InterfacedBase & i, Type newValue) const
  {
  if ( InterfaceBase::readOnly() ) throw InterExReadOnly(*this, i);
  T * t = dynamic_cast<T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( ( ParameterBase::lowerLimit() && newValue < tminimum(i) ) ||
       ( ParameterBase::upperLimit() && newValue > tmaximum(i) ) )
    throw ParExSetLimit(*this, i, newValue);
  Type oldValue = tget(i);
  if ( theSetFn ) {
    try { (t->*theSetFn)(newValue); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw ParExSetUnknown(*this, i, newValue); }
  } else {
    if ( theMember ) t->*theMember = newValue;
    else throw InterExSetup(*this, i);
  }
  if ( !InterfaceBase::dependencySafe() && oldValue != tget(i)) i.touch();
}

template <typename T>
void Parameter<T,string>::tset(InterfacedBase & i, string newValue) const
  {
  if ( InterfaceBase::readOnly() ) throw InterExReadOnly(*this, i);
  T * t = dynamic_cast<T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  string oldValue = tget(i);
  if ( theSetFn ) {
    try { (t->*theSetFn)(newValue); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw ParExSetUnknown(*this, i, newValue); }
  } else {
    if ( theMember ) t->*theMember = newValue;
    else throw InterExSetup(*this, i);
  }
  if ( !InterfaceBase::dependencySafe() && oldValue != tget(i)) i.touch();
}

template <typename T, typename Type>
Type Parameter<T,Type>::tget(const InterfacedBase & i) const
  {
  const T * t = dynamic_cast<const T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( theGetFn ) {
    try { return (t->*theGetFn)(); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw ParExGetUnknown(*this, i, "current"); }
  }
  if ( theMember ) return t->*theMember;
  else throw InterExSetup(*this, i);
}

template <typename T>
string Parameter<T,string>::tget(const InterfacedBase & i) const
  {
  const T * t = dynamic_cast<const T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( theGetFn ) {
    try { return (t->*theGetFn)(); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw ParExGetUnknown(*this, i, "current"); }
  }
  if ( theMember ) return t->*theMember;
  else throw InterExSetup(*this, i);
}

template <typename T, typename Type>
Type Parameter<T,Type>::tminimum(const InterfacedBase & i) const
  {
  if ( theMinFn ) {
    const T * t = dynamic_cast<const T *>(&i);
    if ( !t ) throw InterExClass(*this, i);
    try { return max(theMin, (t->*theMinFn)()); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw ParExGetUnknown(*this, i, "minimum"); }
  }
  return theMin;
}

template <typename T, typename Type>
Type Parameter<T,Type>::tmaximum(const InterfacedBase & i) const
  {
  if ( theMaxFn ) {
    const T * t = dynamic_cast<const T *>(&i);
    if ( !t ) throw InterExClass(*this, i);
    try { return min(theMax, (t->*theMaxFn)()); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw ParExGetUnknown(*this, i, "maximum"); }
  }
  return theMax;
}

template <typename T, typename Type>
Type Parameter<T,Type>::tdef(const InterfacedBase & i) const
  {
  if ( theDefFn ) {
    const T * t = dynamic_cast<const T *>(&i);
    if ( !t ) throw InterExClass(*this, i);
    try { return (t->*theDefFn)(); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw ParExGetUnknown(*this, i, "default"); }
  }
  return theDef;
}

template <typename T>
string Parameter<T,string>::tdef(const InterfacedBase & i) const
  {
  if ( theDefFn ) {
    const T * t = dynamic_cast<const T *>(&i);
    if ( !t ) throw InterExClass(*this, i);
    try { return (t->*theDefFn)(); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw ParExGetUnknown(*this, i, "default"); }
  }
  return theDef;
}

template <typename T, typename Type>
void Parameter<T,Type>::doxygenDescription(ostream & os) const {
  ParameterTBase<Type>::doxygenDescription(os);
  os << "<b>Default value:</b> ";
  this->putUnit(os, theDef);
  if ( theDefFn ) os << " (May be changed by member function.)";
  if ( ParameterBase::lowerLimit() ) {
    os << "<br>\n<b>Minimum value:</b> ";
    this->putUnit(os, theMin);
    if ( theMinFn ) os << " (May be changed by member function.)";
  }
  if ( ParameterBase::upperLimit() ) {
    os << "<br>\n<b>Maximum value:</b> ";
    this->putUnit(os, theMax);
    if ( theMaxFn ) os << " (May be changed by member function.)";
  }
  os << "<br>\n";
}

template <typename T>
void Parameter<T,string>::doxygenDescription(ostream & os) const {
  ParameterTBase<string>::doxygenDescription(os);
  os << "<b>Default value:</b> " << theDef;
  if ( theDefFn ) os << " (May be changed by member function.)";
  os << "<br>\n";
}

namespace {
  template <typename T>
  inline
  void ostreamInsert(ostream & os, T v, DimensionT) {
    os << ounit(v, T::baseunit());
  }
  
  template <typename T>
  inline
  void ostreamInsert(ostream & os, T v, StandardT) {
    os << v;
  }
  
  template <typename T>
  inline
  void ostreamInsert(ostream & os, T v, EnumT) {
    os << v;
  }
}

template <typename T>
ParExSetLimit::ParExSetLimit(const InterfaceBase & i,
			     const InterfacedBase & o, T v) {
  theMessage << "Could not set the parameter \"" << i.name()
	     << "\" for the object \"" << o.name() << "\" to ";
  ostreamInsert(theMessage,v,typename TypeTraits<T>::DimType() );
  theMessage << " because the value is outside the specified limits.";
  severity(setuperror);
}

template <typename T>
ParExSetUnknown::ParExSetUnknown(const InterfaceBase & i,
				 const InterfacedBase & o, T v) {
  theMessage << "Could not set the parameter \"" << i.name()
	     << "\" for the object \"" << o.name() << "\" to ";
  ostreamInsert(theMessage,v,typename TypeTraits<T>::DimType() );
  theMessage << " because the set function threw an unknown exception.";
  severity(setuperror);
}

}

