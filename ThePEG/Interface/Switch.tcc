// -*- C++ -*-
//
// Switch.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the Switch class.
//

namespace ThePEG {

template <typename T, typename Int>
void Switch<T,Int>::set(InterfacedBase & i, long newValue) const
  {
  T * t = dynamic_cast<T *>(&i);
  if ( readOnly() ) throw InterExReadOnly(*this, i);
  if ( !t ) throw InterExClass(*this, i);
  if ( !check(newValue) ) throw SwExSetOpt(*this, i, newValue);
  long oldValue = get(i);
  if ( theSetFn ) {
    try { (t->*theSetFn)(Int(newValue)); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw SwExSetUnknown(*this, i, newValue); }
  } else {
    if ( theMember ) t->*theMember = Int(newValue);
    else throw InterExSetup(*this, i);
  }
  if ( !InterfaceBase::dependencySafe() && oldValue != get(i) ) i.touch();
}


template <typename T, typename Int>
long Switch<T,Int>::get(const InterfacedBase & i) const
  {
  const T * t = dynamic_cast<const T *>(&i);
  if ( !t ) throw InterExClass(*this, i);
  if ( theGetFn ) {
    try { return (t->*theGetFn)(); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw SwExGetUnknown(*this, i, "current"); }
  }
  if ( theMember ) return t->*theMember;
  throw InterExSetup(*this, i);
}

template <typename T, typename Int>
long Switch<T,Int>::def(const InterfacedBase & i) const
  {
  if ( theDefFn ) {
    const T * t = dynamic_cast<const T *>(&i);
    if ( !t ) throw InterExClass(*this, i);
    try { return (t->*theDefFn)(); }
    catch (InterfaceException & e) { throw e; }
    catch ( ... ) { throw SwExGetUnknown(*this, i, "default"); }
  }
  return theDef;
}

template <typename T, typename Int>
void Switch<T,Int>::doxygenDescription(ostream & os) const {
  SwitchBase::doxygenDescription(os);
  os << "<b>Registered options:</b>\n<dl>\n";
  for ( OptionMap::const_iterator it = options().begin();
	it != options().end(); ++it )
    os << "<dt>" << it->first << "(<code>" << it->second.name() << "</code>)</dt>"
       << "<dd>" << it->second.description() << "\n";
  os << "</dl>\n<b>Default value:</b> " << theDef;
  if ( theDefFn ) os << " (May be changed by member function.)";
  os << "\n\n";
}

}
