// -*- C++ -*-
//
// Parameter.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ParameterBase class.
//

#include "InterfacedBase.h"
#include "Parameter.h"
#include "Parameter.xh"

namespace ThePEG {

ParameterBase::~ParameterBase() {}

string ParameterBase::exec(InterfacedBase & i, string action,
			   string arguments) const  {
  if ( action == "get" ) {
    return get(i);
  }
  else if ( action == "min" ) {
    return minimum(i);
  }
  else if ( action == "max" ) {
    return maximum(i);
  }
  else if ( action == "def" ) {
    return def(i);
  }
  else if ( action == "setdef" ) {
    if ( objectDefaults(i).find(name()) == objectDefaults(i).end() )
      setDef(i);
    else
      set(i, objectDefaults(i)[name()]);
  }
  else if ( action == "set" || action == "newdef" ) {
    set(i, arguments);
    if ( action == "newdef" ) objectDefaults(i)[name()] = get(i);
  } else if ( action == "notdef" ) {
    string deflt = def(i);
    if ( objectDefaults(i).find(name()) != objectDefaults(i).end() )
      deflt = objectDefaults(i)[name()];
    else if ( !hasDefault ) return "";
    if ( deflt != get(i) ) return get(i) + " (" + deflt + ")";
  }
  else
    throw InterExUnknown(*this, i);
  return "";
}

string ParameterBase::fullDescription(const InterfacedBase & ib) const {
  string min = minimum(ib);
  if ( min.empty() ) min = "-inf";
  string max = maximum(ib);
  if ( max.empty() ) max = "inf";
  
  return InterfaceBase::fullDescription(ib) + get(ib) + '\n' +
    min + '\n' + def(ib) + '\n' + max + '\n';
}

ParExGetUnknown::ParExGetUnknown(const InterfaceBase & i,
				 const InterfacedBase & o, const char * s) {
  theMessage << "Could not get the " << s << " value of parameter \""
	     << i.name() << "\" for the object \"" << o.name()
	     << "\" because the get function threw an unknown exception.";
  severity(setuperror);
}

}

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
#include "Parameter.tcc"
#endif
