// -*- C++ -*-
//
// ParMap.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ParMapBase class.
//

#include "InterfacedBase.h"
#include "ParMap.h"
#include "ParMap.xh"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
#include "ParMap.tcc"
#endif

namespace ThePEG {

string ParMapBase::
exec(InterfacedBase & i, string action, string arguments) const
  {
  istringstream arg(arguments.c_str());
  int place = 0;
  if ( !(arg >> place) ) place = -1;
  ostringstream ret;
  if ( action == "get" ) {
    StringMap v = get(i);
    if ( place >= 0 ) return v[place];
    for ( StringMap::const_iterator it = v.begin(); it != v.end(); ++it ) {
      if ( it != v.begin() ) ret << ", ";
       ret << "(" << it->first << "," << it->second << ")";
    }
  }
  else if ( action == "erase" ) {
    erase(i, place);
  }
  else if ( action == "clear" ) {
    clear(i);
  }
  else if ( action == "min" ) {
    return minimum(i, place);
  }
  else if ( action == "max" ) {
    return maximum(i, place);
  }
  else if ( action == "def" ) {
    return def(i, place);
  }
  else if ( action == "setdef" ) {
    if ( objectDefaults(i).find(tag(place)) == objectDefaults(i).end() )
      setDef(i, place);
    else
      set(i, objectDefaults(i)[tag(place)], place);
  }
  else if ( action == "set" || action == "insert" || action == "newdef") {
    string val;
    arg >> val;
    if ( action == "insert" ) insert(i, val, place);
    else set(i, val, place);
    if ( action == "newdef" ) objectDefaults(i)[tag(place)] = get(i)[place];
  }
  else if ( action == "notdef" ) {
    StringMap v = get(i);
    for ( place = 0; unsigned(place) < v.size(); ++place ) {
      string deflt = def(i, place);
      if ( objectDefaults(i).find(tag(place)) != objectDefaults(i).end() )
	deflt = objectDefaults(i)[tag(place)];
      else if ( !hasDefault ) continue;
      if ( v[place] == deflt ) continue;
      ret << "[" << place << "] " << v[place] << " (" << deflt << ") ";
    }
  }
  else
    throw InterExUnknown(*this, i);
  return ret.str();
}

string ParMapBase::fullDescription(const InterfacedBase & ib) const {
  ostringstream os;
  StringMap vals = get(ib);
  os << InterfaceBase::fullDescription(ib)
     << size() << "\n" << vals.size() << "\n";
  for ( int i = 0, N = vals.size(); i < N; ++i ) {
    string min = minimum(ib, i);
    if ( min.empty() ) min = "-inf";
    string max = maximum(ib, i);
    if ( max.empty() ) max = "inf";
    os << vals[i] << "\n"
       << min << "\n"
       << def(ib, i) << "\n"
       << max << "\n";
  }
  return os.str();
}


ParMExIndex::ParMExIndex(const InterfaceBase & i, const InterfacedBase & o,
			 int j) {
  theMessage << "Could not access element " << j
	     << " of the parameter map \"" << i.name()
	     << "\" for the object \"" << o.name()
	     << "\" because the index was outside of the allowed range.";
  severity(setuperror);
}

ParMExFixed::ParMExFixed(const InterfaceBase & i, const InterfacedBase & o) {
  theMessage << "Cannot insert or delete in the parameter map \""
	     << i.name() << "\" for the object \"" << o.name()
	     << "\" since the map is of fixed size.";
  severity(setuperror);
}

ParMExDelUnknown::ParMExDelUnknown(const InterfaceBase & i,
				   const InterfacedBase & o, int j) {
  theMessage << "Could not delete the value at position " << j
	     << " from the parameter map \"" << i.name()
	     << "\" for the object \"" << o.name()
	     << "\" because the delete function threw an unknown exception.";
  severity(setuperror);
}

ParMExGetUnknown::ParMExGetUnknown(const InterfaceBase & i,
				   const InterfacedBase & o, const char * s) {
  theMessage << "Could not get the " << s
	     << " values from the parameter map\"" << i.name()
	     << "\" for the object \"" << o.name()
	     << "\" because the get function threw an unknown exception.";
  severity(setuperror);
}

}

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
#include "ParMap.tcc"
#endif

