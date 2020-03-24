// -*- C++ -*-
//
// ParVector.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ParVectorBase class.
//

#include "InterfacedBase.h"
#include "ParVector.h"
#include "ParVector.xh"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
#include "ParVector.tcc"
#endif

namespace ThePEG {

string ParVectorBase::
exec(InterfacedBase & i, string action, string arguments) const
  {
  istringstream arg(arguments.c_str());
  int place = 0;
  if ( !(arg >> place) ) place = -1;
  ostringstream ret;
  if ( action == "get" ) {
    StringVector v = get(i);
    if ( place >= 0 ) return v[place];
    for ( StringVector::const_iterator it = v.begin(); it != v.end(); ++it ) {
      if ( it != v.begin() ) ret << ", ";
      ret << *it;
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
    StringVector v = get(i);
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

string ParVectorBase::fullDescription(const InterfacedBase & ib) const {
  ostringstream os;
  StringVector vals = get(ib);
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


ParVExIndex::ParVExIndex(const InterfaceBase & i, const InterfacedBase & o,
			 int j) {
  theMessage << "Could not access element " << j
	     << " of the parameter vector \"" << i.name()
	     << "\" for the object \"" << o.name()
	     << "\" because the index was outside of the allowed range.";
  severity(setuperror);
}

ParVExFixed::ParVExFixed(const InterfaceBase & i, const InterfacedBase & o) {
  theMessage << "Cannot insert or delete in the parameter vector \""
	     << i.name() << "\" for the object \"" << o.name()
	     << "\" since the vector is of fixed size.";
  severity(setuperror);
}

ParVExDelUnknown::ParVExDelUnknown(const InterfaceBase & i,
				   const InterfacedBase & o, int j) {
  theMessage << "Could not delete the value at position " << j
	     << " from the parameter vector \"" << i.name()
	     << "\" for the object \"" << o.name()
	     << "\" because the delete function threw an unknown exception.";
  severity(setuperror);
}

ParVExGetUnknown::ParVExGetUnknown(const InterfaceBase & i,
				   const InterfacedBase & o, const char * s) {
  theMessage << "Could not get the " << s
	     << " values from the parameter vector\"" << i.name()
	     << "\" for the object \"" << o.name()
	     << "\" because the get function threw an unknown exception.";
  severity(setuperror);
}

}

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
#include "ParVector.tcc"
#endif

