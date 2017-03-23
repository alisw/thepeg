// -*- C++ -*-
//
// RefVector.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RefVector class.
//

#include "InterfacedBase.h"
#include "RefVector.h"
#include "RefVector.xh"
#include "ThePEG/Utilities/HoldFlag.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Repository/BaseRepository.h"
#include "ThePEG/Repository/EventGenerator.h"

using namespace ThePEG;

RefVectorBase::
RefVectorBase(string newName, string newDescription,
	      string newClassName, const type_info & newTypeInfo,
	      string newRefClassName, const type_info & newRefTypeInfo,
	      int newSize, bool depSafe, bool readonly,
	      bool norebind, bool nullable, bool defnull)
  : RefInterfaceBase(newName, newDescription, newClassName, newTypeInfo,
		     newRefClassName, newRefTypeInfo, depSafe,
		     readonly, norebind, nullable, defnull),
    theSize(newSize) {}

IVector RefVectorBase::getReferences(const InterfacedBase & i) const {
  return get(i);
}

string RefVectorBase::exec(InterfacedBase & i, string action,
			   string arguments) const {
  istringstream arg(arguments.c_str());
  int place = -1;
  if ( !( arg >> place ) ) place = -1;
  ostringstream ret;
  if ( action == "get" ) {
    IVector refvec = get(i);
    if ( place >= 0 && unsigned(place) < refvec.size() ) {
      if ( refvec[place] ) return refvec[place]->fullName();
      else return "*** NULL Reference ***";
    }
    for ( IVector::size_type j = 0; j < refvec.size(); ++j ) {
      if ( j != 0 ) ret << ", ";
      if ( refvec[j] ) ret << refvec[j]->fullName();
      else ret << "*** NULL Reference ***";
    }
  }
  else if ( action == "erase" ) {
    erase(i, place);
  }
  else if ( action == "clear" ) {
    clear(i);
  }
  else if ( action == "set" || action == "insert" ||
	      action == "newdef" || action == "setdef" ) {
    string refname;
    arg >> refname;
    if ( action == "setdef" ) {
      if ( objectDefaults(i).find(tag(place)) == objectDefaults(i).end() )
	return "Error: No default value defined for this object.";
      refname = objectDefaults(i)[tag(place)];
    }
    IBPtr ip;
    if ( refname.size() && refname != "NULL") {
      Interfaced * ii = dynamic_cast<Interfaced *>(&i);
      if ( ii && ii->generator() )
	ip = ii->generator()->getObject<Interfaced>(refname);
      else
	ip = BaseRepository::TraceObject(refname);
    }
    if ( action == "insert" ) insert(i, ip, place);
    else set(i, ip, place);
    if ( action == "newdef" ) {
      IVector refvec = get(i);
      if ( place >= 0 && unsigned(place) < refvec.size() )
	objectDefaults(i)[tag(place)] =
	  refvec[place]? refvec[place]->fullName(): string("NULL");
    }
  }
  else if ( action == "notdef" ) {
    IVector refvec = get(i);
    for ( place = 0; unsigned(place) < refvec.size(); ++place ) {
      if ( objectDefaults(i).find(tag(place)) == objectDefaults(i).end() )
	continue;
      string refname = refvec[place]? refvec[place]->fullName(): string("NULL");
      if ( refname == objectDefaults(i)[tag(place)] ) continue;
      ret << "[" << place << "] " << refname
	  << " (" << objectDefaults(i)[tag(place)] << ") ";
    }
  }
  else
    throw InterExUnknown(*this, i);
  return ret.str();
}

string RefVectorBase::fullDescription(const InterfacedBase & ib) const {
  ostringstream os;
  os << InterfaceBase::fullDescription(ib)
     << ( noNull()? "nevernull\n": "nullable\n" )
     << ( defaultIfNull()? "defnull\n": "nodefnull\n" );
  IVector refs = get(ib);
  os << size() << '\n' << refs.size() << '\n';
  for ( int i = 0, N = refs.size(); i < N; ++i )
  if ( !refs[i] ) os << "NULL\n";
  else os << refs[i]->fullName() << '\n';
  return os.str();
}

string RefVectorBase::type() const {
  return string("V<") + refClassName() + ">";
}

string RefVectorBase::doxygenType() const {
  ostringstream os;
  if ( size() <= 0 ) os << "Varying size ";
  else os << "Fixed size (" << size() << ") ";
  os << "vector of references to objects of class " << refClassName();
  return os.str();
}

void RefVectorBase::
rebind(InterfacedBase & i, const TranslationMap & trans,
       const IVector & defs) const {
  if ( noRebind() ) return;
  IVector oldrefs =  get(i);
  IVector newRefs;
  for ( unsigned int ir = 0; ir < oldrefs.size(); ++ir ) {
    IBPtr oldref = oldrefs[ir];;
    IBPtr newref;
    if ( oldref ) {
      newref = trans.translate(oldref);
      if ( !dependencySafe() && newref->fullName() != oldref->fullName() )
	i.touch();
    }
    else if ( defaultIfNull() ) {
      for ( IVector::const_iterator p = defs.begin(); p != defs.end(); ++p ) {
	try {
	  if ( *p && check(i, *p, ir) ) {
	    newref = *p;
	    i.touch();
	    break;
	  }
	} catch ( ... ) {}
      }
    }
    newRefs.push_back(newref);
  }

  HoldFlag<> depflag(isDependencySafe);
  HoldFlag<> roflag(isReadOnly, false);

  if ( size() <= 0 )
    for ( IVector::size_type j = oldrefs.size(); j > 0; --j )
      erase(i, j - 1);

  for ( IVector::size_type j = 0; j < oldrefs.size(); ++j )
    if ( size() > 0 )
      set(i, newRefs[j], j, false);
    else
      insert(i, newRefs[j], j, false);
}

RefVExRefClass::RefVExRefClass(const RefInterfaceBase & i,
			       const InterfacedBase & o,
			       cIBPtr r, const char * s) {
  theMessage << "Could not " << s << " the object \""
	     << (r? r->name().c_str(): "<NULL>")
	     << "\" in the reference vector \""
	     << i.name() << "\" for the object \"" << o.name()
	     << "\" because it is not of the required class ("
	     << i.refClassName() << ").";
  severity(setuperror);
}

RefVExSetUnknown::RefVExSetUnknown(const RefInterfaceBase & i,
				   const InterfacedBase & o,
				   cIBPtr r, int j, const char * s) {
  theMessage << "Could not " << s
	     << " the object \"" << (r? r->name().c_str(): "<NULL>")
	     << " at position " << j << " in the reference vector \""
	     << i.name() << "\" for the object \"" << o.name()
	     << "\" because the " << s
	     << " function threw an unknown exception.";
  severity(setuperror);
}

RefVExSetUnfound::RefVExSetUnfound(const InterfaceBase & i,
				   const InterfacedBase & o,
				   string n) {
  theMessage << "Could not set the object named \""
	     <<  n
	     << " in the reference vector \""
	     << i.name() << "\" of \"" << o.name()
	     << "\"because the object was not found.";
  severity(setuperror);
}

RefVExIndex::RefVExIndex(const InterfaceBase & i,
			 const InterfacedBase & o, int j) {
  theMessage << "Could not access element " << j
	     << " of the reference vector \"" << i.name()
	     << "\" for the object \"" << o.name()
	     << "\" because the index was outside of the allowed range.";
  severity(setuperror);
}

RefVExFixed::RefVExFixed(const InterfaceBase & i, const InterfacedBase & o) {
  theMessage << "Cannot insert or delete in the reference vector \""
	     << i.name() << "\" for the object \"" << o.name()
	     << "\" since the vector is of fixed size.";
  severity(setuperror);
}

RefVExDelUnknown::RefVExDelUnknown(const InterfaceBase & i,
				   const InterfacedBase & o, int j) {
  theMessage << "Could not delete the value at position " << j
	     << " from the reference vector \"" << i.name()
	     << "\" for the object \"" << o.name()
	     << "\" because the delete function threw an unknown exception.";
  severity(setuperror);
}

RefVExGetUnknown::RefVExGetUnknown(const InterfaceBase & i,
				   const InterfacedBase & o) {
  theMessage << "Could not get the reference vector \"" << i.name()
	     << "\" for the object \"" << o.name()
	     << "\" because the get function threw an  unknown exception.";
  severity(setuperror);
}

RefVExNoSet::RefVExNoSet(const InterfaceBase & i,  const InterfacedBase & o) {
  theMessage << "Could not set an object in the reference vector \"" << i.name()
	     << "\" for the object \"" << o.name() 
	     << "\" because no set function has been specified.";
  severity(setuperror);
}


RefVExNoIns::RefVExNoIns(const InterfaceBase & i,  const InterfacedBase & o) {
  theMessage << "Could not insert an object in the reference vector \""
	     << i.name() << "\" for the object \"" << o.name() 
	     << "\" because no insert function has been specified.";
  severity(setuperror);
}


RefVExNoDel::RefVExNoDel(const InterfaceBase & i,  const InterfacedBase & o) {
  theMessage << "Could not erase an object in the reference vector \""
	     << i.name() << "\" for the object \"" << o.name() 
	     << "\" because no erase function has been specified.";
  severity(setuperror);
}

