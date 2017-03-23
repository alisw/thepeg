// -*- C++ -*-
//
// InterfaceBase.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the InterfaceBase class.
//

#include "InterfaceBase.h"
#include "InterfacedBase.h"
#include "ThePEG/Repository/BaseRepository.h"

namespace ThePEG {

InterfaceBase::InterfaceBase(string newName,
			     string newDescription,
			     string newClassName,
			     const type_info & newTypeInfo, bool depSafe,
			     bool readonly)
  : Named(newName), theDescription(newDescription),
    theClassName(newClassName), theRank(-1.0), hasDefault(true),
    isDependencySafe(depSafe), isReadOnly(readonly) {
  BaseRepository::Register(*this, newTypeInfo);
}

string InterfaceBase::tag(int pos) const {
  if ( pos == -1 ) return name();
  ostringstream os;
  os << name() << "[" << pos << "]";
  return os.str();
}


bool InterfaceBase::notDefault(InterfacedBase & ib) const {
  return exec(ib, "notdef", "") != "" ;
}

map<string,string> & InterfaceBase::objectDefaults(InterfacedBase & ib) const {
  return ib.objectDefaults;
}

string InterfaceBase::fullDescription(const InterfacedBase &) const {
  return type() + '\n' + name() + '\n' + description() +
    ( readOnly()? "\n-*-readonly-*-\n": "\n-*-mutable-*-\n" );
}

void InterfaceBase::doxygenDescription(ostream & os) const {
  os << "\n<hr><b>Name: <a name=\"" << name() << "\"><code>"
     << name() << "</code></a></b><br>\n"
     << "<b>Type:</b> " << doxygenType();
  if ( readOnly() ) os << " (read-only)";
  os << " <br>\n"
     << "\\par Description:\n"
     << description() << "<br>\n";
}

bool InterfaceBase::NoReadOnly = false;

InterExClass::InterExClass(const InterfaceBase & i, const InterfacedBase & o) {
  theMessage << "Could not access the interface \"" << i.name()
	     << "\" of the object \"" << o.name() << "\" because the object "
	     << "is not of the required class (" << i.className() << ").";
  severity(setuperror);
}

InterExSetup::InterExSetup(const InterfaceBase & i, const InterfacedBase & o) {
  theMessage << "Could not access the interface \"" << i.name()
	     << "\" for the object \"" << o.name()
	     << "\" since no get/set member function or variable was found.";
  severity(setuperror);
}

InterExUnknown::InterExUnknown(const InterfaceBase & i,
			       const InterfacedBase & o) {
  theMessage << "Could not perform action on the interface  \""
	     << i.name() << "\" for the object \"" << o.name()
	     << "\" because the requested action was not recognized";
  severity(setuperror);
}

InterExReadOnly::InterExReadOnly(const InterfaceBase & i,
				 const InterfacedBase & o) {
  theMessage << "Could not perform action on the interface  \""
	     << i.name() << "\" for the object \"" << o.name()
	     << "\" because this interface is read-only.";
  severity(setuperror);
}

InterExNoNull::InterExNoNull(const InterfaceBase & i,
			     const InterfacedBase & o) {
  theMessage << "Could not set reference \""
	     << i.name() << "\" for the object \"" << o.name()
	     << "\" to <Null> because null pointers are explicitly "
	     << "disallowed.";
  severity(setuperror);
}


RefInterfaceBase::
RefInterfaceBase(string newName, string newDescription, string newClassName,
		 const type_info & newTypeInfo, string newRefClassName,
		 const type_info & newRefTypeInfo, bool depSafe,
		 bool readonly, bool norebind, bool nullable, bool defnull)
  : InterfaceBase(newName, newDescription, newClassName, newTypeInfo, depSafe,
		  readonly), theRefClassName(newRefClassName),
  theRefTypeInfo(newRefTypeInfo), dontRebind(norebind),
  isNullable(nullable), theDefaultIfNull(defnull) {
  hasDefault = false;
}


}
