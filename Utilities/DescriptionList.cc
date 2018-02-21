// -*- C++ -*-
//
// DescriptionList.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DescriptionList class.
//

#include "DescriptionList.h"
#include "ClassDescription.h"

using namespace ThePEG;

void DescriptionList::hookup() {
  for ( DescriptionMap::iterator it = descriptionMap().begin();
	it != descriptionMap().end(); ++it ) it->second->setup();
}

void DescriptionList::Register(ClassDescriptionBase & pd) {
  if ( find(pd.info()) ) return;
  insert(pd);
  hookup();
}

DescriptionList::DescriptionMap & DescriptionList::descriptionMap() {
  static DescriptionMap theDescriptionMap;
  return theDescriptionMap;
}

DescriptionList::StringMap & DescriptionList::stringMap() {
  static StringMap theStringMap;
  return theStringMap;
}

string DescriptionList::className(const type_info & ti) {
  const ClassDescriptionBase * d = find(ti);
  return d? d->name(): string();
}

int DescriptionList::version(const type_info & ti) {
  const ClassDescriptionBase * d = find(ti);
  return d? d->version(): 0;
}

string DescriptionList::library(const type_info & ti) {
  const ClassDescriptionBase * d = find(ti);
  return d? d->library(): string();
}

void DescriptionList::insert(ClassDescriptionBase & pb) {
#ifndef THEPEG_DYNAMIC_TYPE_INFO_BUG
  descriptionMap()[&(pb.info())] = &pb;
#else
  descriptionMap()[pb.info().name()] = &pb;
#endif
  stringMap()[pb.name()] = &pb;
  // The following is for backward compatibility and will eventually
  // be removed.
  string name = pb.name();
  if ( name.find('/') == string::npos ) {
    if ( name.substr(0, 2) == "::" ) name = name.substr(2);
    if ( name.find("::") == string::npos ) return;
    while ( name.find("::") != string::npos )
      name.replace(name.find("::"), 2, "/");
    name = "/" + name;
    stringMap()[name] = &pb;
  } else {
    if ( name[0] == '/' ) name = name.substr(1);
    while ( name.find('/') != string::npos )
      name.replace(name.find('/'), 1, "::");
    stringMap()[name] = &pb;
  }
}

void DescriptionList::printHierarchies(ostream & os) {
  for ( DescriptionMap::iterator it = descriptionMap().begin();
	it != descriptionMap().end(); ++it ) {
    os << "Class Name '" << it->second->name() << "'\n ("
#ifndef THEPEG_DYNAMIC_TYPE_INFO_BUG
       << it->first->name()
#else
       << it->first
#endif
       << "," <<  it->second << ") version " << it->second->version()
       << endl << "  Base classes:" << endl;
    for ( unsigned int i = 0; i < it->second->descriptions().size(); ++i )
      os << "   " << i << " '" << it->second->descriptions()[i]->name()
	 << "' (" << it->second->descriptions()[i] << ")" << endl;
  }
}
