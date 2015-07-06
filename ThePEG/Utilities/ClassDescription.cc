// -*- C++ -*-
//
// ClassDescription.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ClassDescription class.
//

#include "ClassDescription.h"

using namespace ThePEG;

ClassDescriptionBase::~ClassDescriptionBase() {}

bool ClassDescriptionBase::isA(const ClassDescriptionBase & base) const {
  if ( base.info() == info() ) return true;
  for ( DescriptionVector::const_iterator i = descriptions().begin();
	i != descriptions().end(); ++i ) 
    if ( (**i).isA(base) ) return true;
  return false;
}

