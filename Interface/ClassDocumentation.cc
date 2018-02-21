// -*- C++ -*-
//
// ClassDocumentation.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ClassDocumentation class.
//

#include "ClassDocumentation.h"
#include "ThePEG/Repository/BaseRepository.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ClassDocumentation.tcc"
#endif

using namespace ThePEG;


ClassDocumentationBase::
ClassDocumentationBase(string newDocumentation,
		       string newModelDescription,
		       string newModelReferences,
		       const type_info & newTypeInfo)
  : theDocumentation(newDocumentation),
    theModelDescription(newModelDescription),
    theModelReferences(newModelReferences) {
  BaseRepository::Register(*this, newTypeInfo);
}

