// -*- C++ -*-
//
// HadronizationHandler.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HadronizationHandler class.
//

#include "HadronizationHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;

AbstractNoPIOClassDescription<HadronizationHandler>
HadronizationHandler::initHadronizationHandler;

void HadronizationHandler::Init() {

  static ClassDocumentation<HadronizationHandler> documentation
    ("This is the base class to be used by all models for hadronization.");

}

