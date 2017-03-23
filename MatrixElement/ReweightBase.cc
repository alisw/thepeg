// -*- C++ -*-
//
// ReweightBase.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ReweightBase class.
//

#include "ReweightBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

ReweightBase::~ReweightBase() {}

void ReweightBase::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb;
}

void ReweightBase::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb;
}

void ReweightBase::setXComb(tXCombPtr xc) {
  theLastXComb = xc;
}

AbstractClassDescription<ReweightBase> ReweightBase::initReweightBase;
// Definition of the static class description member.

void ReweightBase::Init() {

  static ClassDocumentation<ReweightBase> documentation
    ("There is no documentation for the ThePEG::ReweightBase class");

}

