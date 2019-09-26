// -*- C++ -*-
//
// ColourBase.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ColourBase class.
//

#include "ColourBase.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

EIPtr ColourBase::clone() const {
  return new_ptr(*this);
}

bool ColourBase::hasColourLine(tcColinePtr line, bool anti) const {
  return ( anti? ( antiColourLine() == line ): ( colourLine() == line ) );
}

vector<tcColinePtr> ColourBase::antiColourLines() const {
  return antiColourLine()? vector<tcColinePtr>(1, antiColourLine()):
    vector<tcColinePtr>();
}

vector<tcColinePtr> ColourBase::colourLines() const {
  return colourLine()? vector<tcColinePtr>(1, colourLine()):
    vector<tcColinePtr>();
}

void ColourBase::rebind(const EventTranslationMap & trans) {
  theAntiColourLine = trans.translate(theAntiColourLine);
  theColourLine = trans.translate(theColourLine);
}

void ColourBase::persistentOutput(PersistentOStream & os) const {
  os << theAntiColourLine << theColourLine;
}

void ColourBase::persistentInput(PersistentIStream & is, int) {
  is >> theAntiColourLine >> theColourLine;
}

ClassDescription<ColourBase> ColourBase::initColourBase;

void ColourBase::Init() {}

