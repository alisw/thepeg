// -*- C++ -*-
//
// HelicityVertex.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2011 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HelicityVertex class.
//
// Author: Peter Richardson
//

#include "SpinInfo.h"
#include "HelicityVertex.h"

using namespace ThePEG;

AbstractNoPIOClassDescription<HelicityVertex>
HelicityVertex::initHelicityVertex;
// Definition of the static class description member.

void HelicityVertex::Init() {}

void HelicityVertex::rebind(const EventTranslationMap & trans) {
  EventInfoBase::rebind(trans);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    _incoming[ix]=trans.translate(_incoming[ix]);
  }
  for(unsigned int ix=0;ix<_outgoing.size();++ix) {
    _outgoing[ix]=trans.translate(_outgoing[ix]);
  }
}
