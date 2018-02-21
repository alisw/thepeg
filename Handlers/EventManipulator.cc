// -*- C++ -*-
//
// EventManipulator.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EventManipulator class.
//

#include "EventManipulator.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;

AbstractNoPIOClassDescription<EventManipulator>
EventManipulator::initEventManipulator;

void EventManipulator::Init() {

  static ClassDocumentation<EventManipulator> documentation
    ("There is no documentation for the ThePEG::EventManipulator class");

}

