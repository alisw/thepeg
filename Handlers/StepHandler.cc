// -*- C++ -*-
//
// StepHandler.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StepHandler class.
//

#include "StepHandler.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;

StepHandler::~StepHandler() {}

void StepHandler::eventHandler(tEHPtr eh) {
  theEventHandler = eh;
  theNewStep = tStepPtr();
  theCurrentStep = eh->currentStep();
}

void StepHandler::createNewStep() {
  useMe();
  theNewStep = eventHandler()->newStep(this);
}

AbstractNoPIOClassDescription<StepHandler> StepHandler::initStepHandler;

void StepHandler::Init() {

  static ClassDocumentation<StepHandler> documentation
    ("There is no documentation for the ThePEG::StepHandler class");

}

