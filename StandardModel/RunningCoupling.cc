// -*- C++ -*-
//
// RunningCoupling.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RunningCoupling class.
//

#include "RunningCoupling.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

AbstractClassDescription<RunningCoupling> RunningCoupling::initRunningCoupling;

void RunningCoupling::persistentOutput(PersistentOStream & os) const {
  os << theScaleFactor;
}

void RunningCoupling::persistentInput(PersistentIStream & is, int) {
  is >> theScaleFactor;
}

void RunningCoupling::Init() {

  static ClassDocumentation<RunningCoupling> documentation
    ("An abstract base class used to implement running couplings.");

  static Parameter<RunningCoupling,double> interfaceScaleFactor
    ("ScaleFactor",
     "The scale factor used to globally rescale the argument of the running coupling",
     &RunningCoupling::theScaleFactor, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  interfaceScaleFactor.rank(-1);

}

