// -*- C++ -*-
//
// SubProcessGroup.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
// Copyright (C) 2009-2017 Simon Platzer
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SubProcessGroup class.
//
#include "SubProcessGroup.h"
#include "ThePEG/EventRecord/Collision.h"
#include "ThePEG/Config/algorithm.h"
#include "ThePEG/EventRecord/ParticleTraits.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <iostream>

using namespace ThePEG;

SubProcessGroup::
SubProcessGroup(const PPair & newIncoming,
 	   tCollPtr newCollision, tcEventBasePtr newHandler)
  : SubProcess(newIncoming,newCollision,newHandler) {}

SubProcessGroup::~SubProcessGroup() {}

SubProPtr SubProcessGroup::clone() const {
  return ptr_new<Ptr<SubProcessGroup>::ptr>(*this);
}

void SubProcessGroup::rebind(const EventTranslationMap & trans) {

  SubProcess::rebind(trans);

  for ( SubProcessVector::iterator sub = dependent().begin();
	sub != dependent().end(); ++sub )
    (*sub = trans.translate(*sub))->rebind(trans);

}

void SubProcessGroup::transform(const LorentzRotation & r) {

  SubProcess::transform(r);

  for ( SubProcessVector::iterator sub = dependent().begin();
	sub != dependent().end(); ++sub )
    (**sub).transform(r);

}

void SubProcessGroup::printMe(ostream& os) const {
  os << "head sub-process of this group with relative weight "
     << groupWeight() << ":\n";
  SubProcess::printMe(os);
  os << "dependent sub-processes in this group:\n";
  for ( SubProcessVector::const_iterator sub = dependent().begin();
	sub != dependent().end(); ++sub ) {
    os << "performed by " << EventConfig::nameHandler((**sub).handler())
       << " with relative weight " << (**sub).groupWeight() << "\n";
    (**sub).printMe(os);
  }
}

void SubProcessGroup::persistentOutput(PersistentOStream & os) const {
  os << theDependent;

}

void SubProcessGroup::persistentInput(PersistentIStream & is, int) {
  is >> theDependent;
}

ClassDescription<SubProcessGroup> SubProcessGroup::initSubProcessGroup;

void SubProcessGroup::Init() {}

