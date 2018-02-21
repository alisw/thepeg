// -*- C++ -*-
//
// Interfaced.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Interfaced class.
//

#include "Interfaced.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Interface/Command.h"

using namespace ThePEG;

Interfaced::~Interfaced() {}

void Interfaced::registerRepository(IBPtr i) {
  Repository::Register(i);
}

void Interfaced::registerRepository(IBPtr i, string newName) {
  Repository::Register(i, newName);
}

void Interfaced::reporeg(IBPtr object, string name) const {
  Repository::CreateDirectory(fullName());
  string full = fullName() + "/" + name;
  IBPtr old = Repository::GetPointer(full);
  if ( old ) {
    if ( Repository::GetObjectsReferringTo(old).empty() )
      Repository::remove(old);
    else
      Repository::rename(old, fullName() + "/old-" + name);
  }
  Repository::Register(object, full);
}

bool Interfaced::defaultInit() {
  return true;
}

void Interfaced::setUsed() const {
  theUseFlag = true;
  if ( generator() ) generator()->use(*this);
}

PPtr Interfaced::getParticle(PID newId) const {
  PPtr p(generator()? generator()->getParticle(newId): PPtr());
  return p;
}

PDPtr Interfaced::getParticleData(PID newId) const {
  PDPtr p(generator()? generator()->getParticleData(newId):
	  Repository::defaultParticle(newId));
  return p;
}

void Interfaced::persistentOutput(PersistentOStream & os) const {
  os << theGenerator << theUseFlag;
}

void Interfaced::persistentInput(PersistentIStream & is , int) {
  is >> theGenerator >> theUseFlag;
}

AbstractClassDescription<Interfaced> Interfaced::initInterfaced;

string Interfaced::doDefaultInit(string) {
  if ( !defaultInit() ) return "Default initialization failed.";
  return "";
}

void Interfaced::Init() {

  static Command<Interfaced> interfaceDefaultInit
    ("DefaultInit",
     "Perform a default initialization of this object. This typically "
     "involves creating sub-objects which are needed. In this case the "
     "objects can be added to the repository in a sub-directory with the "
     "same name as this object.",
     &Interfaced::doDefaultInit, true);

}

