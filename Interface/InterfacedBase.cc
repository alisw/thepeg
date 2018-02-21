// -*- C++ -*-
//
// InterfacedBase.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the InterfacedBase class.
//

#include "InterfacedBase.h"
#include "ThePEG/Repository/BaseRepository.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/Utilities/DescriptionList.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"

using namespace ThePEG;

InterfacedBase::~InterfacedBase() {}

void InterfacedBase::readSetup(istream &) {}

bool InterfacedBase::preInitialize() const {
  return false;
}

void InterfacedBase::persistentOutput(PersistentOStream & os) const {
  os << fullName() << isLocked << isTouched << oenum(initState) << theComment
     << objectDefaults;
}

void InterfacedBase::persistentInput(PersistentIStream & is, int) {
  string n;
  is >> n >> isLocked >> isTouched >> ienum(initState) >> theComment
     >> objectDefaults;
  name(n);
}

string InterfacedBase::addComment(string c) {
  if ( theComment.length() ) theComment += "\n";
  theComment += StringUtils::stripws(c);
  return "";
}

AbstractClassDescription<InterfacedBase> InterfacedBase::initInterfacedBase;

void InterfacedBase::debugme() const {
  cerr << name() << " ["
       << DescriptionList::find(typeid(*this))->name() << "] ";
  PersistentBase::debugme();
}

void InterfacedBase::Init() {

  static Parameter<InterfacedBase,string> interfaceComment
    ("Comment",
     "A comment assigned to this object.",
     &InterfacedBase::theComment, "", true, false);
  interfaceComment.setHasDefault(false);

  static Command<InterfacedBase> interfaceAddComment
    ("AddComment",
     "Add a comment to this object. Will be concatenated with the exixting "
     "comment.",
     &InterfacedBase::addComment, true);

}

void InterfacedBase::UpdateChecker::check(tIBPtr ip, bool & touch) {
  if ( !ip ) return;
  ip->update();
  if ( ip->touched() ) touch = true;
}

