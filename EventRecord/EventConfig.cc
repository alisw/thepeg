// -*- C++ -*-
//
// EventConfig.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EventConfig class.
//

#include "EventConfig.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EventConfig.tcc"
#endif


using namespace ThePEG;

tcEventBasePtr EventConfig::currentGenerator;

void EventConfig::putHandler(PersistentOStream & os, tcEventBasePtr h) {
  if ( !currentGenerator ) os << "" << h;
  else {
    tcIPtr handler = dynamic_ptr_cast<tcIPtr>(h);
    if ( handler ) os << handler->fullName();
    else os << "";
    os << tcEventBasePtr();
  }
}

void EventConfig::getHandler(PersistentIStream & is, tcEventBasePtr & h) {
  string pxh;
  is >> pxh >> h;
  if ( currentGenerator ) {
    tcEGPtr eg = dynamic_ptr_cast<tcEGPtr>(currentGenerator);
    h = eg->getObject<EventRecordBase>(pxh);
  }
}

void EventConfig::putParticleData(PersistentOStream & os, tcEventPDPtr pd) {
  if ( !currentGenerator ) os << 0 << pd;
  else os << pd->id() << tcEventPDPtr();
}

void EventConfig::getParticleData(PersistentIStream & is, cEventPDPtr & pd) {
  long pid;
  is >> pid >> pd;
  if ( !pd && currentGenerator ) {
    tcEGPtr eg = dynamic_ptr_cast<tcEGPtr>(currentGenerator);
    pd = eg->getParticleData(pid);
  }
}

string EventConfig::nameHandler(tcEventBasePtr h) {
  tcIPtr handler = dynamic_ptr_cast<tcIPtr>(h);
  if ( handler ) return handler->name();
  return "";
}

