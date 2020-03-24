// -*- C++ -*-
//
// Exception.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Exception class.
//

#include "Exception.h"
#include <iostream>
#include <cstdlib>
#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Utilities/Debug.h"

void breakThePEG() {
  return;
}

extern "C" {
  void debugThePEG(const ThePEG::Interfaced * i) {
    i->debug();
  }

  void debugEvent() {
    using namespace ThePEG;
    if ( !CurrentGenerator::isVoid() &&
	 CurrentGenerator::current().currentEvent() )
      cerr << *CurrentGenerator::current().currentEvent();
  }

  long debugEventNumber() {
    using namespace ThePEG;
    if ( !CurrentGenerator::isVoid() )
      return CurrentGenerator::current().currentEventNumber();
    return 0;
  }

  void debugDump() {
    using namespace ThePEG;
    if ( !CurrentGenerator::isVoid() ) CurrentGenerator::current().dump();
  }

  void debugParticle(const ThePEG::Particle * p) {
    using namespace ThePEG;
    cerr << *p;
  }
  void debugParticles(int n, const ThePEG::Particle ** p) {
    using namespace ThePEG;
    LorentzMomentum sum;
    for ( int i = 0; i < n; i++ ) {
      cerr << **p;
      sum += (**p).momentum();
      ++p;
    }
    cerr << ounit(sum,GeV) << "GeV \t" << ounit(sum.m(),GeV) << " GeV\n";
  }

}

namespace ThePEG {

Veto::Veto() {
  if ( ThePEG_DEBUG_LEVEL ) breakThePEG();
}

Exception::Exception(const string & newMessage, Severity newSeverity)
  : theMessage(newMessage), handled(false), theSeverity(newSeverity) {
  breakThePEG();
  if ( noabort && ( theSeverity == abortnow || theSeverity == maybeabort ) )
       theSeverity = runerror;
  if ( theSeverity == abortnow ) {
    writeMessage();
    abort();
  }
}

Exception::~Exception() noexcept {
  if ( !handled ) {
    writeMessage();
    if ( theSeverity == maybeabort ) abort();
  }
}

void Exception::severity(Severity newSeverity) {
  theSeverity = newSeverity;
  if ( noabort && ( theSeverity == abortnow || theSeverity == maybeabort ) )
       theSeverity = runerror;
  if ( theSeverity == abortnow ) {
    writeMessage(cerr);
    abort();
  }
}

void Exception::writeMessage(ostream & os) const {
  switch ( severity() ) {
  case unknown:
    os << "unknown error type: ";
    break;
  case info:
    os << "Informational exception: ";
    break;
  case warning:
    os << "Warning: ";
    break;
  case setuperror:
  case eventerror:
  case runerror:
  case maybeabort:
  case abortnow:
    os << "Error: ";
    break;
  }
  os << message() << endl;
  switch ( severity() ) {
  case eventerror:
    os << "The generated event will be discarded." << endl;
    break;
  case runerror:
    os << "This run will be aborted." << endl;
    break;
  case maybeabort:
  case abortnow:
    os << "The program will now abort and dump core." << endl;
    break;
  case unknown:
  case info:
  case warning:
  case setuperror:
    break;
  }
}

ostream * Exception::errstream = &cerr;

bool Exception::noabort = false;

}
