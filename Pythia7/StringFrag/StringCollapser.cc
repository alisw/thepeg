// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StringCollapser class.
//

#include "StringCollapser.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "StringCollapser.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Pythia7;

StringCollapser::~StringCollapser() {}

void StringCollapser::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void StringCollapser::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<StringCollapser> StringCollapser::initStringCollapser;
// Definition of the static class description member.

void StringCollapser::Init() {

  static ClassDocumentation<StringCollapser> documentation
    ("Used by the LundFragHandler class to collapse strings which are "
     "deemed too small to fragment into one or two hadrons.");

}

