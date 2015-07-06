// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RndmEngine class.
//

#include "RndmEngine.h"
#include "ThePEG/Repository/UseRandom.h"

using namespace TheP8I;

RndmEngine::RndmEngine() {}

RndmEngine::~RndmEngine() {}

double RndmEngine::flat() {
  return UseRandom::rnd();
}

