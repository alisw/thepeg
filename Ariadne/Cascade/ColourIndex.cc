// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ColourIndex class.
//

#include "ColourIndex.h"
#include "ThePEG/Repository/UseRandom.h"
#include "AriadneHandler.h"
#include "ThePEG/Utilities/Current.h"

using namespace Ariadne5;

void ColourIndex::generate(const ColourIndex & c1,
			   const ColourIndex & c2,
			   const ColourIndex & c3) {
  do {
    index = UseRandom::rnd(Current<AriadneHandler>()->nCol()) + 1;
  } while ( index == c1.index || index == c2.index || index == c3.index );
}

