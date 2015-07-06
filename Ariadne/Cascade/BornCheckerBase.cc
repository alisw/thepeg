// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BornCheckerBase class.
//

#include "BornCheckerBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Ariadne5;

BornCheckerBase::BornCheckerBase() {}

BornCheckerBase::~BornCheckerBase() {}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

// The following static variable is needed for the type description
// system in ThePEG.
DescribeAbstractNoPIOClass<BornCheckerBase,Interfaced>
  describeAriadne5BornCheckerBase("Ariadne5::BornCheckerBase", "libAriadne5.so");

void BornCheckerBase::Init() {

  static ClassDocumentation<BornCheckerBase> documentation
    ("The BornCheckerBase is the base class for any class implementing "
     "the checking of a DipoleState reclustered in the CKKW-L algorithm, "
     "to make sure it corresponds to a reasonable Born-level stat (the "
     "lowest jet-multiplicity state in a CKKW-L merging). Besided the "
     "main virtual function check() also the scale() function must be "
     "overridden in sub-classes.");

}

