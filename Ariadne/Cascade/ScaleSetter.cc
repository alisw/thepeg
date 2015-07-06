// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScaleSetter class.
//

#include "ScaleSetter.h"
#include "DipoleState.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Ariadne5;

ScaleSetter::ScaleSetter() {}

ScaleSetter::~ScaleSetter() {}

IBPtr ScaleSetter::clone() const {
  return new_ptr(*this);
}

IBPtr ScaleSetter::fullclone() const {
  return new_ptr(*this);
}

Energy ScaleSetter::scale(const DipoleState & state) const {
  return sqrt(state.totalMomentum().m2());
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

// The following static variable is needed for the type description
// system in ThePEG.
DescribeNoPIOClass<ScaleSetter,Interfaced>
  describeAriadne5ScaleSetter("Ariadne5::ScaleSetter", "libAriadne5.so");

void ScaleSetter::Init() {

  static ClassDocumentation<ScaleSetter> documentation
    ("The only task of the ScaleSetter class is to calculate a starting "
     "scale for the dipole shower from a given DipoleState. This base class "
     "only gives the total collision energy, but other strategies can be "
     "implemented in sub-classes overriding the virtual scale( function.");

}

