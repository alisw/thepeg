// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MinPTScale class.
//

#include "MinPTScale.h"
#include "DipoleState.h"
#include "Junction.h"
#include "EmitterBase.h"
#include "QCDDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/MaxCmp.h"

using namespace Ariadne5;

MinPTScale::MinPTScale() {}

MinPTScale::~MinPTScale() {}

IBPtr MinPTScale::clone() const {
  return new_ptr(*this);
}

IBPtr MinPTScale::fullclone() const {
  return new_ptr(*this);
}

Energy MinPTScale::scale(const DipoleState & state) const {
  MinCmp<Energy2> minscale(state.totalMomentum().m2());

  for ( set<tDBPtr>::const_iterator it = state.activeDipoles().begin();
	it != state.activeDipoles().end(); ++it )
    if ( tcQCDPtr d = dynamic_ptr_cast<tcQCDPtr>(*it) )
      if ( d->next() ) minscale(EmitterBase::invPT2(d->iPart(), d->oPart(),
						    d->next()->oPart()));

  return sqrt(minscale.value());
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

// The following static variable is needed for the type description
// system in ThePEG.
DescribeNoPIOClass<MinPTScale,ScaleSetter>
  describeAriadne5MinPTScale("Ariadne5::MinPTScale", "MinPTScale.so");

void MinPTScale::Init() {

  static ClassDocumentation<MinPTScale> documentation
    ("The only task of the MinPTScale class is to calculate a starti scale "
     "for the dipole shower from a given DipoleState. This class gives the "
     "minimum invariant transverse momentum of all gluons in the sub "
     "process as the starting scale (or the total collision energy if no "
     "gluons are found).");

}

