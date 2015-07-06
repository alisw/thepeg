// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MaxPTScale class.
//

#include "MaxPTScale.h"
#include "DipoleState.h"
#include "Junction.h"
#include "EmitterBase.h"
#include "QCDDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/MaxCmp.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

MaxPTScale::MaxPTScale(): strategy(invariantPT) {}

MaxPTScale::~MaxPTScale() {}

IBPtr MaxPTScale::clone() const {
  return new_ptr(*this);
}

IBPtr MaxPTScale::fullclone() const {
  return new_ptr(*this);
}

Energy MaxPTScale::scale(const DipoleState & state) const {
  MaxCmp<Energy2> maxscale(ZERO);
  MaxCmp<Energy> maxpt(ZERO);
  Energy sumpt = ZERO;

  for ( set<tDBPtr>::const_iterator it = state.activeDipoles().begin();
	it != state.activeDipoles().end(); ++it ) {
    tcQCDPtr d = dynamic_ptr_cast<tcQCDPtr>(*it);
    if ( !d ) continue;
    Energy pt = d->iPart()->momentum().perp();
    sumpt += pt;
    maxpt(pt);
    if ( d->next() ) maxscale(EmitterBase::invPT2(d->iPart(), d->oPart(),
						    d->next()->oPart()));
    else {
      pt = d->oPart()->momentum().perp();
      sumpt += pt;
      maxpt(pt);
    }
  }

  switch ( strategy ) {
  case invariantPT:
    sumpt =  sqrt(maxscale.value());
  case actualPT:
    sumpt = maxpt.value();
  case scalarSum:;
  }

  return sumpt > ZERO? sumpt: state.totalMomentum().m();
  
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MaxPTScale::persistentOutput(PersistentOStream & os) const {
  os << oenum(strategy);
}

void MaxPTScale::persistentInput(PersistentIStream & is, int) {
  is >> ienum(strategy);
}

// The following static variable is needed for the type description
// system in ThePEG.
DescribeClass<MaxPTScale,ScaleSetter>
  describeAriadne5MaxPTScale("Ariadne5::MaxPTScale", "MaxPTScale.so");

void MaxPTScale::Init() {

  static ClassDocumentation<MaxPTScale> documentation
    ("The only task of the MaxPTScale class is to calculate a starti scale "
     "for the dipole shower from a given DipoleState. This class gives the "
     "maximum invariant transverse momentum of all gluons in the sub "
     "process as the starting scale. Alternatively It may give "
     "the maximum actual transverse momentum of all partons or the scalar "
     "sum of all transverse momenta (or the total collision energy if no "
     "gluons are found).");

  
  static Switch<MaxPTScale,PTStrategy> interfaceStrategy
    ("Strategy",
     "The strategy used to determine the maximum transverse momentum for the "
     "starting scale.",
     &MaxPTScale::strategy, invariantPT, true, false);
  static SwitchOption interfaceStrategyInvariantPT
    (interfaceStrategy,
     "InvariantPT",
     "Use the invariant transverse momentum of any gluon.",
     invariantPT);
  static SwitchOption interfaceStrategyActualPT
    (interfaceStrategy,
     "ActualPT",
     "Use the actual transverse momentum of any parton.",
     actualPT);
  static SwitchOption interfaceStrategyScalarSum
    (interfaceStrategy,
     "ScalarSum",
     "Use the scalar sum of all transverse momenta.",
     scalarSum);


}

