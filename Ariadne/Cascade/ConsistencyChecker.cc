// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ConsistencyChecker class.
//

#include "ConsistencyChecker.h"
#include "DipoleState.h"
#include "Emission.h"
#include "Parton.h"
#include "QCDDipole.h"
#include "AriadneHandler.h"
#include "Models/RemnantModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Utilities/DebugItem.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

ConsistencyChecker::~ConsistencyChecker() {}

IBPtr ConsistencyChecker::clone() const {
  return new_ptr(*this);
}

IBPtr ConsistencyChecker::fullclone() const {
  return new_ptr(*this);
}

bool ConsistencyChecker::
isConsistent(const DipoleState & state, tcEmPtr) const {
  if ( checkGluonPT && !gluonPT(state) ) return false;
  if ( checkStrings && !strings(state) ) return false;
  if ( !remnantPhaseSpace(state) ) return false;
  return true;
}

bool ConsistencyChecker::gluonPT(const DipoleState & state) const {
  Energy2 pt2cut = sqr(Current<AriadneHandler>()->pTCut());

  for ( set<tDBPtr>::const_iterator it = state.activeDipoles().begin();
	it != state.activeDipoles().end(); ++it )
    if ( tcQCDPtr d = dynamic_ptr_cast<tcQCDPtr>(*it) )
      if ( d->next() &&
	   !gluonPT(d->iPart(), d->oPart(), d->next()->oPart(), pt2cut) )
	return false;
  return true;
}

bool ConsistencyChecker::remnantPhaseSpace(const DipoleState & state) const {
  static DebugItem debugrem("Ariadne5::ConsistencyChecker::Remnant");
  const RemnantModel & remmod = Current<AriadneHandler>()->remnantModel();
  Lorentz5Momentum ph = state.hadronicMomentum();
  for ( int i = 0, N = state.remnants().first.size(); i < N; ++i ) {
    tcRemParPtr rem = state.remnants().first[i];
    if ( rem->hard() ) continue;
    if ( rem->recoilWeight() <= 0.0 ) continue;
    Energy effMr = remmod.effectiveMass(rem);
    if ( effMr <= ZERO ) continue;
    if ( (ph + rem->momentum()).m2() <
	 sqr(ph.mass() + effMr + rem->momentum().perp()) ) {
      if ( debugrem ) {
	cerr << "ConsistencyChecker '" << name()
	   << "' found remnant with too small momentum:" << endl;
	rem->debug();
	cerr << endl;
      }
      return false;
    }
  }
  for ( int i = 0, N = state.remnants().second.size(); i < N; ++i ) {
    tcRemParPtr rem = state.remnants().second[i];
    if ( rem->hard() ) continue;
    Energy effMr = remmod.effectiveMass(rem);
    if ( effMr <= ZERO ) continue;
    if ( (ph + rem->momentum()).m2() <
	 sqr(ph.mass() + effMr + rem->momentum().perp()) ) {
      if ( debugrem ) {
	cerr << "ConsistencyChecker '" << name()
	   << "' found remnant with too small momentum:" << endl;
	rem->debug();
	cerr << endl;
      }
      return false;
    }
  }
  return true;
}

bool ConsistencyChecker::
gluonPT(tcParPtr p1, tcParPtr p2, tcParPtr p3, Energy2 pt2cut) const {
  static DebugItem debugpt("Ariadne5::ConsistencyChecker::PT");
  if ( p1->failsafe() || p2->failsafe() || p3->failsafe() ) return true;
  if ( p1->touched() || p2->touched() || p3->touched() ) {
    if ( EmitterBase::invPT2(p1, p2, p3) > pt2cut ) return true;
    if ( debugpt ) {
      cerr << "ConsistencyChecker '" << name()
	   << "' found pt below cut:" << endl;
      p1->debug();
      cerr << endl;
      p2->debug();
      cerr << endl;
      p3->debug();
      cerr << endl;
    }
    return false;
  }
  return true;
}

bool ConsistencyChecker::strings(const DipoleState &) const {
  // *** ATTENTION *** Implement this.
  return false;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ConsistencyChecker::persistentOutput(PersistentOStream & os) const {
  os << vetoInitial << checkGluonPT << checkStrings;
}

void ConsistencyChecker::persistentInput(PersistentIStream & is, int) {
  is >> vetoInitial >> checkGluonPT >> checkStrings;
}


// The following static variable is needed for the type description
// system in ThePEG.
DescribeClass<ConsistencyChecker,HandlerBase>
  describeAriadne5ConsistencyChecker("Ariadne5::ConsistencyChecker",
				     "libAriadne5.so");

void ConsistencyChecker::Init() {

  static ClassDocumentation<ConsistencyChecker> documentation
    ("The ConsistencyChecker class is used by the AriadneHandler to "
     "check a DipoleState for consistency. The check is done on the "
     "initial DipoleState, and then after each Emission, to veto any "
     "emission which gives an inconsistent state. If the initial "
     "DipoleState is inconsistent, the whole event may be vetoed, "
     "otherwise further checking is suspended. This base class "
     "optionally checks that all gluons have an invariant transverse "
     "momentum above the cutoff and that all stings are large enough "
     "to be able to produce a minimum set of hadrons. Sub-classes may "
     "introduce further checks.");


  static Switch<ConsistencyChecker,bool> interfaceInitialVeto
    ("InitialVeto",
     "If the initial DipoleState is found inconsistent, veto the event. "
     "If false, the event is accepted and further checks are suspended.",
     &ConsistencyChecker::vetoInitial, true, true, false);
  static SwitchOption interfaceInitialVetoVeto
    (interfaceInitialVeto,
     "Veto",
     "The event is vetoed if found inconsistent.",
     true);
  static SwitchOption interfaceInitialVetoSuspend
    (interfaceInitialVeto,
     "Suspend",
     "If the initial state is inconsistent, suspend further checks.",
     false);

  static Switch<ConsistencyChecker,bool> interfaceCheckGluonPT
    ("CheckGluonPT",
     "Checks that all gluons have an invariant transverse momentum "
     "above the cutoff.",
     &ConsistencyChecker::checkGluonPT, false, true, false);
  static SwitchOption interfaceCheckGluonPTOn
    (interfaceCheckGluonPT,
     "On",
     "The check is performed.",
     true);
  static SwitchOption interfaceCheckGluonPTOff
    (interfaceCheckGluonPT,
     "Off",
     "The check is not performed.",
     false);

  static Switch<ConsistencyChecker,bool> interfaceCheckStrings
    ("CheckStrings",
     "Check that all strings have an invariant mass big enough to "
     "produce a minimum set of hadrons.",
     &ConsistencyChecker::checkStrings, false, true, false);
  static SwitchOption interfaceCheckStringsOn
    (interfaceCheckStrings,
     "On",
     "Strings are checked.",
     true);
  static SwitchOption interfaceCheckStringsOff
    (interfaceCheckStrings,
     "Off",
     "Strings are not checked.",
     false);

  static Switch<ConsistencyChecker,bool> interfaceTryFixing
    ("TryFixing",
     "If the initial state is inconsistent but not vetoed, a subclass may implement a function to possibly fix the state to become consistent.",
     &ConsistencyChecker::tryFixing, false, true, false);
  static SwitchOption interfaceTryFixingOn
    (interfaceTryFixing,
     "On",
     "Will try to fix a broken intial dipole state. ",
     true);
  static SwitchOption interfaceTryFixingOff
    (interfaceTryFixing,
     "Off",
     "No fixing will be attempted.",
     false);
}

