// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RecoilGluonEmitter class.
//

#include "RemnantGluonEmission.h"
#include "RecoilGluonEmitter.h"
#include "RemnantModel.h"
#include "Ariadne/Cascade/QCDDipole.h"
#include "Ariadne/Cascade/DipoleState.h"
#include "Ariadne/Cascade/AriadneHandler.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"


#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

RecoilGluonEmitter::RecoilGluonEmitter() {}

RecoilGluonEmitter::~RecoilGluonEmitter() {}

IBPtr RecoilGluonEmitter::clone() const {
  return new_ptr(*this);
}

IBPtr RecoilGluonEmitter::fullclone() const {
  return new_ptr(*this);
}


bool RecoilGluonEmitter::canHandle(const DipoleBase & e) const {
  if ( const QCDDipole * d = dynamic_cast<const QCDDipole *>(&e) ) {
    if ( tRemParPtr r = dynamic_ptr_cast<tRemParPtr>(d->iPart()) )
      if ( !r->hard() ) return true;
    if ( tRemParPtr r = dynamic_ptr_cast<tRemParPtr>(d->oPart()) )
      if ( !r->hard() ) return true;
  }
  return false;
}

bool RecoilGluonEmitter::
overrides(const EmitterBase & em, DipoleBase &) const {
  return false;
}

EmPtr RecoilGluonEmitter::
generate(const DipoleBase & dipole, Energy rhomin, Energy rhomax) const {
  const QCDDipole & d = dynamic_cast<const QCDDipole &>(dipole);

  EmSel sel(rhomin);

  if ( tRemParPtr r = dynamic_ptr_cast<tRemParPtr>(d.iPart()) )
    if ( !r->hard() ) 
      sel = r->model().generateRecoilGluon(*this, d, r, rhomin, rhomax);
  if ( tRemParPtr r = dynamic_ptr_cast<tRemParPtr>(d.oPart()) )
    if ( !r->hard() )
      sel = r->model().generateRecoilGluon(*this, d, r, rhomin, rhomax);

  return sel;
}


bool RecoilGluonEmitter::
perform(const Emission & emission) const {
  const RemnantGluonEmission & e =
    dynamic_cast<const RemnantGluonEmission &>(emission);
  return Current<AriadneHandler>()->remnantModel().performRecoilGluon(e);
}

void RecoilGluonEmitter::revert(const Emission & emission) const {
  const RemnantGluonEmission & e =
    dynamic_cast<const RemnantGluonEmission &>(emission);
  Current<AriadneHandler>()->remnantModel().revertRecoilGluon(e);
  return;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void RecoilGluonEmitter::persistentOutput(PersistentOStream &) const {}

void RecoilGluonEmitter::persistentInput(PersistentIStream &, int) {}

DescribeClass<RecoilGluonEmitter,FSGluonEmitter>
describeAriadne5RecoilGluonEmitter("Ariadne5::RecoilGluonEmitter",
				    "libAriadne5.so");

void RecoilGluonEmitter::Init() {

  static ClassDocumentation<RecoilGluonEmitter> documentation
    ("The RecoilGluonEmitter class implements the emission of recloil gluons "
     "from a remnant dipole.");

}

