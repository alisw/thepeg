// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EmissionGenerator class.
//

#include "EmissionGenerator.h"
#include "AriadneHandler.h"
#include "EmitterBase.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Utilities/MaxCmp.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"


using namespace Ariadne5;

bool EmissionGenerator::init() const {
  emitters.clear();
  emission = EmPtr();
  exhausted = false;
  const vector<EmitterPtr> & potential = Current<AriadneHandler>()->emitters();
  for ( int i = potential.size() - 1; i >= 0; --i ) {
    bool ok = true;
    if ( !potential[i]->canHandle(*dipole) ) continue;
    list<tEmitterPtr>::iterator it = emitters.begin();
    while ( it != emitters.end() ) {
      if ( (**it).overrides(*potential[i], *dipole) ) {
	ok = false;
	break;
      }
      if ( potential[i]->overrides(**it, *dipole) )
	it = emitters.erase(it);
      else
	++it;
    }
    if ( ok ) emitters.push_back(potential[i]);
  }
  return !emitters.empty();
}

bool EmissionGenerator::reinit() const {
  if ( dipole->touched() ) return init();
  for ( list<tEmitterPtr>::iterator it = emitters.begin();
	it != emitters.end(); ++it )
    if ( (**it).touched(*dipole) ) return init();
  return false;
}

Energy EmissionGenerator::generate(Energy rhomin, Energy rhomax) const {
  if ( emission || exhausted ) return rho();
  //  rhomin = max(rhomin, dipole->rhoCut());
  MinCmp<Energy> mincut;
  for ( list<tEmitterPtr>::iterator it = emitters.begin();
	it != emitters.end(); ++it ) {
    Energy rhocut = (**it).rhoCut();
    mincut(rhocut);
    EmPtr em = (**it).generate(*dipole, max(rhomin, rhocut), rhomax);
    if ( em && em->rho > max(rhomin, rhocut) ) {
      emission = em;
      rhomin = em->rho;
      em->dipole = dipole;
    }
  }
  if ( !emission && rhomin <= mincut.value() ) exhausted = true;
  return rho();
}

void EmissionGenerator::persistentOutput(PersistentOStream & os) const {
  os << dipole << emitters << emission << exhausted;
}

void EmissionGenerator::persistentInput(PersistentIStream & is) {
  is >> dipole >> emitters >> emission >> exhausted;
}

