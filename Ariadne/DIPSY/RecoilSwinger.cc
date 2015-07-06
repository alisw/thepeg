// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RecoilSwinger class.
//

#include "RecoilSwinger.h"
#include "Dipole.h"
#include "DipoleState.h"
#include "DipoleEventHandler.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Utilities/Current.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

RecoilSwinger::RecoilSwinger() {}

RecoilSwinger::~RecoilSwinger() {}

IBPtr RecoilSwinger::clone() const {
  return new_ptr(*this);
}

IBPtr RecoilSwinger::fullclone() const {
  return new_ptr(*this);
}

void RecoilSwinger::
generate(Dipole & dipole, double miny, double maxy, bool force) const {
  Swinger::generate(dipole, miny, maxy, force);
//   generateExchange(dipole, miny, maxy, force);
}

void RecoilSwinger::recombine(Dipole & dip) const {
  DipolePtr d1 = & dip;
  DipolePtr d2 = d1->swingDipole();
//   if ( !kinematicsVeto(d1->partons(), d2->partons()) &&
//        exchangeAmp(d1->partons(), d2->partons())/
//        swingAmp(d1->partons(), d2->partons()) > UseRandom::rnd() ) {
//     cout << "recombining. swingamp: " << swingAmp(d1->partons(), d2->partons())
// 	 << ", exchangeamp: " << exchangeAmp(d1->partons(), d2->partons()) << endl;
//     exchange(*d1);
//   }
  if ( 1.0/3.0 > UseRandom::rnd() ) {
    exchange(*d1);
  }
  else
    Swinger::recombine(*d1);
}

void RecoilSwinger::exchange(Dipole & d1) const {

//   cout << "doing gluon exchange!!" << endl;
//   d1.dipoleState().diagnosis(true);

  Dipole & d2 = *d1.swingDipole();
  DipoleXSec::InteractionRecoil recoil =
    Current<DipoleEventHandler>()->xSecFn().recoil(d1.partons(), d2.partons(), ImpactParameters());
  d1.partons().first->pT()  += recoil.first.first;
  d1.partons().second->pT() += recoil.first.second;
  d2.partons().first->pT()  += recoil.second.first;
  d2.partons().second->pT() += recoil.second.second;
  d1.partons().first->updateYMinus();
  d1.partons().second->updateYMinus();
  d2.partons().first->updateYMinus();
  d2.partons().second->updateYMinus();

  Swinger::recombine(d1);

//   cout << "done gluon exchange for a pt of " << (recoil.first.first.pt() + recoil.first.second.pt() + recoil.second.first.pt() + recoil.second.second.pt())/GeV << endl;
//   d1.dipoleState().diagnosis(true);
}

double RecoilSwinger::
exchangeAmp(const pair<tPartonPtr, tPartonPtr> first,
	    const pair<tPartonPtr, tPartonPtr> second) const {
  return Current<DipoleEventHandler>()->xSecFn().fij(first, second, ImpactParameters())*9.0;
}

bool RecoilSwinger::
kinematicsVeto(const pair<tPartonPtr, tPartonPtr> first,
	       const pair<tPartonPtr, tPartonPtr> second) const {
  return false; //the veto is included in fij
//   return Current<DipoleEventHandler>()->xSecFn().kinematicsVeto(d1, d2, ImpactParameters());
}

void RecoilSwinger::
generateExchange(Dipole & dipole, double miny, double maxy, bool force) const {
  cout << "generating exchange for dip at "
       << max(max(dipole.partons().first->y(), dipole.partons().second->y()), miny) << endl;

  const list<DipolePtr> & candidates = dipole.dipoleState().getDipoles();
  list<DipolePtr>::const_iterator it = candidates.begin();
  while( it != candidates.end() && *(it++) != & dipole );
  for (; it != candidates.end(); it++ ) {

    if ( !force && (*it)->hasGen() ) continue;
    if ( dipole.partons().second == (*it)->partons().first ||
	 (*it)->partons().second == dipole.partons().first )
      continue;

    InvEnergy2 a = dipole.partons().first->dist2(*dipole.partons().second);
    InvEnergy2 b = (*it)->partons().first->dist2(*((*it)->partons().second));
    InvEnergy2 c = dipole.partons().second->dist2(*(*it)->partons().first);
    InvEnergy2 d = dipole.partons().first->dist2(*(*it)->partons().second);

    //calculate recoils and kinematics
    if ( false ) //veto from kinematics
      continue;

    a = sqr(Current<DipoleEventHandler>()->rMax())/
      (Current<DipoleEventHandler>()->alphaS(sqrt(a)))*
      sqr(exp(sqrt(a)/Current<DipoleEventHandler>()->rMax()) - 1.0);
    b = sqr(Current<DipoleEventHandler>()->rMax())/
      (Current<DipoleEventHandler>()->alphaS(sqrt(b)))*
      sqr(exp(sqrt(b)/Current<DipoleEventHandler>()->rMax()) - 1.0);
    c = sqr(Current<DipoleEventHandler>()->rMax())/
      (Current<DipoleEventHandler>()->alphaS(sqrt(c)))*
      sqr(exp(sqrt(c)/Current<DipoleEventHandler>()->rMax()) - 1.0);
    d = sqr(Current<DipoleEventHandler>()->rMax())/
      (Current<DipoleEventHandler>()->alphaS(sqrt(d)))*
      sqr(exp(sqrt(d)/Current<DipoleEventHandler>()->rMax()) - 1.0);
    InvEnergy scale = sqrt(min(min(a, b),min(c, d)));
    double A = c*d/(a*b);
    double y = max(max(dipole.partons().first->y(), dipole.partons().second->y()), miny)
      - A*log(UseRandom::rnd())/sqr(Current<DipoleEventHandler>()->alphaS(scale));
    //decide at what rate this happens!!! alpha_S should be there...

    cout << "generated exchange: " << y << ", yGen " << dipole.generatedY() << endl;

    if ( y < dipole.generatedY() ) {
      cout << "generated exchange for y = " << y << endl;
      dipole.swingDipole(*it);
      dipole.generatedY(y);
      dipole.recoilSwing(true);
     }
  }
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void RecoilSwinger::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void RecoilSwinger::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<RecoilSwinger,DIPSY::Swinger>
  describeDIPSYRecoilSwinger("DIPSY::RecoilSwinger", "RecoilSwinger.so");

void RecoilSwinger::Init() {

  static ClassDocumentation<RecoilSwinger> documentation
    ("There is no documentation for the RecoilSwinger class");

}

