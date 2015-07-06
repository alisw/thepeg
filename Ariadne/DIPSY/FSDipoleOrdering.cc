// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FSDipoleOrdering class.
//

#include "FSDipoleOrdering.h"
#include "Ariadne/DipoleCascade/Parton.h"
#include "Ariadne/DipoleCascade/DipoleState.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

FSDipoleOrdering::FSDipoleOrdering()
  : isGenerous(0), useOnlyOriginal(false) {}

FSDipoleOrdering::~FSDipoleOrdering() {}

bool FSDipoleOrdering::
finalVeto(tcEmiPtr dipole, Ariadne5::tcDipoleStatePtr state, tcParPtr parton,
	  Energy2 pt2, const EmissionType & type) const {
  pair<tcParPtr,tcParPtr> parents = parton->parents();
  bool o1 = parents.first->orig();
  bool o3 = parents.second->orig();
  while ( !parents.first->orig() )
    parents.first = parents.first->parents().first;
  while ( !parents.second->orig() )
    parents.second = parents.second->parents().second;
  LorentzMomentum p1 = parents.first->orig()->momentum();
  LorentzMomentum p2 = parton->momentum();
  LorentzMomentum p3 = parents.second->orig()->momentum();
  if ( useOnlyOriginal ) {
    if ( !o1 && !o3 ) return false;
    if ( o1 && ( p2.perp2() > p1.perp2() ||
		 p2.plus() > p1.plus() ||
		 p2.minus() > p1.minus() ) ) return true;
    if ( o3 && ( p2.perp2() > p3.perp2() ||
		 p2.plus() > p3.plus() ||
		 p2.minus() > p3.minus() ) ) return true;
							      
    if ( o1 || o3 ) return false;
  }

  Energy2 pt2cut = min(p1.perp2(), p3.perp2());
  if ( isGenerous > 0 ) pt2cut = sqrt(p1.perp2()*p3.perp2());
  if ( isGenerous < 0 ) pt2cut = ZERO;
  if ( p2.perp2() < pt2cut ) return false;
  if ( p2.plus() < p1.plus() && p2.minus() < p1.minus() ) return false;
  if ( p2.plus() < p3.plus() && p2.minus() < p3.minus() ) return false;

  return true;
}

IBPtr FSDipoleOrdering::clone() const {
  return new_ptr(*this);
}

IBPtr FSDipoleOrdering::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FSDipoleOrdering::persistentOutput(PersistentOStream & os) const {
  os << isGenerous << useOnlyOriginal;
}

void FSDipoleOrdering::persistentInput(PersistentIStream & is, int) {
  is >> isGenerous >> useOnlyOriginal;
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<FSDipoleOrdering,Ariadne::ReweightBase>
  describeDIPSYFSDipoleOrdering("DIPSY::FSDipoleOrdering",
				"FSDipoleOrdering.so");

void FSDipoleOrdering::Init() {

  static ClassDocumentation<FSDipoleOrdering> documentation
    ("There is no documentation for the FSDipoleOrdering class");


  static Switch<FSDipoleOrdering,int> interfaceGenerous
    ("Generous",
     "Determines whether the phase space constraints for new emissions is "
     "relaxed from the minimum of the transverse momenta of the emittors to "
     "the geometric mean. Alternatively the phase space can be more "
     "constrained giving pure light-cone momentum ordering.",
     &FSDipoleOrdering::isGenerous, 0, true, false);
  static SwitchOption interfaceGenerousNormal
    (interfaceGenerous,
     "Normal",
     "The transverse momentum of a gluon is restricted to be less than any "
     "of its mothers, if outside their light-cone triangle.",
     0);
  static SwitchOption interfaceGenerousGenerous
    (interfaceGenerous,
     "Generous",
     "The transverse momentum of a gluon is restricted to be less than the "
     "geometrical mean of that of its original parents, if outside their "
     "light-cone triangle.",
     1);
  static SwitchOption interfaceGenerousRestrictive
    (interfaceGenerous,
     "Restrictive",
     "No gluons are allowed outside the light-cone triangle of the "
     "original parents.",
     -1);


  static Switch<FSDipoleOrdering,bool> interfaceOnlyOriginal
    ("OnlyOriginal",
     "Determines whether all emissions are restricted by the kinematics of "
     "the original partons, or only those where original partons are involved",
     &FSDipoleOrdering::useOnlyOriginal, false, true, false);
  static SwitchOption interfaceOnlyOriginalAllEmissions
    (interfaceOnlyOriginal,
     "AllEmissions",
     "All emissions are restricted by the kinematics of the original partons.",
     false);
  static SwitchOption interfaceOnlyOriginalOnlyOriginalEmissions
    (interfaceOnlyOriginal,
     "OnlyOriginalEmissions",
     "Emissions are restricted by the kinematics of the original partons only "
     "if original partons are involved in the emission.",
     true);


}

