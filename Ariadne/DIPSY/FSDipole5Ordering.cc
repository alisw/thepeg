// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FSDipole5Ordering class.
//

#include "FSDipole5Ordering.h"
#include "Ariadne/Cascade/AriadneHandler.h"
#include "Ariadne/Cascade/Parton.h"
#include "Ariadne/Cascade/DipoleState.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

FSDipole5Ordering::FSDipole5Ordering()
  : ReweightBase(true), isGenerous(0), useOnlyOriginal(false), f(1.0),
    ptmin(ZERO), hardSuppression(0), fudge(1.0) {}

FSDipole5Ordering::~FSDipole5Ordering() {}

bool FSDipole5Ordering::
finalVeto(const Emission & emission) const {
  if ( emission.rho < ptmin || emission.partons.empty() ) return false;

  pair<tcParPtr,tcParPtr> parents =
    make_pair(emission.colourParent, emission.antiColourParent);

  if ( hardSuppression ) {
    LorentzRotation Rcm = Utilities::getBoostToCM(emission.pold);
    double alpha = ( hardSuppression == 1? Current<AriadneHandler>()->hardAlpha():
		     Current<AriadneHandler>()->softAlpha());
    double beta = Current<AriadneHandler>()->beta();
    Energy W = (emission.pold.first + emission.pold.second).m();
    double frac =
      (fudge*emission.rho*(Rcm*(parents.first->vertex() - parents.second->vertex())).perp())/hbarc;
    double fmax = pow(1.0 - (Rcm*parents.first->momentum()).plus()/W, -1.0/alpha);
    if ( frac > fmax && pow(fmax/frac, beta) < UseRandom::rnd() ) return true;
    fmax = pow(1.0 - (Rcm*parents.first->momentum()).minus()/W, -1.0/alpha);
    if ( frac > fmax && pow(fmax/frac, beta) < UseRandom::rnd() ) return true;
  }

  bool o1 = parents.first->orig();
  bool o3 = parents.second->orig();
  while ( !parents.first->orig() )
    parents.first = parents.first->emission()->colourParent;
  while ( !parents.second->orig() )
    parents.second = parents.second->emission()->antiColourParent;
  LorentzMomentum p1 = parents.first->orig()->momentum();
  LorentzMomentum p3 = parents.second->orig()->momentum();

  for ( int i = 0, N = emission.partons.size(); i < N; ++i ) {
    tcParPtr parton = emission.partons[i];
    LorentzMomentum p2 = parton->momentum();
    if ( useOnlyOriginal ) {
      if ( !o1 && !o3 ) return false;
      if ( o1 && ( p2.perp2() > f*p1.perp2() ||
		   p2.plus() > f*p1.plus() ||
		   p2.minus() > f*p1.minus() ) ) return true;
      if ( o3 && ( p2.perp2() > f*p3.perp2() ||
		   p2.plus() > f*p3.plus() ||
		   p2.minus() > f*p3.minus() ) ) return true;
      
      if ( o1 || o3 ) return false;
    }

    Energy2 pt2cut = min(p1.perp2(), p3.perp2());
    if ( isGenerous > 0 ) pt2cut = sqrt(p1.perp2()*p3.perp2());
    if ( isGenerous < 0 ) pt2cut = ZERO;
    if ( p2.perp2() < f*pt2cut ) return false;
    if ( p2.plus() < f*p1.plus() && p2.minus() < f*p1.minus() ) return false;
    if ( p2.plus() < f*p3.plus() && p2.minus() < f*p3.minus() ) return false;

  }
  return true;
}

IBPtr FSDipole5Ordering::clone() const {
  return new_ptr(*this);
}

IBPtr FSDipole5Ordering::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FSDipole5Ordering::persistentOutput(PersistentOStream & os) const {
  os << isGenerous << useOnlyOriginal << f << ounit(ptmin, GeV) << hardSuppression << fudge;
}

void FSDipole5Ordering::persistentInput(PersistentIStream & is, int) {
  is >> isGenerous >> useOnlyOriginal >> f >> iunit(ptmin, GeV) >> hardSuppression >> fudge;
}

DescribeClass<FSDipole5Ordering,Ariadne5::ReweightBase>
describeDIPSYFSDipole5Ordering("DIPSY::FSDipole5Ordering",
			       "FSDipole5Ordering.so");

void FSDipole5Ordering::Init() {

  static ClassDocumentation<FSDipole5Ordering> documentation
    ("There is no documentation for the FSDipole5Ordering class");


  static Switch<FSDipole5Ordering,int> interfaceGenerous
    ("Generous",
     "Determines whether the phase space constraints for new emissions is "
     "relaxed from the minimum of the transverse momenta of the emittors to "
     "the geometric mean. Alternatively the phase space can be more "
     "constrained giving pure light-cone momentum ordering.",
     &FSDipole5Ordering::isGenerous, 0, true, false);
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


  static Switch<FSDipole5Ordering,bool> interfaceOnlyOriginal
    ("OnlyOriginal",
     "Determines whether all emissions are restricted by the kinematics of "
     "the original partons, or only those where original partons are involved",
     &FSDipole5Ordering::useOnlyOriginal, false, true, false);
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


  static Parameter<FSDipole5Ordering,double> interfaceFudge
    ("Fudge",
     "Fudge parameter to allow less strict ordering if larger than unity.",
     &FSDipole5Ordering::f, 1.0, 0.0, 2.0,
     true, false, Interface::limited);

  static Parameter<FSDipole5Ordering,Energy> interfacePTMin
    ("PTMin",
     "Limit on invariant transverse momentum. Any emission below this is "
     "not checked.",
     &FSDipole5Ordering::ptmin, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     true, false, Interface::lowerlim);

  static Switch<FSDipole5Ordering,int> interfaceHardSuppression
    ("HardSuppression",
     "Options for suppression of hard emissions due to the distance between the partons in the radiating dipole.",
     &FSDipole5Ordering::hardSuppression, 0, true, false);
  static SwitchOption interfaceHardSuppressionOff
    (interfaceHardSuppression,
     "Off",
     "No suppression.",
     0);
  static SwitchOption interfaceHardSuppressionSuppression 
    (interfaceHardSuppression,
     "HardParameters",
     "Suppression using the parameters for hard remnants in the soft radiation model.",
     1);
  static SwitchOption interfaceHardSuppressionSoftParameters
    (interfaceHardSuppression,
     "SoftParameters",
     "Suppression using the parameters for soft remnants in the soft radiation model.",
     2);

  static Parameter<FSDipole5Ordering,double> interfaceHardFudge
    ("HardFudge",
     "Fudge factor multiplying transverse momentum of an emission when calculating hard suppression.",
     &FSDipole5Ordering::fudge, 1.0, 0.0, 0,
     true, false, Interface::lowerlim);

}

