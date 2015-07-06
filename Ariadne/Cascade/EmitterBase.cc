// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EmitterBase class.
//

#include "EmitterBase.h"
#include "Emission.h"
#include "QCDDipole.h"
#include "Parton.h"
#include "Junction.h"
#include "RemnantParton.h"
#include "ResonanceParton.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Repository/UseRandom.h"
#include "AriadneHandler.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"


using namespace Ariadne5;

/*

  We will assume that an external alphas is of second order where the
  first order is

  alphas_s(mu)=1/(b0 x)  [with x = log(mu^2/Lambda^2)]

  giving the second order

  alpha_s(mu)= 1/(b0 x){1 - (b1/(b0**2))*(log(x)/x) +
              ((b1**2)*(log(x)**2-log(x)-1) + b0*b2)/(b0*x)**2}

  with
  b0 = (33-2*nf)/(12*pi)
  b1 = (153-19*nf)/(24*pi**2)
  b2 = (2857 - 5033*nf/9 + 325*(nf**2)/27)/(128*pi**3)

  the second order effects are approximated by a shift in Lambda and a factor

  alpha1_s(mu)=C/(b0*(x-shift))

  as2(x1)==C/(b0*(x1-shift))
  as2(x2)==C//b0*(x2-shift))

  as2(x1)/as2(x2)==(x1-shift)/(x2-shift)

 */

Energy EmitterBase::rhoCut() const {
  return Current<AriadneHandler>()->pTCut();
}


bool EmitterBase::touched(const DipoleBase & dipole) const {
  if ( dipole.touched() ) return true;
  if ( tcQCDPtr d = dynamic_ptr_cast<tcQCDPtr>(&dipole) )
    if ( d->iPart()->touched() || d->oPart()->touched() ) return true;
  return false;
}

Energy2 EmitterBase::rndsud(double C, Energy2 pt2max, Energy2 pt2min) {
  if ( Current<AriadneHandler>()->runningCoupling() !=
       AriadneHandler::noRunning ) return runrndsud(C, pt2max, pt2min);
  double CN = 1.0/(C*Current<AriadneHandler>()->alpha0());
  double R = UseRandom::rnd();
  return CN*log(R) < log(pt2min/pt2max)? -1.0*GeV2: pt2max*pow(R, CN);
}

Energy2 EmitterBase::
runrndsud(double C, Energy2 pt2max, Energy2 pt2min) {

  // Check which flavour thresholds we are between.
  set<Energy2>::iterator fit =
    Current<AriadneHandler>()->flavourThresholds().lower_bound(pt2max);

  Energy2 pt2cut = pt2min;
  if ( fit != Current<AriadneHandler>()->flavourThresholds().begin() )
    pt2cut = max(*(--fit), pt2min);

  Energy2 pt2m = sqrt(pt2max*pt2cut);
  Energy2 pt2u = sqrt(pt2max*pt2m);
  Energy2 pt2l = sqrt(pt2cut*pt2m);

  // Sample alpha0 away from flavour thresholds to be on the safe side.
  double aup = Current<AriadneHandler>()->alphaS(pt2u);
  double alo = Current<AriadneHandler>()->alphaS(pt2l);
  if ( abs(aup - alo) < 1.0e-12 ) {
    pt2u *= 2.0;
    aup = Current<AriadneHandler>()->alphaS(pt2u);
  }
  double a0 = aup*alo*log(pt2l/pt2u)/(aup - alo);
  Energy2 Lam2 = pt2u/exp(a0/aup);
  // Add a little margin to the derived alpha0 to avoid rounding errors.
  double CN = 0.9/(C*a0);
  double logpt2min = log(pt2min/Lam2);

  while ( true ) {
    double RtoCN = pow(UseRandom::rnd(), CN);
    if ( RtoCN <= logpt2min/log(pt2max/Lam2) )
      pt2max = ZERO;
    else
      pt2max = Lam2*pow(pt2max/Lam2, RtoCN);
    if ( pt2max < pt2cut )
      return pt2cut > pt2min? runrndsud(C, pt2cut, pt2min): pt2max;
    double weight = Current<AriadneHandler>()->alphaS(pt2max)*
      log(pt2max/Lam2)*0.9/a0;
    if ( weight > 1.0 )
      Throw<WeightException>()
	<< "The Ariadne5::CascadeHandler '"
	<< Current<AriadneHandler>()->name()
	<< "' failed to overestimate the alpha_S specified by the StandardModel "
	<< "object. If this hapens too often you should contact the authors."
	<< Exception::warning;			   
    if ( UseRandom::rnd() < weight ) return pt2max;
  }
}

bool EmitterBase::
check(double x1, double x3, double y1, double y2, double y3) {
  double x2 = 2.0 - x1 - x3;
  if ( x1 < 0.0 || x2 < 0.0 || x3 < 0.0 ) return false;
  double pp1 = sqr(0.5*x1) - y1;
  double pp2 = sqr(0.5*x2) - y2;
  double pp3 = sqr(0.5*x3) - y3;
  if ( pp1 < 0.0 || pp2 < 0.0 || pp3 < 0.0 ||
       2.0*(pp1*pp2 + pp2*pp3 + pp3*pp1) <= sqr(pp1) + sqr(pp2) + sqr(pp3) )
    return false;
  return true;
}

double EmitterBase::
preweight(const Emission & emission) const {
  double w = 1.0;
  for ( int i = 0, N = theReweighters.size(); i < N; ++i )
    w *= theReweighters[i]->preweight(emission);
  const vector<DipoleRWPtr> & auxw =
    Current<AriadneHandler>()->reweighters();
  for ( int i = 0, N = auxw.size(); i < N; ++i )
    w *= auxw[i]->preweight(emission);
  return w;
}

double EmitterBase::
reweight(const Emission & emission) const {
  double w = 1.0;
  for ( int i = 0, N = theReweighters.size(); i < N; ++i )
    w *= theReweighters[i]->reweight(emission);
  const vector<DipoleRWPtr> & auxw =
    Current<AriadneHandler>()->reweighters();
  for ( int i = 0, N = auxw.size(); i < N; ++i )
    w *= auxw[i]->reweight(emission);
  return w;
}

bool EmitterBase::hasFinalVeto() const {
  for ( int i = 0, N = theReweighters.size(); i < N; ++i )
    if ( theReweighters[i]->hasFinalVeto() ) return true;
  const vector<DipoleRWPtr> & auxw =
    Current<AriadneHandler>()->reweighters();
  for ( int i = 0, N = auxw.size(); i < N; ++i )
    if ( auxw[i]->hasFinalVeto() ) return true;
  return false;
}

bool EmitterBase::
finalVeto(const Emission & emission) const {
  for ( int i = 0, N = theReweighters.size(); i < N; ++i )
    if ( theReweighters[i]->hasFinalVeto() &&
	 theReweighters[i]->finalVeto(emission) ) return true;
  const vector<DipoleRWPtr> & auxw =
    Current<AriadneHandler>()->reweighters();
  for ( int i = 0, N = auxw.size(); i < N; ++i )
    if ( auxw[i]->hasFinalVeto() &&
	 auxw[i]->finalVeto(emission) ) return true;
  return false;
}

Energy2 EmitterBase::invPT2(tcParPtr p1, tcParPtr p2, tcParPtr p3) {
  Energy2 s123 = (p1->momentum() + p2->momentum() + p3->momentum()).m2();
  Energy2 s12 = (p1->momentum() + p2->momentum()).m2();
  Energy2 s23 = (p2->momentum() + p3->momentum()).m2();
  Energy m1 = p1->momentum().mass();
  Energy m2 = p2->momentum().mass();
  Energy m3 = p3->momentum().mass();
  return (s12 - sqr(m1 + m2)) * (s23 - sqr(m2 + m3)) / s123;
}

double EmitterBase::invY(tcParPtr p1, tcParPtr p2, tcParPtr p3) {
  Energy2 s12 = (p1->momentum() + p2->momentum()).m2();
  Energy2 s23 = (p1->momentum() + p2->momentum()).m2();
  Energy m1 = p1->momentum().mass();
  Energy m2 = p2->momentum().mass();
  Energy m3 = p3->momentum().mass();
  return 0.5*log((s12 - sqr(m1 + m2))/(s23 - sqr(m2 + m3)));
}

Triplet<Lorentz5Momentum,Lorentz5Momentum,Lorentz5Momentum>
EmitterBase::getMomenta(Energy2 s, double x1, double x3,
			Energy m1, Energy m2, Energy m3,
			bool nr1, bool nr3, bool userho) {
  Trip p =
    Trip(Lorentz5Momentum(m1), Lorentz5Momentum(m2), Lorentz5Momentum(m3));
  SimplePhaseSpace::CMS(p.first, p.second, p.third, s, x1, x3, 0.0, 0.0, 0.0);
  double Psi = Constants::pi - p.third.theta();
  double beta = 0.0;
  if ( userho ) {
    x1 = p.first.rho()/sqrt(s);
    x3 = p.third.rho()/sqrt(s);
  }
  if ( nr1 && nr3 ) beta = Psi*sqr(x3)/(sqr(x1) + sqr(x3)); // minimize pt
  else if ( nr3 ) beta = Psi;
  else if ( !nr1 && sqr(x3) > UseRandom::rnd()*(sqr(x1) + sqr(x3)) )
    beta = Psi;
  LorentzRotation R;
  R.rotateY(-beta);
  p.first.transform(R);
  p.second.transform(R);
  p.third.transform(R);
  return p;
}

Triplet<Lorentz5Momentum,Lorentz5Momentum,Lorentz5Momentum>
EmitterBase::getMomenta(Energy2 s, double x1, double x3,
			Energy m1, Energy m2, Energy m3,
			bool nr1, bool nr3, double phi, bool userho) {
  Trip p =
    Trip(Lorentz5Momentum(m1), Lorentz5Momentum(m2), Lorentz5Momentum(m3));
  SimplePhaseSpace::CMS(p.first, p.second, p.third, s, x1, x3, 0.0, 0.0, 0.0);
  double Psi = Constants::pi - p.third.theta();
  double beta = 0.0;
  if ( userho ) {
    x1 = p.first.rho()/sqrt(s);
    x3 = p.third.rho()/sqrt(s);
  }
  if ( nr1 && nr3 ) beta = Psi*sqr(x3)/(sqr(x1) + sqr(x3)); // minimize pt
  else if ( nr3 ) beta = Psi;
  else if ( !nr1 && sqr(x3) > UseRandom::rnd()*(sqr(x1) + sqr(x3)) )
    beta = Psi;
  LorentzRotation R;
  R.rotateY(-beta);
  R.rotateZ(phi);
  p.first.transform(R);
  p.second.transform(R);
  p.third.transform(R);
  return p;
}

Triplet<Lorentz5Momentum,Lorentz5Momentum,Lorentz5Momentum>
EmitterBase::
getMomenta(Energy2 s, double x1, double x3,
	   Energy m1, Energy m2, Energy m3, bool nr1, bool nr3, double phi,
	   const Lorentz5Momentum & p1, const Lorentz5Momentum & p3,
	   bool userho) {
  Trip p =
    Trip(Lorentz5Momentum(m1), Lorentz5Momentum(m2), Lorentz5Momentum(m3));
  SimplePhaseSpace::CMS(p.first, p.second, p.third, s, x1, x3, 0.0, 0.0, 0.0);
  double Psi = Constants::pi - p.third.theta();
  double beta = 0.0;
  if ( userho ) {
    x1 = p.first.rho()/sqrt(s);
    x3 = p.third.rho()/sqrt(s);
  }
  if ( nr1 && nr3 ) beta = Psi*sqr(x3)/(sqr(x1) + sqr(x3)); // minimize pt
  else if ( nr3 ) beta = Psi;
  else if ( !nr1 && sqr(x3) > UseRandom::rnd()*(sqr(x1) + sqr(x3)) ) beta = Psi;
  LorentzRotation R;
  R.rotateY(-beta);
  R.rotateZ(phi);
  R.transform(Utilities::getBoostFromCM(make_pair(p1, p3)));
  p.first.transform(R);
  p.second.transform(R);
  p.third.transform(R);
  return p;
}

void EmitterBase::revert(const Emission & emission) const {
  Throw<Exception>() 
    << "The emitter model " << fullName()
    << " has not implemented a revert() method even though it reports "
    << "that e performed emission was reversible. Please contact the author "
    << "to have this error corrected." << Exception::runerror;
}

vector<EmPtr> EmitterBase::inverseEmissions(const DipoleState &) const {
  vector<EmPtr> ret;
  return ret;
}

bool EmitterBase::overrideInverse(const Emission &) const {
  return false;
}

bool EmitterBase::performInverse(const Emission &, DipoleState &) const {
  return false;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void EmitterBase::persistentOutput(PersistentOStream & os) const {
  os << theReweighters << canReconstruct << willReconstruct;
}

void EmitterBase::persistentInput(PersistentIStream & is, int) {
  is >> theReweighters >> canReconstruct >> willReconstruct;
}

void EmitterBase::setReconstruct(bool t) {
  if ( t && ! canReconstruct )
    Throw<InterfaceException>()
      << name() << " cannot be used to reconstruct emissions."
      << Exception::warning;
  willReconstruct = t;
}

bool EmitterBase::defReconstruct() const {
  return canReconstruct;
}

DescribeAbstractClass<EmitterBase,HandlerBase>
describeAriadne5EmitterBase("Ariadne5::EmitterBase", "libAriadne5.so");

void EmitterBase::Init() {

  static ClassDocumentation<EmitterBase> documentation
    ("EmitterBase is the base class of all Ariadne classes implementing "
     "a specific model for emission from different kinds of dipoles. A "
     "sub-class must implement a sub-class of Emission to store the "
     "result of an emission.");

  static RefVector<EmitterBase,Ariadne5::ReweightBase> interfaceReweighters
    ("Reweighters",
     "A vector of objects implementing reweightings of basic dipole "
     "emissions. Each emission modelled by this emission model "
     "will be reweighted. Also global reweightors given in the "
     "Ariadne5::CascadeHandler will be included.",
     &EmitterBase::theReweighters, -1, true, false, true, false, false);


  static Switch<EmitterBase,bool> interfaceReconstruct
    ("Reconstruct",
     "Determines whether this emitter will be used to reconstruct emissions in the CKKW-L algorithm. May not be possible for all emitters.",
     &EmitterBase::willReconstruct, false, true, false,
     &EmitterBase::setReconstruct,
     (bool(EmitterBase::*)()const)(0), &EmitterBase::defReconstruct);
  static SwitchOption interfaceReconstructon
    (interfaceReconstruct,
     "on",
     "This emitter will be used for CKKW-L reconstruction.",
     true);
  static SwitchOption interfaceReconstructOff
    (interfaceReconstruct,
     "Off",
     "This emitter will not be used for CKKW-L reconstruction.",
     false);

}

