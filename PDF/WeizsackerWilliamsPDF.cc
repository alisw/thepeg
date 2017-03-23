// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WeizsackerWilliamsPDF class.
//

#include "WeizsackerWilliamsPDF.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"

using namespace ThePEG;

WeizsackerWilliamsPDF::WeizsackerWilliamsPDF() 
  : _q2min(ZERO), _q2max(4.*GeV2), _a(0.8)
{}

bool WeizsackerWilliamsPDF::canHandleParticle(tcPDPtr particle) const {
  return ( abs(particle->id()) == ParticleID::eminus ||
	   abs(particle->id()) == ParticleID::muminus );
}

cPDVector WeizsackerWilliamsPDF::partons(tcPDPtr) const {
  // only photon
  return cPDVector(1,getParticleData(ParticleID::gamma));
}

double WeizsackerWilliamsPDF::xfl(tcPDPtr, tcPDPtr parton, Energy2,
                      double l, Energy2 ) const {
  if(parton->id()!=ParticleID::gamma) return 0.;
  double x(exp(-l));
  // remember this is x fx so no x in the denominator!!
  return 0.5*SM().alphaEM()/Constants::pi*(1.+sqr(1.-x));
}

double WeizsackerWilliamsPDF::xfvl(tcPDPtr, tcPDPtr, Energy2, double,
				   Energy2) const {
  // valence density is zero
  return 0.0;
}


ClassDescription<WeizsackerWilliamsPDF> 
WeizsackerWilliamsPDF::initWeizsackerWilliamsPDF;
// Definition of the static class description member.

void WeizsackerWilliamsPDF::Init() {

  static ClassDocumentation<WeizsackerWilliamsPDF> documentation
    ("The WeizsackerWilliamsPDF provides the PDF for a photon inside"
     " an incoming lepton in the Weisszacker-Williams approximation");

  static Parameter<WeizsackerWilliamsPDF,Energy2> interfaceQ2Min
    ("Q2Min",
     "Minimum value of the magnitude of Q^2 for the photon",
     &WeizsackerWilliamsPDF::_q2min, GeV2, ZERO, ZERO, 100.0*GeV2,
     false, false, Interface::limited);

  static Parameter<WeizsackerWilliamsPDF,Energy2> interfaceQ2Max
    ("Q2Max",
     "Maximum value of the magnitude of Q^2 for the photon",
     &WeizsackerWilliamsPDF::_q2max, GeV2, 4.0*GeV2, ZERO, 100.0*GeV2,
     false, false, Interface::limited);


  static Parameter<WeizsackerWilliamsPDF,double> interfaceJacobianParameter
    ("JacobianParameter",
     "Relative contribution for the 1/x and flat term in sampling x ",
     &WeizsackerWilliamsPDF::_a, 0.8, 0.0, 1.0,
     false, false, Interface::limited);

}

double WeizsackerWilliamsPDF::
flattenScale(tcPDPtr a, tcPDPtr, const PDFCuts & c,
	     double l, double z, double & jacobian) const {
  double x = exp(-l);
  Energy2 qqmax = min(_q2max,0.25*sqr(x)*c.sMax());
  Energy2 qqmin = max(_q2min,sqr(x*a->mass())/(1.-x));
  if(qqmin>=qqmax) {
    jacobian = 0.;
    return 0.;
  }
  double low = log(qqmin/c.scaleMaxL(l));
  double upp = log(qqmax/c.scaleMaxL(l));
  double ret = exp(low+z*(upp-low));
  // jacobian factor
  jacobian *= log(qqmax/qqmin);
  return ret;
}

double WeizsackerWilliamsPDF::flattenL(tcPDPtr, tcPDPtr, const PDFCuts & c,
				       double z, double & jacobian) const {
  double xmin = exp(-c.lMax());
  double xmax = exp(-c.lMin());
  double x;
  if(z<=_a) {
    z /=_a;
    x = exp(-(c.lMin() + z*(c.lMax() - c.lMin())));
  }
  else {
    z = (z-_a)/(1.-_a);
    x = xmin+z*(xmax-xmin);
  }
  jacobian /= _a/(c.lMax() - c.lMin()) + (1.-_a)*x/(xmax-xmin);
  return -log(x);
}

void WeizsackerWilliamsPDF::persistentOutput(PersistentOStream & os) const {
  os << ounit(_q2min,GeV2) << ounit(_q2max,GeV2) << _a;
}

void WeizsackerWilliamsPDF::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_q2min,GeV2) >> iunit(_q2max,GeV2) >> _a;
}

