// -*- C++ -*-
//
// GRVBase.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GRVBase class.
//

#include "GRVBase.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Utilities/Maths.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;

GRVBase::GRVBase()
  : theLx(-1.0), thex(-1.0), theEps(-1.0), theRootx(-1.0),
    Q2(-GeV2), theLam2(-GeV2), theMu2(-GeV2),
    theS(-1.0), theS2(-1.0), theS3(-1.0), theRootS(-1.0),
    uvSave(-1.0), dvSave(-1.0), delSave(-1.0), udbSave(-1.0), sbSave(-1.0),
    cbSave(-1.0), bbSave(-1.0), glSave(-1.0) {}

GRVBase::~GRVBase() {}

bool GRVBase::canHandleParticle(tcPDPtr particle) const {
  return ( abs(particle->id()) == abs(long(ParticleID::pplus)) ||
	   abs(particle->id()) == abs(long(ParticleID::n0)) );
}

cPDVector GRVBase::partons(tcPDPtr p) const {
  cPDVector ret;
  if ( canHandleParticle(p) ) {
    ret.push_back(getParticleData(ParticleID::g));
    for ( int i = 1; i <= 5; ++i ) {
      ret.push_back(getParticleData(i));
      ret.push_back(getParticleData(-i));
    }
  }
  return ret;
}

double GRVBase::xfl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		    double l, Energy2) const {
  setup(l, partonScale);
  if ( S() < 0.0 ) return 0.0;
  using namespace ParticleID;
  bool anti = particle->id() < 0;
  bool neutron = abs(particle->id()) == n0;
  switch ( parton->id() ) {
  case b:
  case bbar:
    return max(fbb(), 0.0);
  case c:
  case cbar:
    return max(fcb(), 0.0);
  case s:
  case sbar:
    return max(fsb(), 0.0);
  case u:
    return max(neutron? (fudb() + fdel() + (anti? 0.0: fdv())):
	       (fudb() - fdel() + (anti? 0.0: fuv())), 0.0);
  case ubar:
    return max(neutron? (fudb() + fdel() + (anti? fdv(): 0.0)):
	       (fudb() - fdel() + (anti? fuv(): 0.0)), 0.0);
  case d:
    return max(neutron? (fudb() - fdel() + (anti? 0.0: fuv())):
	       (fudb() + fdel() + (anti? 0.0: fdv())), 0.0);
  case dbar:
    return max(neutron? (fudb() - fdel() + (anti? fuv(): 0.0)):
	       (fudb() + fdel() + (anti? fdv(): 0.0)), 0.0);
  case g:
    return max(fgl(), 0.0);
  }
  return 0.0;
}

double GRVBase::xfvl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double l, Energy2) const {
  setup(l, partonScale);
  if ( S() < 0.0 ) return 0.0;
  using namespace ParticleID;
  bool anti = particle->id() < 0;
  bool neutron = abs(particle->id()) == n0;
  switch ( parton->id() ) {
  case u:
    return max(neutron? (anti? 0.0: fdv()): (anti? 0.0: fuv()), 0.0);
  case ubar:
    return max(neutron? (anti? fdv(): 0.0): (anti? fuv(): 0.0), 0.0);
  case d:
    return max(neutron? (anti? 0.0: fuv()): (anti? 0.0: fdv()), 0.0);
  case dbar:
    return max(neutron? (anti? fuv(): 0.0): (anti? fdv(): 0.0), 0.0);
  }
  return 0.0;
}

double GRVBase::valens(double N, double ak, double bk,
		       double a, double b, double c, double d) const {
  return N*pow(x(), ak)*pow(eps(), d)*
    (1.0 + a*pow(x(), bk) + x()*(b + c*sqrt(x())));
}

double GRVBase::
lightsea(double al, double be, double ak, double bk,
	 double a, double b, double c, double d, double e, double es) const {
  return (pow(x(), ak)*(a + x()*(b + x()*c))*pow(lx(), bk) +
	  pow(S(), al)*exp(sqrt(es*pow(S(), be)*lx()) - e))*pow(eps(), d);
}

double GRVBase::
heavysea(double sth, double al, double be, double ak, double ag,
	 double b, double d, double e, double es) const {
  return S() <= sth? 0.0:
    pow(S() - sth, al)*(1.0 + rootx()*ag + x()*b)*pow(eps(), d)*
    exp(sqrt(es*pow(S(), be)*lx()) - e)/pow(lx(),ak);
}

void GRVBase::
setup(double l, Energy2 scale, Energy2 mu2, Energy2 lam2) const {
  if ( l == lx() && scale == Q2 && mu2 == theMu2 && lam2 == theLam2 ) return;
  if ( l != lx() || scale != Q2 || mu2 != theMu2 || lam2 != theLam2 ) {
    uvSave = dvSave = delSave = udbSave = -1.0;
    sbSave = cbSave = bbSave = glSave = -1.0;
  }
  if ( l != lx() ) {
    theLx = l;
    thex = exp(-l);
    if ( l < 0.0 ) throw PDFRange(name(), "momentum fraction", thex, 1.0);
    theEps = Math::exp1m(-l);
    theRootx = sqrt(x());
  }
  if ( scale != Q2 || mu2 != theMu2 || lam2 != theLam2 ) {
    Q2 = scale;
    theMu2 = mu2;
    theLam2 = lam2;
    if ( scale <= mu2 ) {
      switch ( rangeException ) {
      case rangeThrow:
	throw PDFRange(name(), "scale (in GeV^2)", scale/GeV2, mu2/GeV2);
      case rangeZero:
	theS = -1.0;
	return;
      case rangeFreeze:
	theS = 0.0;
      }
    } else
      theS = log(log(scale/lam2)/log(mu2/lam2));
    theS2 = sqr(S());
    theS3 = S()*S2();
    theRootS = sqrt(S());
  }
}

DescribeAbstractNoPIOClass<GRVBase,PDFBase>
describeGRVBase("ThePEG::GRVBase", "GRVBase.so");

void GRVBase::Init() {

  static ClassDocumentation<GRVBase> documentation
    ("This is the base class used by different GRV PDF parameterizations.");

}

