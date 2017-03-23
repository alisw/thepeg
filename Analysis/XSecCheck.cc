// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the XSecCheck class.
//

#include "XSecCheck.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/EventRecord/Event.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "XSecCheck.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

XSecCheck::~XSecCheck() {}

IBPtr XSecCheck::clone() const {
  return new_ptr(*this);
}

IBPtr XSecCheck::fullclone() const {
  return new_ptr(*this);
}

void XSecCheck::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  sumw += event->weight();
}

void XSecCheck::dofinish() {
  AnalysisHandler::dofinish();
  CrossSection xsec = sumw*generator()->histogramScale();
  if ( abs(xsec - target) > tol*abs(xsec + target) )
    Throw<UnexpectedXSec>()
      << "The total cross section of this run, " << xsec/picobarn
      << " picobarn, dit not match the target value, " << target/picobarn
      << " picobarn." << Exception::warning;
}

void XSecCheck::persistentOutput(PersistentOStream & os) const {
  os << ounit(target, picobarn) << tol << sumw;
}

void XSecCheck::persistentInput(PersistentIStream & is, int) {
  is >> iunit(target, picobarn) >> tol >> sumw;
}

ClassDescription<XSecCheck> XSecCheck::initXSecCheck;
// Definition of the static class description member.

void XSecCheck::Init() {

  static ClassDocumentation<XSecCheck> documentation
    ("The XSecCheck class is a simple analysis class used for testing "
     "purposes. If the total cross section does not match the one specified "
     "by <interface>TargetXSec</interface>, an exception will be thrown.");

  static Parameter<XSecCheck,CrossSection> interfaceTargetXSec
    ("TargetXSec",
     "The expected total cross section in units of picobarn.",
     &XSecCheck::target, picobarn, ZERO, ZERO, ZERO,
     true, false, Interface::lowerlim);
  interfaceTargetXSec.setHasDefault(false);

  static Parameter<XSecCheck,double> interfaceTolerance
    ("Tolerance",
     "The relative tolerance accepted when comparing the total cross "
     "section with <interface>TargetXSec</interface>.",
     &XSecCheck::tol, 0.01, 0.0, 1.0,
     true, false, Interface::limited);
  interfaceTolerance.setHasDefault(false);

}

