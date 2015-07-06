// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleAnalysisHandler class.
//

#include "DipoleAnalysisHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"



using namespace DIPSY;

DipoleAnalysisHandler::DipoleAnalysisHandler() {}

DipoleAnalysisHandler::~DipoleAnalysisHandler() {}


void DipoleAnalysisHandler::analyze(const DipoleState &, const DipoleState &,
				    const ImpactParameters &, const DipoleXSec &,
				    double, CrossSection) {}


void DipoleAnalysisHandler::analyze(const vector<DipoleStatePtr> &,
				    const vector<DipoleStatePtr> &,
				    const vector<ImpactParameters> &, const DipoleXSec &,
				    const Vec3D & probs, double) {}


ostream & DipoleAnalysisHandler::stub(string text) const {
    generator()->log().setf(ios::left, ios::adjustfield);
    generator()->log() << setw(50) << name() +text;
    return generator()->log();
}



// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).



// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeAbstractNoPIOClass<DipoleAnalysisHandler,HandlerBase>
describeDIPSYDipoleAnalysisHandler("DIPSY::DipoleAnalysisHandler",
				   "libAriadne5.so libDIPSY.so");

void DipoleAnalysisHandler::Init() {

  static ClassDocumentation<DipoleAnalysisHandler> documentation
    ("The DipoleAnalysisHandler class can be used as a base class for "
     "sub-classes which can be used in the initialization of a "
     "DipoleEventHandler, where statistics on, eg. total and elastic "
     "cross sections, can be collected during the presampling phase.");

}

