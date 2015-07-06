// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TotalXSecAnalysis class.
//

#include "TotalXSecAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "DipoleXSec.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

TotalXSecAnalysis::TotalXSecAnalysis(): sum(0.0*picobarn), sum2(0.0*sqr(picobarn)) {}

TotalXSecAnalysis::~TotalXSecAnalysis() {}

void TotalXSecAnalysis::initialize() {
  sum = 0.0*picobarn;
  sum2 = sqr(sum);
}

void TotalXSecAnalysis::
analyze(const DipoleState & dl, const DipoleState & dr,
	const ImpactParameters & b, const DipoleXSec & xsec,
	double fsum, CrossSection weight) {
  sum += 2.0*xsec.unitarize(fsum)*weight;
  sum2 += sqr(2.0*xsec.unitarize(fsum)*weight);
}

void TotalXSecAnalysis::finalize(long neve) {
  if ( neve <= 0 ) return;
  sum /= neve;
  sum2 /= neve;
  CrossSection err = sqrt((sum2 - sqr(sum))/neve);
  generator()->log().setf(ios::left, ios::adjustfield);
  generator()->log() << setw(50) << name() + ": Total cross section: "
		     << ouniterr(sum, err, nanobarn) << " nb." << endl;
}

IBPtr TotalXSecAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr TotalXSecAnalysis::fullclone() const {
  return new_ptr(*this);
}


void TotalXSecAnalysis::persistentOutput(PersistentOStream & os) const {
  os << ounit(sum, picobarn) << ounit(sum2, sqr(picobarn));
}

void TotalXSecAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> iunit(sum, picobarn) >> iunit(sum2, sqr(picobarn));
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<TotalXSecAnalysis,DIPSY::DipoleAnalysisHandler>
  describeDIPSYTotalXSecAnalysis("DIPSY::TotalXSecAnalysis", "libAriadne5.so libDIPSY.so");


void TotalXSecAnalysis::Init() {

  static ClassDocumentation<TotalXSecAnalysis> documentation
    ("TotalXSecAnalysis class calculates the total cross section in the "
     "presampling phase of the DipoleEventHandler.");

}

