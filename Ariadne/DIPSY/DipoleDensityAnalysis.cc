// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleDensityAnalysis class.
//

#include "DipoleDensityAnalysis.h"
#include "DipoleXSec.h"
#include "DipoleState.h"
#include "DipoleEventHandler.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

DipoleDensityAnalysis::DipoleDensityAnalysis() {}

DipoleDensityAnalysis::~DipoleDensityAnalysis() {}

void DipoleDensityAnalysis::initialize() {
  generator()->histogramFactory()->initrun();
  generator()->histogramFactory()->registerClient(this);
  dipsizes = generator()->histogramFactory()->createHistogram1D
    ("dipsizes",60, -3.0, 3.0);
  gluonpts = generator()->histogramFactory()->createHistogram1D
    ("gluonpts",30, 0.0, 3.0);
  sumw = 0.0;
}

void DipoleDensityAnalysis::
analyze(const DipoleState & dl, const DipoleState & dr,
	const ImpactParameters & b, const DipoleXSec & xsec,
	double fsum, CrossSection weight) {

  fill(dl);
  fill(dr);

}

void DipoleDensityAnalysis::fill(const DipoleState & ds) {

  sumw += ds.weight();
  vector<tDipolePtr> final;
  ds.extract(back_inserter(final));
  for ( int i = 0, N = final.size(); i < N; ++i ) {
    double size = final[i]->size()/Current<DipoleEventHandler>()->rMax();
    dipsizes->fill(log10(size), ds.weight());
    double pt = final[i]->partons().first->pT().pt()/GeV;
    gluonpts->fill(log10(pt), ds.weight());
    if ( final[i]->neighbors().second ) {
      pt = final[i]->partons().second->pT().pt()/GeV;
      gluonpts->fill(log10(pt), ds.weight());
    }
  }

}

void DipoleDensityAnalysis::finalize(long neve) {
  if ( sumw <= 0.0 ) return;

  dipsizes->scale(10.0/sumw);
  gluonpts->scale(10.0/sumw);

}

IBPtr DipoleDensityAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr DipoleDensityAnalysis::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void DipoleDensityAnalysis::persistentOutput(PersistentOStream & os) const {
}

void DipoleDensityAnalysis::persistentInput(PersistentIStream & is, int) {
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<DipoleDensityAnalysis,DIPSY::DipoleAnalysisHandler>
  describeDIPSYDipoleDensityAnalysis("DIPSY::DipoleDensityAnalysis",
				   "DipoleDensityAnalysis.so");

void DipoleDensityAnalysis::Init() {


  static ClassDocumentation<DipoleDensityAnalysis> documentation
    ("There is no documentation for the DipoleDensityAnalysis class");

}

