// -*- C++ -*-
//
// PartonBin.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PartonBin class.
//

#include "PartonBin.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDF/PDFBase.h"
#include "ThePEG/PDF/NoPDF.h"
#include "ThePEG/PDF/RemnantHandler.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/Maths.h"

using namespace ThePEG;

PartonBin::PartonBin() : thePDFDim(), theRemDim() {}

PartonBin::
PartonBin(tcPDPtr p, tPBPtr inc, tcPDPtr pi,
	  tcPDFPtr pdf, const PDFCuts & newCuts)
  : theParticle(p), theIncomingBin(inc), theParton(pi), thePDF(pdf),
    thePDFDim(0), theRemDim(0), theCuts(newCuts) {
  if ( pdf ) theRemnantHandler = pdf->remnantHandler();
}

PartonBin::~PartonBin() {}

int PartonBin::nDim(bool doscale) {
  if ( !incoming() ) return 0;
  if ( dynamic_ptr_cast<Ptr<NoPDF>::tcp>(pdf()) ) thePDFDim = 0;
  else if ( doscale ) thePDFDim = 2;
  else thePDFDim = 1;
  theRemDim = remnantHandler()->nDim(*this, !doscale);
  return pdfDim() + remDim() + incoming()->nDim(true);
}

tPBPtr PartonBin::getFirst() {
  return incoming()? incoming()->getFirst(): tPBPtr(this);
}

void PartonBin::persistentOutput(PersistentOStream & os) const {
  os << theParticle << theIncomingBin << theOutgoing << theParton << thePDF
     << theRemnantHandler << thePDFDim << theRemDim
     << cuts().lMin() << cuts().lMax() << ounit(cuts().scaleMin(), GeV2)
     << ounit(cuts().scaleMax(), GeV2) << ounit(cuts().sMax(), GeV2);
}

void PartonBin::persistentInput(PersistentIStream & is, int) {
  double lmin = 0.0;
  double lmax = 0.0;
  Energy2 scmin = ZERO;
  Energy2 scmax = ZERO;
  Energy2 smax = ZERO;
  is >> theParticle >> theIncomingBin >> theOutgoing >> theParton >> thePDF
     >> theRemnantHandler >> thePDFDim >> theRemDim
     >> lmin >> lmax >> iunit(scmin, GeV2) >> iunit(scmax, GeV2)
     >> iunit(smax, GeV2);
  theCuts = PDFCuts(Interval<double>(lmin, lmax),
		    SInterval(scmin, scmax), smax);
}

ClassDescription<PartonBin> PartonBin::initPartonBin;

void PartonBin::Init() {}
