// -*- C++ -*-
//
// LWHFactory.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LWHFactory class.
//

#include "LWHFactory.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"

#ifndef LWH 
#define LWH ThePEGLWH
#endif

#include "AnalysisFactory.h"

using namespace ThePEG;

void LWHFactory::doinitrun() {
  analysisFactory(new LWH::AnalysisFactory);
  FactoryBase::doinitrun();
}

void LWHFactory::normalizeToXSec(tH1DPtr histogram, CrossSection unit) const {
  LWH::Histogram1D * h = dynamic_cast<LWH::Histogram1D *>(histogram);
  if ( h )
    h->normalize(h->sumAllBinHeights()*generator()->integratedXSec()/
		 (generator()->sumWeights()*unit));
}

void LWHFactory::normalizeToXSecFraction(tH1DPtr histogram) const {
  LWH::Histogram1D * h = dynamic_cast<LWH::Histogram1D *>(histogram);
  if ( h ) h->normalize(h->sumAllBinHeights()/generator()->sumWeights());
}

void LWHFactory::normalizeToUnity(tH1DPtr histogram) const {
  LWH::Histogram1D * h = dynamic_cast<LWH::Histogram1D *>(histogram);
  if ( h ) h->normalize(1.0);
}

void LWHFactory::normalizeToXSec(tH2DPtr histogram, CrossSection unit) const {
  LWH::Histogram2D * h = dynamic_cast<LWH::Histogram2D *>(histogram);
  if ( h )
    h->normalize(h->sumAllBinHeights()*generator()->integratedXSec()/
		 (generator()->sumWeights()*unit));
}

void LWHFactory::normalizeToXSecFraction(tH2DPtr histogram) const {
  LWH::Histogram2D * h = dynamic_cast<LWH::Histogram2D *>(histogram);
  if ( h ) h->normalize(h->sumAllBinHeights()/generator()->sumWeights());
}

void LWHFactory::normalizeToUnity(tH2DPtr histogram) const {
  LWH::Histogram2D * h = dynamic_cast<LWH::Histogram2D *>(histogram);
  if ( h ) h->normalize(1.0);
}

void LWHFactory::persistentOutput(PersistentOStream &) const {}

void LWHFactory::persistentInput(PersistentIStream &, int) {}

ClassDescription<LWHFactory> LWHFactory::initLWHFactory;
// Definition of the static class description member.

void LWHFactory::Init() {

  static ClassDocumentation<LWHFactory> documentation
    ("This class represents the Light-Weight Histogram package which "
     "implements the most rudimentary histogramming facilities according "
     "to the <a href=\"http://aida.freehep.org\">AIDA</a> interface "
     "specifications. Currently the only thing that is supported is "
     "simple, equally binned, one dimensional histograms. If you are "
     "using AnalysisHandlers which accesses other features in the AIDA "
     "interface you may end up with an ungraceful crash.");

}

