// -*- C++ -*-
//
// PDFCuts.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PDFCuts class.
//

#include "PDFCuts.h"
#include "ThePEG/Cuts/Cuts.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PDFCuts.tcc"
#endif

using namespace ThePEG;

PDFCuts::
PDFCuts(const Cuts & kc, bool first, const SInterval & S,
	const SInterval & sHatIn) {
  SInterval sHat(max(sHatIn.lower(), kc.sHatMin()),
		 min(sHatIn.upper(), kc.sHatMax()));
  double xmin = sqrt(sHat.lower()/S.upper());
  double xmax = sqrt(sHat.upper()/S.lower());
  theScale = SInterval(kc.scaleMin(), kc.scaleMax());
  double x1max = min(xmax*exp( 2.*kc.yHatMax()), kc.x1Max());
  double x2max = min(xmax*exp(-2.*kc.yHatMin()), kc.x2Max());
  x1max = min(xmax,x1max);
  x2max = min(xmax,x2max);
  if ( first ) {
    double x1min = max(xmin*exp(2.*kc.yHatMin()), kc.x1Min());
    x1min = max(sqr(xmin),x1min);
    theSMax = x2max*S.upper();
    theL = Interval<double>(log(1.0/x1max), log(1.0/x1min));
  } 
  else {
    double x2min = max(xmin*exp(-2.*kc.yHatMax()), kc.x2Min());
    x2min = max(sqr(xmin),x2min);
    theSMax = x1max*S.upper();
    theL = Interval<double>(-log(x2max), -log(x2min));
  }
}

PDFCuts::
PDFCuts(const Cuts & kc, bool first, Energy maxEnergy) {
  SInterval sHat(kc.sHatMin(), kc.sHatMax());
  double xmin = sqrt(sHat.lower())/maxEnergy;
  double xmax = 1.0;
  theScale = SInterval(kc.scaleMin(), kc.scaleMax());
  double x1max = min(xmax*exp( 2.*kc.yHatMax()), kc.x1Max());
  double x2max = min(xmax*exp(-2.*kc.yHatMin()), kc.x2Max());
  x1max = min(xmax,x1max);
  x2max = min(xmax,x2max);
  if ( first ) {
    double x1min = max(xmin*exp(2.*kc.yHatMin()), kc.x1Min());
    x1min = max(sqr(xmin),x1min);
    theSMax = x2max*sqr(maxEnergy);
    theL = Interval<double>(-log(x1max), -log(x1min));
  } 
  else {
    double x2min = max(xmin*exp(-2.*kc.yHatMax()), kc.x2Min());
    x2min = max(sqr(xmin),x2min);
    theSMax = x1max*sqr(maxEnergy);
    theL = Interval<double>(-log(x2max), -log(x2min));
  }
}

