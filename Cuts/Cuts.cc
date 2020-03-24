// -*- C++ -*-
//
// Cuts.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Cuts class.
//

#include "Cuts.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/EventRecord/Collision.h"
#include "ThePEG/EventRecord/TmpTransform.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/HoldFlag.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/algorithm.h"

using namespace ThePEG;

Cuts::Cuts(Energy MhatMin)
  : theSMax(ZERO), theY(0), theCurrentSHat(-1.0*GeV2),
    theCurrentYHat(0), theMHatMin(MhatMin), theMHatMax(Constants::MaxEnergy),
    theYHatMin(-Constants::MaxRapidity), theYHatMax(Constants::MaxRapidity),
    theX1Min(0.0), theX1Max(1.0), theX2Min(0.0), theX2Max(1.0),
    theScaleMin(ZERO), theScaleMax(Constants::MaxEnergy2),
    theSubMirror(false), theCutWeight(1.0), theLastCutWeight(1.0) {}

Cuts::~Cuts() {}

IBPtr Cuts::clone() const {
  return new_ptr(*this);
}

IBPtr Cuts::fullclone() const {
  return new_ptr(*this);
}

void Cuts::doinitrun() {
  Interfaced::doinitrun();
  if ( Debug::level ) {
    describe();
    for_each(theOneCuts,   mem_fn(&OneCutBase::describe));
    for_each(theTwoCuts,   mem_fn(&TwoCutBase::describe));
    for_each(theMultiCuts, mem_fn(&MultiCutBase::describe));
  }
}

void Cuts::describe() const {
  CurrentGenerator::log() 
    << fullName() << ":\n"
    << "MHat  = " << theMHatMin/GeV << " .. " << theMHatMax/GeV << " GeV\n"
    << "Scale = " << theScaleMin/GeV2 << " .. " << theScaleMax/GeV2 << " GeV2\n"
    << "YHat  = " << theYHatMin << " .. " << theYHatMax << '\n'
    << "X1    = "   << theX1Min << " .. " << theX1Max << '\n'
    << "X2    = "   << theX2Min << " .. " << theX2Max << "\n\n";
  if ( theJetFinder )
    theJetFinder->describe();
}

void Cuts::initialize(Energy2 smax, double Y) {
  theSMax = smax;
  theMHatMax = min(theMHatMax, sqrt(smax));
  theY = Y;
  theSubMirror = false;
}

void Cuts::initEvent() {
  theCurrentSHat = -1.0*GeV2;
  theCurrentYHat = 0.0;
  theSubMirror = false;
}

bool Cuts::initSubProcess(Energy2 shat, double yhat, bool mirror) const {
  theCutWeight = 1.0;
  theLastCutWeight = 1.0;
  theSubMirror = mirror;
  theCurrentSHat = shat;
  theCurrentYHat = yhat;
  if ( shat <= sHatMin() || shat > sHatMax()*(1.0 + 1000.0*Constants::epsilon) ) {
    theCutWeight = 0.0;
    return false;
  }
  if ( yhat <= yHatMin() || yhat >= yHatMax() ) {
    theCutWeight = 0.0;
    return false;
  }
  double x1 = min(1.0, sqrt(shat/SMax())*exp(yhat));
  if ( x1 <= x1Min() || x1 > x1Max() ) {
    theCutWeight = 0.0;
    return false;
  }
  double x2 = min(1.0, sqrt(shat/SMax())*exp(-yhat));
  if ( x2 <= x2Min() || x2 > x2Max() ) {
    theCutWeight = 0.0;
    return false;
  }
  return true;
}

bool Cuts::passCuts(const tcPDVector & ptype, const vector<LorentzMomentum> & p,
		    tcPDPtr t1, tcPDPtr t2) const {
  if ( subMirror() ) {
    vector<LorentzMomentum> pmir = p;
    for ( int i = 0, N = pmir.size(); i < N; ++i ) pmir[i].setZ(-pmir[i].z());
    swap(t1,t2);
    HoldFlag<> nomir(theSubMirror, false);
    return passCuts(ptype, pmir, t1, t2);
  }

  bool pass = true;
  theCutWeight = 1.0;
  theLastCutWeight = 1.0;

  if ( jetFinder() ) {

    if ( ptype.size() > jetFinder()->minOutgoing() ) {


      vector<LorentzMomentum> jets;
      tcPDVector jettype;

      if ( !jetFinder()->restrictConsitutents() ) {
	jets = p;
	jettype = ptype;
      } else {
	tcPDVector::const_iterator pd = ptype.begin();
	vector<LorentzMomentum>::const_iterator pm = p.begin();
	for ( ; pd != ptype.end(); ++pd, ++pm ) {
	  if ( pm->rapidity() > jetFinder()->constituentRapidityRange().first &&
	       pm->rapidity() < jetFinder()->constituentRapidityRange().second ) {
	    jets.push_back(*pm);
	    jettype.push_back(*pd);
	  }
	}
      }

      if ( jetFinder()->cluster(jettype,jets,this,t1,t2) ){
	return passCuts(jettype,jets,t1,t2);
      }
    }

  }

  for ( int i = 0, N = p.size(); i < N; ++i )
    for ( int j = 0, M = theOneCuts.size(); j < M; ++j ) {
      pass &= theOneCuts[j]->passCuts(this, ptype[i], p[i]);
      theCutWeight *= theLastCutWeight;
      theLastCutWeight = 1.0;
      if ( !pass ) {
	theCutWeight = 0.0;
	return false;
      }
    }

  for ( int i1 = 0, N1 = p.size() - 1; i1 < N1; ++i1 )
    for ( int i2 = i1 + 1, N2 = p.size(); i2 < N2; ++i2 )
      for ( int j = 0, M = theTwoCuts.size(); j < M; ++j ) {
	pass &= theTwoCuts[j]->passCuts(this, ptype[i1], ptype[i2],
					p[i1], p[i2]);
	theCutWeight *= theLastCutWeight;
	theLastCutWeight = 1.0;
	if ( !pass ) {
	  theCutWeight = 0.0;
	  return false;
	}
      }

  for ( int j = 0, M = theMultiCuts.size(); j < M; ++j ) {
    pass &= theMultiCuts[j]->passCuts(this, ptype, p);
    theCutWeight *= theLastCutWeight;
    theLastCutWeight = 1.0;
    if ( !pass ) {
      theCutWeight = 0.0;
      return false;
    }
  }

  if ( t1 ) {
    LorentzMomentum p1(ZERO, ZERO, 0.5*sqrt(currentSHat()),
		       0.5*sqrt(currentSHat()));
    for ( int i = 0, N = p.size(); i < N; ++i )
      for ( int j = 0, M = theTwoCuts.size(); j < M; ++j ) {
	pass &= theTwoCuts[j]->passCuts(this, t1, ptype[i], p1, p[i],
					true, false);
	theCutWeight *= theLastCutWeight;
	theLastCutWeight = 1.0;
	if ( !pass ) {
	  theCutWeight = 0.0;
	  return false;
	}
      }
  }

  if ( t2 ) {
    LorentzMomentum p2(ZERO, ZERO,
		       -0.5*sqrt(currentSHat()), 0.5*sqrt(currentSHat()));
    for ( int i = 0, N = p.size(); i < N; ++i )
      for ( int j = 0, M = theTwoCuts.size(); j < M; ++j ) {
	pass &= theTwoCuts[j]->passCuts(this, ptype[i], t2, p[i], p2,
					false, true);
	theCutWeight *= theLastCutWeight;
	theLastCutWeight = 1.0;
	if ( !pass ) {
	  theCutWeight = 0.0;
	  return false;
	}
      }
  }

  return pass;

}

bool Cuts::passCuts(const tcPVector & p, tcPDPtr t1, tcPDPtr t2) const {
  tcPDVector ptype(p.size());
  vector<LorentzMomentum> mom(p.size());
  for ( int i = 0, N = p.size(); i < N; ++i ) {
    ptype[i] = p[i]->dataPtr();
    mom[i] = p[i]->momentum();
  }
  return passCuts(ptype, mom, t1, t2);
}

bool Cuts::passCuts(const SubProcess & sub) const {
  if ( !passCuts(tcPVector(sub.outgoing().begin(), sub.outgoing().end()),
		 sub.incoming().first->dataPtr(),
		 sub.incoming().second->dataPtr()) ) return false;
  return true;
}

bool Cuts::passCuts(const Collision & coll) const {
  tSubProPtr sub = coll.primarySubProcess();
  LorentzMomentum phat = sub->incoming().first->momentum() +
    sub->incoming().second->momentum();
  if ( !initSubProcess(phat.m2(), phat.rapidity()) ) return false;
  TmpTransform<tSubProPtr> tmp(sub, Utilities::getBoostToCM(sub->incoming()));
  if ( !passCuts(*sub) ) return false;
  return true;
}

Energy2 Cuts::minS(const tcPDVector & pv) const {
  Energy2 mins = ZERO;
  for ( int i = 0, N = theMultiCuts.size(); i < N; ++i )
    mins = max(mins, theMultiCuts[i]->minS(pv));
  return mins;
}

Energy2 Cuts::maxS(const tcPDVector & pv) const {
  Energy2 maxs = SMax();
  for ( int i = 0, N = theMultiCuts.size(); i < N; ++i )
    maxs = min(maxs, theMultiCuts[i]->maxS(pv));
  return maxs;
}

Energy2 Cuts::minSij(tcPDPtr pi, tcPDPtr pj) const {
  Energy2 mins = ZERO;
  for ( int i = 0, N = theTwoCuts.size(); i < N; ++i )
    mins = max(mins, theTwoCuts[i]->minSij(pi, pj));
  if ( mins > ZERO ) return mins;
  mins = sqr(pi->massMin() + pj->massMin());
  mins = max(mins, sqr(minKTClus(pi, pj))/4.0);
  mins = max(mins, minDurham(pi, pj)*currentSHat()/2.0);
  mins = max(mins, minKT(pi)*minKT(pj)*minDeltaR(pi, pj)/4.0);
  return mins;
}

Energy2 Cuts::minTij(tcPDPtr pi, tcPDPtr po) const {
  Energy2 mint = ZERO;
  for ( int i = 0, N = theTwoCuts.size(); i < N; ++i )
    mint = max(mint, theTwoCuts[i]->minTij(pi, po));
  if ( mint > ZERO ) return mint;
  mint = max(mint, sqr(minKT(po)));
  return mint;
}

double Cuts::minDeltaR(tcPDPtr pi, tcPDPtr pj) const {
  double mindr = 0.0;
  for ( int i = 0, N = theTwoCuts.size(); i < N; ++i )
    mindr = max(mindr, theTwoCuts[i]->minDeltaR(pi, pj));
  return mindr;
}

Energy Cuts::minKTClus(tcPDPtr pi, tcPDPtr pj) const {
  Energy minkt = ZERO;
  for ( int i = 0, N = theTwoCuts.size(); i < N; ++i )
    minkt = max(minkt, theTwoCuts[i]->minKTClus(pi, pj));
  return minkt;
}

double Cuts::minDurham(tcPDPtr pi, tcPDPtr pj) const {
  double y = 0.0;
  for ( int i = 0, N = theTwoCuts.size(); i < N; ++i )
    y = max(y, theTwoCuts[i]->minDurham(pi, pj));
  return y;
}

Energy Cuts::minKT(tcPDPtr p) const {
  Energy minkt = ZERO;
  for ( int i = 0, N = theOneCuts.size(); i < N; ++i )
    minkt = max(minkt, theOneCuts[i]->minKT(p));
  if ( minkt > ZERO ) return minkt;
  minkt = minKTClus(p, tcPDPtr());
  return minkt;
}

double Cuts::minEta(tcPDPtr p) const {
  double mineta = -Constants::MaxRapidity;
  for ( int i = 0, N = theOneCuts.size(); i < N; ++i )
    mineta = max(mineta, theOneCuts[i]->minEta(p));
  return mineta;
}

double Cuts::maxEta(tcPDPtr p) const {
  double maxeta = Constants::MaxRapidity;
  for ( int i = 0, N = theOneCuts.size(); i < N; ++i )
    maxeta = min(maxeta, theOneCuts[i]->maxEta(p));
  return maxeta;
}

double Cuts::minYStar(tcPDPtr p) const {
  if ( currentSHat() < ZERO ) return -Constants::MaxRapidity;
  if ( subMirror() ) {
    HoldFlag<>  nomir(theSubMirror, false);
    return -maxYStar(p);
  }
  double etamin = minEta(p);
  double ytot = Y() + currentYHat();
  if ( etamin > 0.0 ) {
    Energy minkt = minKT(p);
    Energy maxm = p->massMax();
    return asinh(minkt*sinh(etamin)/sqrt(sqr(minkt) + sqr(maxm))) - ytot;
  } else {
    return etamin - ytot;
  }
}

double Cuts::maxYStar(tcPDPtr p) const {
  if ( currentSHat() < ZERO ) return Constants::MaxRapidity;
  if ( subMirror() ) {
    HoldFlag<>  nomir(theSubMirror, false);
    return -minYStar(p);
  }
  double etamax = maxEta(p);
  double ytot = Y() + currentYHat();
  if ( etamax > 0.0 ) {
    return etamax - ytot;
  } else {
    Energy minkt = minKT(p);
    Energy maxm = p->massMax();
    return asinh(minkt*sinh(etamax)/sqrt(sqr(minkt) + sqr(maxm))) - ytot;
  }
}

double Cuts::minRapidityMax(tcPDPtr p) const {
  double minRapidityMax = -Constants::MaxRapidity;
  for ( int i = 0, N = theOneCuts.size(); i < N; ++i )
    minRapidityMax = max(minRapidityMax, theOneCuts[i]->minRapidityMax(p));
  return minRapidityMax;
}

double Cuts::maxRapidityMin(tcPDPtr p) const {
  double maxRapidityMin = Constants::MaxRapidity;
  for ( int i = 0, N = theOneCuts.size(); i < N; ++i )
    maxRapidityMin = min(maxRapidityMin, theOneCuts[i]->maxRapidityMin(p));
  return maxRapidityMin;
}

void Cuts::persistentOutput(PersistentOStream & os) const {
  os << ounit(theSMax, GeV2) << theY << ounit(theCurrentSHat, GeV2)
     << theCurrentYHat << ounit(theMHatMin, GeV) << ounit(theMHatMax, GeV)
     << theYHatMin << theYHatMax
     << theX1Min << theX1Max << theX2Min << theX2Max << ounit(theScaleMin, GeV2)
     << ounit(theScaleMax, GeV2) << theOneCuts << theTwoCuts << theMultiCuts
     << theJetFinder << theSubMirror
     << theCutWeight << theLastCutWeight << theFuzzyTheta;
}

void Cuts::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theSMax, GeV2) >> theY >> iunit(theCurrentSHat, GeV2)
     >> theCurrentYHat >> iunit(theMHatMin, GeV) >> iunit(theMHatMax, GeV)
     >> theYHatMin >> theYHatMax
     >> theX1Min >> theX1Max >> theX2Min >> theX2Max >> iunit(theScaleMin, GeV2)
     >> iunit(theScaleMax, GeV2) >> theOneCuts >> theTwoCuts >> theMultiCuts
     >> theJetFinder >> theSubMirror
     >> theCutWeight >> theLastCutWeight >> theFuzzyTheta;
}

ClassDescription<Cuts> Cuts::initCuts;
// Definition of the static class description member.

Energy Cuts::maxMHatMin() const {
  return theMHatMax;
}

Energy Cuts::minMHatMax() const {
  return theMHatMin;
}

Energy2 Cuts::maxScaleMin() const {
  return theScaleMax;
}

Energy2 Cuts::minScaleMax() const {
  return theScaleMin;
}

double Cuts::maxYHatMin() const {
  return theYHatMax;
}

double Cuts::minYHatMax() const {
  return theYHatMin;
}

double Cuts::maxX1Min() const {
  return theX1Max;
}

double Cuts::minX1Max() const {
  return theX1Min;
}

double Cuts::maxX2Min() const {
  return theX2Max;
}

double Cuts::minX2Max() const {
  return theX2Min;
}

void Cuts::Init() {

  typedef double (ThePEG::Cuts::*IGFN)() const;
  typedef void (ThePEG::Cuts::*ISFN)(double);

  static ClassDocumentation<Cuts> documentation
    ("Cuts is a class for implementing kinematical cuts in ThePEG. The "
     "class itself only implements cuts on the total momentum of the hard "
     "sub-process, implemented as minimum and maximum values of \\f$x_1\\f$ "
     "and \\f$x_2\\f$ (or \\f$\\hat{s}\\f$ and \\f$\\hat{y}\\f$. Further cuts "
     "can be implemented either by inheriting from this base class, in which "
     "the virtual cut() function should be overridden, or by assigning "
     "objects of class OneCutBase, TwoCutBase and MultiCutBase defining "
     "cuts on single particles, pairs of particles and groups of "
     "particles respectively.");

  static Parameter<Cuts,Energy> interfaceMHatMin
    ("MHatMin",
     "The minimum allowed value of \\f$\\sqrt{\\hat{s}}\\f$.",
     &Cuts::theMHatMin, GeV, 2.0*GeV, ZERO, Constants::MaxEnergy,
     true, false, Interface::limited,
     0, 0, 0, &Cuts::maxMHatMin, 0);
  interfaceMHatMin.setHasDefault(false);

  static Parameter<Cuts,Energy> interfaceMHatMax
    ("MHatMax",
     "The maximum allowed value of \\f$\\sqrt{\\hat{s}}\\f$.",
     &Cuts::theMHatMax, GeV, 100.0*GeV, ZERO, ZERO,
     true, false, Interface::lowerlim,
     0, 0,
     &Cuts::minMHatMax, 0, 0);
  interfaceMHatMax.setHasDefault(false);

  static Parameter<Cuts,Energy2> interfaceScaleMin
    ("ScaleMin",
     "The minimum allowed value of the scale to be used in PDFs and "
     "coupling constants.",
     &Cuts::theScaleMin, GeV2, ZERO, ZERO, Constants::MaxEnergy2,
     true, false, Interface::limited,
     0, 0, 0, &Cuts::maxScaleMin, 0);
  interfaceScaleMin.setHasDefault(false);

  static Parameter<Cuts,Energy2> interfaceScaleMax
    ("ScaleMax",
     "The maximum allowed value of the scale to be used in PDFs and "
     "coupling constants.",
     &Cuts::theScaleMax, GeV2, 10000.0*GeV2, ZERO, ZERO,
     true, false, Interface::lowerlim,
     0, 0,
     &Cuts::minScaleMax, 0, 0);
  interfaceScaleMax.setHasDefault(false);

  static Parameter<Cuts,double> interfaceYHatMin
    ("YHatMin",
     "The minimum value of the rapidity of the hard sub-process "
     "(wrt. the rest system of the colliding particles).",
     &Cuts::theYHatMin, -10.0, 0.0, Constants::MaxRapidity,
     true, false, Interface::upperlim,
     (ISFN)0, (IGFN)0, (IGFN)0, &Cuts::maxYHatMin, (IGFN)0);
  interfaceYHatMin.setHasDefault(false);

  static Parameter<Cuts,double> interfaceYHatMax
    ("YHatMax",
     "The maximum value of the rapidity of the hard sub-process "
     "(wrt. the rest system of the colliding particles).",
     &Cuts::theYHatMax, 10.0, -Constants::MaxRapidity, 0.0,
     true, false, Interface::lowerlim,
     (ISFN)0, (IGFN)0, &Cuts::minYHatMax, (IGFN)0, (IGFN)0);
  interfaceYHatMax.setHasDefault(false);

  static Parameter<Cuts,double> interfaceX1Min
    ("X1Min",
     "The minimum value of the positive light-cone fraction of the hard "
     "sub-process.",
     &Cuts::theX1Min, 0.0, 0.0, 1.0,
     true, false, Interface::limited,
     (ISFN)0, (IGFN)0, (IGFN)0, &Cuts::maxX1Min, (IGFN)0);
  interfaceX1Min.setHasDefault(false);

  static Parameter<Cuts,double> interfaceX1Max
    ("X1Max",
     "The maximum value of the positive light-cone fraction of the hard "
     "sub-process.",
     &Cuts::theX1Max, 0.0, 0.0, 1.0,
     true, false, Interface::limited,
     (ISFN)0, (IGFN)0, &Cuts::minX1Max, (IGFN)0, (IGFN)0);
  interfaceX1Max.setHasDefault(false);

  static Parameter<Cuts,double> interfaceX2Min
    ("X2Min",
     "The minimum value of the negative light-cone fraction of the hard "
     "sub-process.",
     &Cuts::theX2Min, 0.0, 0.0, 1.0,
     true, false, Interface::limited,
     (ISFN)0, (IGFN)0, (IGFN)0, &Cuts::maxX2Min, (IGFN)0);
  interfaceX2Min.setHasDefault(false);

  static Parameter<Cuts,double> interfaceX2Max
    ("X2Max",
     "The maximum value of the negative light-cone fraction of the hard "
     "sub-process.",
     &Cuts::theX2Max, 0.0, 0.0, 1.0,
     true, false, Interface::limited,
     (ISFN)0, (IGFN)0, &Cuts::minX2Max, (IGFN)0, (IGFN)0);
  interfaceX2Max.setHasDefault(false);

  static RefVector<Cuts,OneCutBase> interfaceOneCuts
    ("OneCuts",
     "The objects defining cuts on single outgoing partons from the "
     "hard sub-process.",
     &Cuts::theOneCuts, -1, true, false, true, false, false);

  static RefVector<Cuts,TwoCutBase> interfaceTwoCuts
    ("TwoCuts",
     "The objects defining cuts on pairs of particles in the "
     "hard sub-process.",
     &Cuts::theTwoCuts, -1, true, false, true, false, false);

  static RefVector<Cuts,MultiCutBase> interfaceMultiCuts
    ("MultiCuts",
     "The objects defining cuts on sets of outgoing particles from the "
     "hard sub-process.",
     &Cuts::theMultiCuts, -1, true, false, true, false, false);

  static Reference<Cuts,JetFinder> interfaceJetFinder
    ("JetFinder",
     "Set a JetFinder object used to define cuts on the"
     "level of reconstructed jets as needed for higher order corrections.",
     &Cuts::theJetFinder, false, false, true, true, false);

  static Reference<Cuts,FuzzyTheta> interfaceFuzzy
    ("Fuzzy",
     "The fuzziness to be applied to cuts (may not be supported by all cut objects).",
     &Cuts::theFuzzyTheta, false, false, true, true, false);


  interfaceX1Min.rank(10);
  interfaceX1Max.rank(9);
  interfaceX2Min.rank(8);
  interfaceX2Max.rank(7);
  interfaceMHatMin.rank(6);
  interfaceMHatMax.rank(5);
  interfaceYHatMin.rank(4);
  interfaceYHatMax.rank(3);
  interfaceOneCuts.rank(2);
  interfaceTwoCuts.rank(1);

}

double Cuts::yHatMin() const {
  return theX1Min > 0.0 && theX2Max > 0.0?
    max(theYHatMin, 0.5*log(theX1Min/theX2Max)): theYHatMin;
}

double Cuts::yHatMax() const {
  return theX1Max > 0.0 && theX2Min > 0.0?
    min(theYHatMax, 0.5*log(theX1Max/theX2Min)): theYHatMax;
}

bool Cuts::yHat(double y) const {
  return y > yHatMin() && y < yHatMax();
}

double Cuts::x1Min() const {
  return max(theX1Min, (theMHatMin/sqrt(SMax()))*exp(theYHatMin));
}

double Cuts::x1Max() const {
  return min(theX1Max, (theMHatMax/sqrt(SMax()))*exp(theYHatMax));
}

bool Cuts::x1(double x) const {
  return x > x1Min() && x <= x1Max();
}

double Cuts::x2Min() const {
  return max(theX2Min, (theMHatMin/sqrt(SMax()))/exp(theYHatMax));
}

double Cuts::x2Max() const {
  return min(theX2Max, (theMHatMax/sqrt(SMax()))/exp(theYHatMin));
}

bool Cuts::x2(double x) const {
  return x > x2Min() && x <= x2Max();
}

template <typename T>
vector<typename Ptr<T>::transient_const_pointer>
Cuts::oneCutObjects() const {
  typedef typename Ptr<T>::transient_const_pointer tcPtr;
  vector<tcPtr> ret;
  for ( int i = 0, N = theOneCuts.size(); i < N; ++i )
    if ( dynamic_ptr_cast<tcPtr>(theOneCuts[i]) )
      ret.push_back(dynamic_ptr_cast<tcPtr>(theOneCuts[i]));
  return ret;
}

template <typename T>
vector<typename Ptr<T>::transient_const_pointer>
Cuts::twoCutObjects() const {
  typedef typename Ptr<T>::transient_const_pointer tcPtr;
  vector<tcPtr> ret;
  for ( int i = 0, N = theTwoCuts.size(); i < N; ++i )
    if ( dynamic_ptr_cast<tcPtr>(theTwoCuts[i]) )
      ret.push_back(dynamic_ptr_cast<tcPtr>(theTwoCuts[i]));
  return ret;
}

template <typename T>
vector<typename Ptr<T>::transient_const_pointer>
Cuts::multiCutObjects() const {
  typedef typename Ptr<T>::transient_const_pointer tcPtr;
  vector<tcPtr> ret;
  for ( int i = 0, N = theMultiCuts.size(); i < N; ++i )
    if ( dynamic_ptr_cast<tcPtr>(theMultiCuts[i]) )
      ret.push_back(dynamic_ptr_cast<tcPtr>(theMultiCuts[i]));
  return ret;
}
