// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RemnantParton class.
//

#include "RemnantParton.h"
#include "Emission.h"
#include "DipoleState.h"
#include "AriadneHandler.h"
#include "Models/RemnantModel.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/EventRecord/Particle.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/DebugItem.h"

using namespace Ariadne5;

RemnantParton::RemnantParton()
  : Parton(true), theModel(&Current<AriadneHandler>()->remnantModel()),
    theMu(0.0*GeV), theAlpha(0.0), theBeta(0.0), theX(1.0),
    theMuF2(ZERO), theY(0.0), theDirection(0) {}


ClonePtr RemnantParton::clone() const {
  return new_ptr(*this);
}

void RemnantParton::setupHard(int indir) {
  dir(indir);
  setBoost();
  theAlpha = Current<AriadneHandler>()->hardAlpha();
  theMu = sqrt(abs(parentMomentum().m2()));
  theBeta = Current<AriadneHandler>()->beta();
  Lorentz5Momentum p = state()->incoming().first->momentum();
  Lorentz5Momentum q = state()->incoming().second->momentum();
  if ( dir() > 0 ) swap(p, q);
  theY = p.dot(parentMomentum())/p.dot(q);
  theX = 0.0;
}

void RemnantParton::setup(const PartonBinInstance & pb, int indir) {
  dir(indir);
  setupParent(pb.particle());
  setupExtracted(pb.parton());
  orig(pb.parton());
  setupPDF(PDF(&pb));
  x(pb.x());
  muF2(pb.scale());
  // This is a soft remnant, so we assign the anti-particle type of
  // the extracted parton. Also we setup the momentum adn colour
  // lines.
  setData();
  momentum() = parentMomentum() - extractedMomentum();
  if ( data().hasAntiColour() && parentData().hasAntiColour() )
    origICol(parent()->antiColourLine());
  else if ( data().hasAntiColour() )
    origICol(originalExtracted()->colourLine());
  else
    origICol(tColinePtr());
  if ( data().hasColour() && parentData().hasColour() )
    origOCol(parent()->colourLine());
  else if ( data().hasColour() )
    origOCol(originalExtracted()->antiColourLine()); 
  else
    origOCol(tColinePtr());
  setBoost();
  theAlpha = Current<AriadneHandler>()->softAlpha();
  theMu = Current<AriadneHandler>()->softMu();
  theBeta = Current<AriadneHandler>()->beta();
}

tPPtr RemnantParton::produceParticle(const LorentzRotation & r) {
  if ( hard() ) return Parton::produceParticle(r);
  if ( untouched() ) return theExtracted = originalExtracted();
  theExtracted = extractedData().produceParticle(extractedMomentum());
  theExtracted->transform(r);
  return theExtracted;
}

const LorentzRotation & RemnantParton::setBoost() {
  Lorentz5Momentum po = dir() > 0?
    state()->incoming().second->momentum():
    state()->incoming().first->momentum();
  //  Lorentz5Momentum tot = parentMomentum() + po;
  theBoost = theInvBoost =
    Utilities::getBoostToCM(make_pair(parentMomentum(), po));
  theInvBoost.invert();
  return theBoost;
}

void RemnantParton::x(double xin) {
  theX = xin;
  if ( xin <= 0.0 ) theX = 1.0 - (getBoost()*momentum()).plus()/
		                 (getBoost()*parentMomentum()).plus();
}

Energy RemnantParton::effectiveMass(tcPDPtr extracted) const {
  return model().effectiveMass(this, extracted);
}

double RemnantParton::xfratio(tcPDPtr newp, Energy2 scale, double z) const {
  static DebugItem debugpdf("Ariadne5::RemnantParton::PDF");
  static double minpdf = 1.0e-5;
  if ( !pdf().pdf() ) return 0.0;
  double den = xfx(scale);
  if ( debugpdf && den > 0.0 && den < minpdf ) {
    cerr << "> Ariadne5::RemnantParton::PDF" << endl
	 << "  Found new small pdf = " << den << " at x = " << x()
	 << ", Q = " << sqrt(scale)/GeV
	 << ", q = " << extractedData().id() << endl;
    minpdf = den;
  }
  if ( den <= 0.0 ) return -1.0;
  if ( z < x() ) return -1.0;
  return xfx(newp, scale, x()/z)/den;
}

double RemnantParton::recoilWeight(const LorentzMomentum & ph,
				   const LorentzMomentum & pr,
				   tcPDPtr extracted) const {
  if ( !extracted ) extracted = theExtractedData;
  return model().recoilWeight(getBoost()*ph, getBoost()*pr, this, extracted);
}

double RemnantParton::recoilWeight() const {
  return recoilWeight(state()->hadronicMomentum(), momentum());
}

double RemnantParton::
reweightFS(Energy rho, LorentzMomentum pem,
	   LorentzMomentum prem, tcPDPtr extracted) const {
  return model().reweightFS(this, rho, pem, prem, extracted);
}

double RemnantParton::softSuppression(Energy rho, double xplus) const {
  return model().softSuppression(rho, xplus, x(), mu(), alpha(), beta());
}
  

LorentzRotation
RemnantParton::getHardTransform(const LorentzMomentum & phold,
				const LorentzMomentum & phnew) const {
  LorentzMomentum k1 = parentMomentum() - momentum();
  LorentzMomentum k = phold - k1;
  LorentzMomentum k3 = k1 + phnew - phold;
  return Utilities::getBoostFromCM(make_pair(k3, k))*
         Utilities::getBoostToCM(make_pair(k1, k));
}

void RemnantParton::setMomentum(const Lorentz5Momentum & p) {
  state()->touchHadronicState();
  momentum() = p;
  setExtracted();
}

void RemnantParton::fillReferences(CloneSet & cset) const {
  Parton::fillReferences(cset);
}

void RemnantParton::rebind(const TranslationMap & trans) {
  Parton::rebind(trans);
}

void RemnantParton::persistentOutput(PersistentOStream & os) const {
  os << ounit(theMu, GeV) << theAlpha << theBeta << theParent << theParentData
     << ounit(theParentMomentum, GeV) << theExtracted << theExtractedData
     << ounit(theExtractedMomentum, GeV) << theX << ounit(theMuF2, GeV2)
     << theY << theDirection << thePDF.pdf() << thePDF.particle()
     << theModel;
}

void RemnantParton::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theMu, GeV) >> theAlpha >> theBeta >> theParent >> theParentData
     >> iunit(theParentMomentum, GeV) >> theExtracted >> theExtractedData
     >> iunit(theExtractedMomentum, GeV) >> theX >> iunit(theMuF2, GeV2)
     >> theY >> theDirection ;
  tcPDFPtr f;
  tcPDPtr pd;
  is >> f >> pd >> theModel;
  thePDF = PDF(f, pd);
  setBoost();
}

DescribeClass<RemnantParton,Parton>
describeAriadne5RemnantParton("Ariadne5::RemnantParton", "libAriadne5.so");

void RemnantParton::Init() {}

void RemnantParton::debugme() const {
  Parton::debugme();
}

