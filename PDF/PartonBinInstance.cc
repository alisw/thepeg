// -*- C++ -*-
//
// PartonBinInstance.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PartonBinInstance class.
//

#include "PartonBinInstance.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

PartonBinInstance::PartonBinInstance()
  : theJacobian(1.0), theXi(-1.0), theEps(-1.0), theLi(-1.0), theX(-1.0),
    theL(-1.0), theScale(ZERO), theRemnantWeight(0.0) {}

PartonBinInstance::PartonBinInstance(const PartonBinInstance & x)
  : Base(x),
    theBin(x.theBin), theBins(x.theBins), theIncoming(x.theIncoming),
    theJacobian(x.theJacobian), theParticle(x.theParticle),
    theParton(x.theParton), thePartons(x.thePartons), theXi(x.theXi),
    theEps(x.theEps), theLi(x.theLi), theX(x.theX), theL(x.theL),
    theScale(x.theScale), theKT(x.theKT), theRemnantWeight(x.theRemnantWeight),
    theRemnants(x.theRemnants), theRemInfo(x.theRemInfo) {}

PartonBinInstance::PartonBinInstance(tcPBPtr pb, tPBIPtr pbi)
  : theBin(pb), theJacobian(1.0), theXi(-1.0), theEps(-1.0), theLi(-1.0),
    theX(-1.0), theL(-1.0), theScale(ZERO), theRemnantWeight(0.0) {
  if ( pbi ) theIncoming = pbi;
  else if ( bin()->incoming() )
    theIncoming = new_ptr(PartonBinInstance(bin()->incoming()));
}

PartonBinInstance::PartonBinInstance(tPPtr part, tcPBPtr pb, Energy2 scale)
  : theBin(pb), theJacobian(1.0), theParton(part), theXi(1.0), theEps(0.0),
    theLi(0.0), theX(1.0), theL(0.0), theScale(scale),
    theRemnantWeight(1.0) {
  if ( !pb->incoming() || part->parents().empty() ) return;
  particle(parton()->parents()[0]);
  Energy2 P2 = max(-particle()->momentum().m2(), ZERO);
  theXi = parton()->momentum().dirPlus()/particle()->momentum().dirPlus();
  theLi = -log(xi());
  theIncoming = new_ptr(PartonBinInstance(particle(), pb->incoming(), P2));
  theX = xi()*incoming()->x();
  theL = li() + incoming()->li();
  theEps =  Math::exp1m(-li());
}

PartonBinInstance::PartonBinInstance(tPPtr Part, tPPtr part, tcPBPtr pb, Energy2 scale)
  : theBin(pb), theJacobian(1.0), theParton(part), theXi(1.0), theEps(0.0),
    theLi(0.0), theX(1.0), theL(0.0), theScale(scale),
    theRemnantWeight(1.0) {
  if ( !pb->incoming() ) return;
  particle(Part);
  Energy2 P2 = max(-particle()->momentum().m2(), ZERO);
  theXi = parton()->momentum().dirPlus()/particle()->momentum().dirPlus();
  theLi = -log(xi());
  theIncoming = new_ptr(PartonBinInstance(particle(), pb->incoming(), P2));
  theX = xi()*incoming()->x();
  theL = li() + incoming()->li();
  theEps =  Math::exp1m(-li());
}

PartonBinInstance::~PartonBinInstance() {}

tPBIPtr PartonBinInstance::getFirst() {
  return incoming()? incoming()->getFirst(): tPBIPtr(this);
}


// TAKE AWAY ?
void PartonBinInstance::reset(double lx, Energy2 Q2) {
  l(lx);
  li(lx);
  scale(Q2);
  particle(tPPtr());
  parton(tPPtr());
  theRemnants.clear();
  thePartons.clear();
  remnantWeight(1.0);
}

void PartonBinInstance::prepare() {
  li(-1.0);
  l(0.0);
  if ( !incoming() ) return;
  //  li(-1.0);
  //  l(-1.0);
  reset(-1.0, ZERO);
  incoming()->prepare();
}

bool PartonBinInstance::hasPoleIn1() const {
  return ( !incoming() || incoming()->hasPoleIn1()) &&
    (!pdf() || pdf()->hasPoleIn1(particleData(), partonData()) );
}

// TAKE AWAY ?
void PartonBinInstance::generate(const double * r) {
  scale(ZERO);
  if ( !incoming() ) return;
  if ( li() >= 0 ) return;
  li(0.0);
  jacobian(1.0);
  remnantWeight(1.0);
  if ( bin()->pdfDim() )
    li(pdf()->flattenL(particleData(), partonData(), bin()->cuts(),
		       *r++, theJacobian));
  if ( bin()->pdfDim() > 1 )
    scale(pdf()->flattenScale(particleData(), partonData(), bin()->cuts(),
			      li(), *r++, theJacobian)
	  *bin()->cuts().scaleMaxL(li()));
  incoming()->generate(r);
  l(li() + incoming()->l());
}

void PartonBinInstance::persistentOutput(PersistentOStream & os) const {
  os << theBin << theBins << theIncoming << theJacobian << theParticle
     << theParton << thePartons << theXi << theEps << theLi << theX << theL
     << ounit(theScale, GeV2) << ounit(theKT, GeV) << theRemnantWeight
     << theRemnants;
}

void PartonBinInstance::persistentInput(PersistentIStream & is, int) {
  is >> theBin >> theBins >> theIncoming >> theJacobian >> theParticle
     >> theParton >> thePartons >> theXi >> theEps >> theLi >> theX >> theL
     >> iunit(theScale, GeV2) >> iunit(theKT, GeV) >> theRemnantWeight
     >> theRemnants;
  theRemInfo = RemIPtr();
}

ClassDescription<PartonBinInstance> PartonBinInstance::initPartonBinInstance;
// Definition of the static class description member.

void PartonBinInstance::Init() {

}

