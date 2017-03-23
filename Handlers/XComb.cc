// -*- C++ -*-
//
// XComb.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the XComb class.
//

#include "XComb.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/SubProcessHandler.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Utilities/Maths.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/Vectors/LorentzRotation.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/MatrixElement/ColourLines.h"
#include "ThePEG/Handlers/LuminosityFunction.h"
#include "ThePEG/Handlers/CascadeHandler.h"
#include "ThePEG/Repository/EventGenerator.h"

using namespace ThePEG;

XComb::XComb()
  : theLastS(Energy2()), theLastSHat(Energy2()), theLastY(0.0),
    theLastP1P2(make_pair(1.0, 1.0)), theLastL1L2(make_pair(1.0, 1.0)),
    theLastX1X2(make_pair(1.0, 1.0)), theLastE1E2(make_pair(0.0, 0.0)),
    theLastScale(ZERO), theLastCentralScale(ZERO), theLastShowerScale(ZERO),
    theLastAlphaS(-1.0), theLastAlphaEM(-1.0), theMaxEnergy(ZERO) {}

XComb::
XComb(Energy newMaxEnergy, const cPDPair & inc, tEHPtr newEventHandler,
      tPExtrPtr newExtractor, tCascHdlPtr newCKKW, const PBPair & newPartonBins,
      tCutsPtr newCuts)
  : theEventHandler(newEventHandler),
    thePartonExtractor(newExtractor), theCKKW(newCKKW), theCuts(newCuts),
    theParticles(inc), thePartonBins(newPartonBins), theLastS(Energy2()),
    theLastSHat(Energy2()), theLastY(0.0), theLastP1P2(make_pair(1.0, 1.0)),
    theLastL1L2(make_pair(1.0, 1.0)), theLastX1X2(make_pair(1.0, 1.0)),
    theLastE1E2(make_pair(0.0, 0.0)), theLastScale(ZERO), theLastCentralScale(ZERO),
    theLastShowerScale(ZERO), theLastAlphaS(-1.0), theLastAlphaEM(-1.0),
  theMaxEnergy(newMaxEnergy) {
  thePartons = cPDPair(partonBins().first->parton(),
		       partonBins().second->parton());
  thePartonBinInstances.first =
    new_ptr(PartonBinInstance(partonBins().first));
  thePartonBinInstances.second =
    new_ptr(PartonBinInstance(partonBins().second));
  theParticleBins.first = thePartonBins.first->getFirst();
  theParticleBins.second = thePartonBins.second->getFirst();
}

XComb::~XComb() {}

void XComb::clean() {
  theLastParticles = PPair();
  theLastPartons = PPair();
  theLastS = theLastSHat = theLastScale = 
    theLastCentralScale = theLastShowerScale = ZERO;
  theLastAlphaS = theLastAlphaEM = -1.0;
  theLastY = 0.0;
  theLastP1P2 = theLastL1L2 = theLastX1X2 = theLastE1E2 = DPair(0.0, 0.0);
  theSub = SubProPtr();
  thePartonBinInstances = PBIPair();
  thePartonBinInstanceMap.clear();
}

void XComb::prepare(const PPair & inc) {
  clean();
  createPartonBinInstances();
  theLastParticles = inc;
  pExtractor()->select(this);
  pExtractor()->prepare(partonBinInstances());
}

void XComb::subProcess(tSubProPtr sp) {
  theSub = sp;
}

void XComb::setPartonBinInfo() {
  partonBinInstances().first->getFirst()->parton(lastParticles().first);
  partonBinInstances().second->getFirst()->parton(lastParticles().second);
}

void XComb::createPartonBinInstances() {
  thePartonBinInstances.first =
    new_ptr(PartonBinInstance(partonBins().first));
  thePartonBinInstances.second =
    new_ptr(PartonBinInstance(partonBins().second));
}

void XComb::setPartonBinInstances(PBIPair pbip, Energy2 scale) {
  clean();
  thePartonBinInstances = pbip;
  theLastParticles = PPair(pbip.first->getFirst()->parton(),
			   pbip.second->getFirst()->parton());
  theLastPartons = PPair(pbip.first->parton(),
			 pbip.second->parton());
  lastS((lastParticles().first->momentum() +
	 lastParticles().second->momentum()).m2());
  lastSHat((lastPartons().first->momentum() +
	    lastPartons().second->momentum()).m2());
  lastP1P2(make_pair(0.0, 0.0));
  lastX1X2(make_pair(lastPartons().first->momentum().plus()/
		     lastParticles().first->momentum().plus(),
		     lastPartons().second->momentum().minus()/
		     lastParticles().second->momentum().minus()));
  lastY(log(lastX1()/lastX2())*0.5);
  lastScale(scale);
}

void XComb::lastL1L2(pair<double,double> ll) {
  theLastL1L2 = ll;
  theLastX1X2 = make_pair(exp(-ll.first), exp(-ll.second));
  theLastE1E2 = make_pair(Math::exp1m(-ll.first), Math::exp1m(-ll.second));
}

void XComb::lastX1X2(pair<double,double> xx) {
  theLastX1X2 = xx;
  theLastL1L2 = make_pair(-log(xx.first), -log(xx.second));
  theLastE1E2 = make_pair(1.0 - xx.first, 1.0 - xx.second);
}

void XComb::lastE1E2(pair<double,double> ee) {
  theLastE1E2= ee;
  theLastL1L2 = make_pair(-Math::log1m(ee.first), -Math::log1m(ee.second));
  theLastX1X2 = make_pair(1.0 - ee.first, 1.0 - ee.second);
}

tPBIPtr XComb::partonBinInstance(tcPPtr p) const {
  return pExtractor()->partonBinInstance(p);
}

void XComb::Init() {}

void XComb::persistentOutput(PersistentOStream & os) const {
  os << theEventHandler << thePartonExtractor << theCKKW
     << theCuts << theParticles << thePartons << thePartonBins
     << theParticleBins << thePartonBinInstances
     << theLastParticles << theLastPartons
     << ounit(theLastS, GeV2) << ounit(theLastSHat, GeV2) << theLastY
     << theLastP1P2 << theLastL1L2 << theLastX1X2 << theLastE1E2
     << ounit(theLastScale, GeV2) << ounit(theLastCentralScale, GeV2) 
      << ounit(theLastShowerScale, GeV2) << theLastAlphaS << theLastAlphaEM
     << ounit(theMaxEnergy, GeV) << theMEInfo << theSub;
}

void XComb::persistentInput(PersistentIStream & is, int) {
  is >> theEventHandler >> thePartonExtractor >> theCKKW
     >> theCuts >> theParticles >> thePartons >> thePartonBins
     >> theParticleBins >> thePartonBinInstances
     >> theLastParticles >> theLastPartons
     >> iunit(theLastS, GeV2) >> iunit(theLastSHat, GeV2) >> theLastY
     >> theLastP1P2 >> theLastL1L2 >> theLastX1X2 >> theLastE1E2
     >> iunit(theLastScale, GeV2) >> iunit(theLastCentralScale, GeV2) 
      >> iunit(theLastShowerScale, GeV2) >> theLastAlphaS >> theLastAlphaEM
     >> iunit(theMaxEnergy, GeV) >> theMEInfo >> theSub;
}

ClassDescription<XComb> XComb::initXComb;

