// -*- C++ -*-
//
// StdXCombGroup.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2007 Leif Lonnblad
// Copyright (C) 2009-2010 Simon Platzer
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StdXCombGroup class.
//

#include "StdXCombGroup.h"
#include "ThePEG/MatrixElement/MEGroup.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
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
#include "ThePEG/EventRecord/SubProcessGroup.h"
#include "ThePEG/Vectors/LorentzRotation.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/MatrixElement/ColourLines.h"
#include "ThePEG/Handlers/LuminosityFunction.h"
#include "ThePEG/Handlers/CascadeHandler.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/TmpTransform.h"

using namespace ThePEG;

StdXCombGroup::StdXCombGroup(Energy newMaxEnergy, const cPDPair & inc,
			     tEHPtr newEventHandler,tSubHdlPtr newSubProcessHandler,
			     tPExtrPtr newExtractor,	tCascHdlPtr newCKKW,
			     const PBPair & newPartonBins, tCutsPtr newCuts, tMEGroupPtr newME,
			     const DiagramVector & newDiagrams, bool mir, tStdXCombPtr newHead)
  : StandardXComb(newMaxEnergy,inc,newEventHandler,newSubProcessHandler,
		  newExtractor, newCKKW, newPartonBins, newCuts,
		  newME, newDiagrams, mir, newHead), 
    theMEGroup(newME), theDependent(), theLastHeadCrossSection(ZERO) {}

StdXCombGroup::StdXCombGroup()
  : StandardXComb(), theDependent() {}

void StdXCombGroup::build(const PartonPairVec& allPBins) {
  for ( MEVector::const_iterator me = theMEGroup->dependent().begin();
	me != theMEGroup->dependent().end(); ++me ) {
    vector<StdXCombPtr> dep = 
      theMEGroup->makeDependentXCombs(this,diagrams().front()->partons(),*me,allPBins);
    copy(dep.begin(),dep.end(),back_inserter(theDependent));
  }
}

StdXCombGroup::~StdXCombGroup() { }

void StdXCombGroup::clean() {
  StandardXComb::clean();
  theLastHeadCrossSection = ZERO;
  for ( vector<StdXCombPtr>::const_iterator dep = theDependent.begin();
	dep != theDependent.end(); ++dep )
    (**dep).clean();
}

int StdXCombGroup::nDim() const { 
  if ( meGroup()->willProject() )
    return StandardXComb::nDim() + 1;
  return StandardXComb::nDim();
}

CrossSection StdXCombGroup::dSigDR(const pair<double,double> ll, int nr, const double * r) {

  if ( matrixElement()->keepRandomNumbers() ) {
    lastRandomNumbers().resize(nDim());
    copy(r,r+nDim(),lastRandomNumbers().begin());
  }

  pExtractor()->select(this);
  setPartonBinInfo();
  lastP1P2(ll);
  lastS(sqr(maxEnergy())/exp(lastP1() + lastP2()));

  meMomenta().resize(mePartonData().size());
  matrixElement()->setXComb(this);

  double rProject = 0.0;
  size_t shift = 0;
  if ( meGroup()->willProject() ) {
    rProject = r[partonDims.first];
    shift = 1;
  }

  PPair partons;

  if ( !matrixElement()->haveX1X2() ) {

    if ( !pExtractor()->generateL(partonBinInstances(),
				  r, r + nr - partonDims.second) ) {
      lastCrossSection(ZERO);
      return ZERO;
    }
    partons = make_pair(partonBinInstances().first->parton(),
			partonBinInstances().second->parton());
    lastSHat(lastS()/exp(partonBinInstances().first->l() +
			 partonBinInstances().second->l()));
    meMomenta()[0] = partons.first->momentum();
    meMomenta()[1] = partons.second->momentum();

  } else {
    if ( !matrixElement()->generateKinematics(r + partonDims.first + shift) ) {
      lastCrossSection(ZERO);
      return ZERO;
    }
    lastSHat((meMomenta()[0]+meMomenta()[1]).m2());
    matrixElement()->setKinematics();

    lastScale(matrixElement()->scale());

    partons.first = mePartonData()[0]->produceParticle(meMomenta()[0]);
    partons.second = mePartonData()[1]->produceParticle(meMomenta()[1]);

    Direction<0> dir(true);
    partonBinInstances().first = 
      new_ptr(PartonBinInstance(lastParticles().first,partons.first,
				partonBins().first,lastScale()));
    dir.reverse();
    partonBinInstances().second = 
      new_ptr(PartonBinInstance(lastParticles().second,partons.second,
				partonBins().second,lastScale()));
  }

  lastPartons(partons);

  if ( lastSHat()  < cuts()->sHatMin() ) {
    lastCrossSection(ZERO);
    return ZERO;
  }

  lastY(0.5*(partonBinInstances().second->l() -
	     partonBinInstances().first->l()));
  if ( !cuts()->initSubProcess(lastSHat(), lastY(), mirror()) ) {
    lastCrossSection(ZERO);
    return ZERO;
  }

  if ( mirror() ) swap(meMomenta()[0], meMomenta()[1]);
  if ( matrixElement()->wantCMS() &&
       !matrixElement()->haveX1X2() ) 
    SimplePhaseSpace::CMS(meMomenta()[0], meMomenta()[1], lastSHat());

  Energy summ = ZERO;
  if ( meMomenta().size() == 3 ) {
    if ( !matrixElement()->haveX1X2() )
      meMomenta()[2] = Lorentz5Momentum(sqrt(lastSHat()));
  } else {
    for ( int i = 2, N = meMomenta().size(); i < N; ++i ) {
      if ( !matrixElement()->haveX1X2() )
	meMomenta()[i] = Lorentz5Momentum(mePartonData()[i]->mass());
      summ += mePartonData()[i]->massMin();
    }
    if ( sqr(summ) >= lastSHat() ) {
      lastCrossSection(ZERO);
      return ZERO;
    }
  }

  if ( !matrixElement()->haveX1X2() )
    lastScale(max(lastSHat()/4.0, cuts()->scaleMin()));

  lastSHat(pExtractor()->generateSHat(lastS(), partonBinInstances(),
				      r, r + nr - partonDims.second,
				      matrixElement()->haveX1X2()));

  if ( !cuts()->sHat(lastSHat()) ) {
    lastCrossSection(ZERO);
    return ZERO;
  }

  r += partonDims.first + shift;

  lastX1X2(make_pair(lastPartons().first->momentum().plus()/
		     lastParticles().first->momentum().plus(),
		     lastPartons().second->momentum().minus()/
		     lastParticles().second->momentum().minus()));

  if ( !cuts()->x1(lastX1()) || !cuts()->x2(lastX2()) ) {
    lastCrossSection(ZERO);
    return ZERO;
  }
  
  lastY((lastPartons().first->momentum() +
	 lastPartons().second->momentum()).rapidity());
  if ( !cuts()->yHat(lastY()) ) {
    lastCrossSection(ZERO);
    return ZERO;
  }
  if ( !cuts()->initSubProcess(lastSHat(), lastY(), mirror()) ) {
    lastCrossSection(ZERO);
    return ZERO;
  }

  meMomenta()[0] = lastPartons().first->momentum();
  meMomenta()[1] = lastPartons().second->momentum();
  if ( mirror() ) swap(meMomenta()[0], meMomenta()[1]);
  if ( matrixElement()->wantCMS() &&
       !matrixElement()->haveX1X2() ) 
    SimplePhaseSpace::CMS(meMomenta()[0], meMomenta()[1], lastSHat());

  if ( meMomenta().size() == 3 ) {
    if ( !matrixElement()->haveX1X2() )
      meMomenta()[2] = Lorentz5Momentum(sqrt(lastSHat()));
  } else {
    if ( sqr(summ) >= lastSHat() ) {
      lastCrossSection(ZERO);
      return ZERO;
    }
  }

  if ( !matrixElement()->haveX1X2() ) {
    if ( !matrixElement()->generateKinematics(r) ) {
      lastCrossSection(ZERO);
      return ZERO;
    }
  }
  lastScale(matrixElement()->scale());
  if ( !cuts()->scale(lastScale()) ) {
    lastCrossSection(ZERO);
    return ZERO;
  }

  pair<bool,bool> evalPDFS = 
    make_pair(matrixElement()->havePDFWeight1(),
	      matrixElement()->havePDFWeight2());
  if ( mirror() )
    swap(evalPDFS.first,evalPDFS.second);
  lastPDFWeight(pExtractor()->fullFn(partonBinInstances(), lastScale(),
				     evalPDFS));
  if ( lastPDFWeight() == 0.0 ) {
    lastCrossSection(ZERO);
    return ZERO;
  }
  matrixElement()->setKinematics();
  CrossSection xsec = matrixElement()->dSigHatDR() * lastPDFWeight();

  bool noHeadPass = !willPassCuts() || xsec == ZERO;
  if ( noHeadPass ) {
    xsec = ZERO;
    lastCrossSection(ZERO);
  }

  xsec *= cutWeight();

  lastAlphaS(matrixElement()->alphaS());
  lastAlphaEM(matrixElement()->alphaEM());

  lastHeadCrossSection(xsec);

  CrossSection depxsec = ZERO;

  vector<tStdXCombPtr> activeXCombs;

  for ( vector<StdXCombPtr>::const_iterator dep = theDependent.begin();
	dep != theDependent.end(); ++dep ) {
    if ( noHeadPass && (**dep).matrixElement()->headCuts() )
      continue;
    depxsec += (**dep).dSigDR(r + theMEGroup->dependentOffset((**dep).matrixElement()));
    if ( theMEGroup->groupReweighted() )
      activeXCombs.push_back(*dep);
  }

  if ( xsec != ZERO &&
       depxsec != ZERO ) {
    if ( abs(xsec) < abs(depxsec) ) {
      volatile double rw = (1.+depxsec/xsec);
      xsec = xsec*rw;
    } else {
      volatile double rw = (1.+xsec/depxsec);
      xsec = depxsec*rw;
    }
  } else if ( xsec == ZERO &&
	      depxsec != ZERO ) {
    xsec = depxsec;
  }

  lastCrossSection(xsec);

  if ( xsec != ZERO )
    theMEGroup->lastEventStatistics();

  matrixElement()->fillProjectors();
  if ( !projectors().empty() ) {
    lastProjector(projectors().select(rProject));
  }

  if ( theMEGroup->groupReweighted() && !activeXCombs.empty() ) {
    xsec = theMEGroup->reweightHead(activeXCombs)*lastHeadCrossSection();
    depxsec = ZERO;
    for ( vector<tStdXCombPtr>::const_iterator dep = activeXCombs.begin();
	  dep != activeXCombs.end(); ++dep ) {
      depxsec += theMEGroup->reweightDependent(*dep,activeXCombs)*(**dep).lastCrossSection();
    }
    if ( xsec != ZERO ) {
      double rw = 1.0 + depxsec/xsec;
      xsec *= rw;
    } else {
      xsec = depxsec;
    }
  }

  subProcess(SubProPtr());
  if ( CKKWHandler() && matrixElement()->maxMultCKKW() > 0 &&
       matrixElement()->maxMultCKKW() > matrixElement()->minMultCKKW() ) {
    newSubProcess(theMEGroup->subProcessGroups());
    CKKWHandler()->setXComb(this);
    xsec *= CKKWHandler()->reweightCKKW(matrixElement()->minMultCKKW(),
					matrixElement()->maxMultCKKW());
  }

  if ( matrixElement()->reweighted() ) {
    newSubProcess(theMEGroup->subProcessGroups());
    xsec *= matrixElement()->reWeight() * matrixElement()->preWeight();
  }

  return xsec;

}

void StdXCombGroup::newSubProcess(bool) {

  StandardXComb::newSubProcess(theMEGroup->subProcessGroups());

  if ( !theMEGroup->subProcessGroups() )
    return;

  subProcess()->groupWeight(lastHeadCrossSection()/lastCrossSection());

  Ptr<SubProcessGroup>::tptr group = 
    dynamic_ptr_cast<Ptr<SubProcessGroup>::tptr>(subProcess());
  assert(group);

  for ( vector<StdXCombPtr>::iterator dep = theDependent.begin();
	dep != theDependent.end(); ++dep ) {
    if ( (**dep).lastCrossSection() == ZERO )
      continue;
    tSubProPtr ds;
    try {
      ds = (**dep).construct();
    } catch(Veto&) {
      throw Exception() << "A veto was encountered while constructing a dependent sub process.\n"
			<< "This situation should not have happened. Will veto the event."
			<< Exception::eventerror;
    }
    if ( ds )
      group->add(ds);
  }

}

tSubProPtr StdXCombGroup::construct() {

  matrixElement()->setXComb(this);

  setPartonBinInfo();
  matrixElement()->setKinematics();

  newSubProcess(theMEGroup->subProcessGroups());

  if ( !theMEGroup->subProcessGroups() ) {
    if ( !cuts()->initSubProcess(lastSHat(), lastY()) ) throw Veto();
    TmpTransform<tSubProPtr>
      tmp(subProcess(), Utilities::getBoostToCM(subProcess()->incoming()));
    if ( !cuts()->passCuts(*subProcess()) ) throw Veto();
  }

  return subProcess();

}


void StdXCombGroup::Init() {}

void StdXCombGroup::persistentOutput(PersistentOStream & os) const {
  os << theDependent << theMEGroup << ounit(theLastHeadCrossSection,nanobarn);
}

void StdXCombGroup::persistentInput(PersistentIStream & is, int) {
  is >> theDependent >> theMEGroup >> iunit(theLastHeadCrossSection,nanobarn);
}

ClassDescription<StdXCombGroup> StdXCombGroup::initStdXCombGroup;

