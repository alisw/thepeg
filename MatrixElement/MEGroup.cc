// -*- C++ -*-
//
// MEGroup.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
// Copyright (C) 2009-2017 Simon Platzer
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEGroup class.
//

#include "MEGroup.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/PDF/PartonBin.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/StdXCombGroup.h"

using namespace ThePEG;

MEGroup::MEGroup()
  : theDependent(), theNDimMap(), theNDim(0) {}

MEGroup::~MEGroup() {}

void MEGroup::doinit() {
  MEBase::doinit();
  head()->init();
  for ( MEVector::iterator me = theDependent.begin();
	me != theDependent.end(); ++me )
    (**me).init();
  use(head());
  theNDim = head()->nDim();
  if (!uniformAdditional()) {
    int off = theNDim;
    for ( MEVector::iterator me = theDependent.begin();
	  me != theDependent.end(); ++me )
      if ( (**me).nDim() > head()->nDim() ) {
	theNDimMap[*me] = off;
	off += ((**me).nDim() - head()->nDim());
	theNDim += ((**me).nDim() - head()->nDim());
      }
  } else {
    int maxadd = 0;
    for ( MEVector::iterator me = theDependent.begin();
	  me != theDependent.end(); ++me )
      if ( (**me).nDim() - head()->nDim() > maxadd ) {
	maxadd = ((**me).nDim() - head()->nDim());
      }
    theNDim += maxadd;
  }
}

void MEGroup::doinitrun() {
  MEBase::doinitrun();
  head()->initrun();
  for ( MEVector::iterator me = theDependent.begin();
	me != theDependent.end(); ++me )
    (**me).initrun();
}

void MEGroup::rebind(const TranslationMap & trans) {
  map<tMEPtr,int> rebound;
  for (map<tMEPtr,int>::iterator it = theNDimMap.begin();
       it != theNDimMap.end(); ++it) {
    rebound.insert(make_pair(trans.translate(it->first),it->second));
  }
  theNDimMap = rebound;
  MEBase::rebind(trans);
}

IVector MEGroup::getReferences() {
  IVector ret = MEBase::getReferences();
  for (map<tMEPtr,int>::iterator it = theNDimMap.begin();
       it != theNDimMap.end(); ++it)
    ret.push_back(it->first);
  return ret;
}

int MEGroup::dependentOffset(tMEPtr dep) const {
  if ( uniformAdditional() )
    return head()->nDim();
  map<tMEPtr,int>::const_iterator it =
    theNDimMap.find(dep);
  if (it == theNDimMap.end())
    return 0;
  return it->second;
}

StdXCombPtr MEGroup::makeXComb(Energy newMaxEnergy, const cPDPair & inc,
			       tEHPtr newEventHandler,tSubHdlPtr newSubProcessHandler,
			       tPExtrPtr newExtractor,	tCascHdlPtr newCKKW,
			       const PBPair & newPartonBins, tCutsPtr newCuts,
			       const DiagramVector & newDiagrams, bool mir,
			       const PartonPairVec& allPBins,
			       tStdXCombPtr newHead,
			       tMEPtr newME) {
  tMEGroupPtr newMEGroup = dynamic_ptr_cast<tMEGroupPtr>(newME);
  if ( !newMEGroup )
    newMEGroup = this;
  StdXCombGroupPtr res =  new_ptr(StdXCombGroup(newMaxEnergy, inc,
						newEventHandler, newSubProcessHandler,
						newExtractor, newCKKW,
						newPartonBins, newCuts, newMEGroup,
						newDiagrams, mir,
						newHead));
  res->build(allPBins);
  return res;
}

vector<StdXCombPtr> MEGroup::makeDependentXCombs(tStdXCombPtr xcHead,
						 const cPDVector& proc,
						 tMEPtr depME,
						 const PartonPairVec& pbs) const {
  
  MEBase::DiagramVector depDiags = dependentDiagrams(proc,depME);

  if ( depDiags.empty() )
    return vector<StdXCombPtr>();

  map<cPDVector,MEBase::DiagramVector> depProcs;

  for ( MEBase::DiagramVector::const_iterator d = depDiags.begin();
	d != depDiags.end(); ++d ) {
    depProcs[(**d).partons()].push_back(*d);
  }

  vector<StdXCombPtr> ret;

  for ( map<cPDVector,MEBase::DiagramVector>::const_iterator pr =
	  depProcs.begin(); pr != depProcs.end(); ++pr ) {


    PartonPairVec::const_iterator ppit = pbs.begin();
    for ( ; ppit != pbs.end(); ++ppit ) {
      if ( ppit->first->parton() == pr->second.front()->partons()[0] &&
	   ppit->second->parton() == pr->second.front()->partons()[1] )
	break;
    }

    if ( ppit == pbs.end() ) {
      generator()->logWarning(
			      Exception() 
			      << "Could not create a dependent XComb object"
			      << " for the MEGroup '"
			      << name() 
			      << "' since the dependent matrix element '"
			      << depME->name() 
			      << "' did not match any of the incoming partons."
			      << Exception::warning);
      continue;
    }

    StdXCombPtr dxc = depME->makeXComb(xcHead,*ppit,pr->second);
    ret.push_back(dxc);

  }

  return ret;

}

bool MEGroup::generateKinematics(const double * r) {
  if (!head()->generateKinematics(r))
    return false;
  return true;
}

void MEGroup::clearKinematics() {
  MEBase::clearKinematics();
  head()->clearKinematics();
  for ( MEVector::iterator me = theDependent.begin();
	me != theDependent.end(); ++me )
    (**me).clearKinematics();
}

void MEGroup::persistentOutput(PersistentOStream & os) const {
  os << theHead << theDependent << theNDimMap << theNDim;
}

void MEGroup::persistentInput(PersistentIStream & is, int) {
  is >> theHead >> theDependent >> theNDimMap >> theNDim;
}

AbstractClassDescription<MEGroup> MEGroup::initMEGroup;
// Definition of the static class description member.

void MEGroup::Init() {

  static ClassDocumentation<MEGroup> documentation
    ("The ThePEG::MEGroup class is the base class for all matrix elements "
     "to be used for generating sub process groups in ThePEG");

  static Reference<MEGroup,MEBase> interfaceHead
    ("Head",
     "The head matrix element for this matrix element group.",
     &MEGroup::theHead, false, false, true, false, false);


  static RefVector<MEGroup,MEBase> interfaceDependent
    ("Dependent",
     "The vector of dependent matrix elements in this matrix element group.",
     &MEGroup::theDependent, -1, false, false, true, false, false);

}

