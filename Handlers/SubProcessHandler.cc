// -*- C++ -*-
//
// SubProcessHandler.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SubProcessHandler class.
//

#include "SubProcessHandler.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/Handlers/Hint.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/CascadeHandler.h"
#include "ThePEG/Handlers/MultipleInteractionHandler.h"
#include "ThePEG/Handlers/HadronizationHandler.h"
#include "ThePEG/Handlers/DecayHandler.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/MatrixElement/ReweightBase.h"

using namespace ThePEG;

SubProcessHandler::~SubProcessHandler() {}

SubProcessHandler::SubProcessHandler() {
  setupGroups();
}

SubProcessHandler::SubProcessHandler(const SubProcessHandler & sph)
  : HandlerBase(sph),
    thePartonExtractor(sph.thePartonExtractor), theMEs(sph.theMEs),
    theCuts(sph.theCuts), theSubprocessGroup(sph.theSubprocessGroup),
    theCascadeGroup(sph.theCascadeGroup), theMultiGroup(sph.theMultiGroup),
    theHadronizationGroup(sph.theHadronizationGroup),
    theDecayGroup(sph.theDecayGroup),
    reweights(sph.reweights), preweights(sph.preweights) {
  setupGroups();
}

void SubProcessHandler::setupGroups() {
  theGroups.clear();
  theGroups.push_back(&theSubprocessGroup);
  theGroups.push_back(&theCascadeGroup);
  theGroups.push_back(&theMultiGroup);
  theGroups.push_back(&theHadronizationGroup);
  theGroups.push_back(&theDecayGroup);
}

IBPtr SubProcessHandler::clone() const {
  return new_ptr(*this);
}

IBPtr SubProcessHandler::fullclone() const {
  return new_ptr(*this);
}

const HandlerGroupBase &
SubProcessHandler::handlerGroup(Group::Handler group) const {
  return *(theGroups[group]);
}

tCascHdlPtr SubProcessHandler::CKKWHandler() const {
  return dynamic_ptr_cast<tCascHdlPtr>(theCascadeGroup.defaultHandler());
}

void SubProcessHandler::persistentOutput(PersistentOStream & os) const {
  os << thePartonExtractor << theCuts << theSubprocessGroup << theCascadeGroup
     << theMultiGroup << theHadronizationGroup << theDecayGroup
     << theMEs << reweights << preweights;
}

void SubProcessHandler::persistentInput(PersistentIStream & is, int) {
  is >> thePartonExtractor >> theCuts >> theSubprocessGroup >> theCascadeGroup
     >> theMultiGroup >> theHadronizationGroup >> theDecayGroup
     >> theMEs >> reweights >> preweights;
}

void SubProcessHandler::doinit() {
  for ( MEVector::iterator me = theMEs.begin(); me != theMEs.end(); ++me ) {
    (**me).init();
    for ( ReweightVector::iterator i = reweights.begin();
	  i != reweights.end(); ++i ) (**me).addReweighter(*i);
    for ( ReweightVector::iterator i = preweights.begin();
	  i != preweights.end(); ++i ) (**me).addPreweighter(*i);
  }
  HandlerBase::doinit();
}

void SubProcessHandler::doinitrun() {
  HandlerBase::doinitrun();
  pExtractor()->initrun();
  for ( MEVector::iterator me = theMEs.begin(); me != theMEs.end(); ++me )
    (**me).initrun();
  for ( ReweightVector::iterator i = reweights.begin();
	i != reweights.end(); ++i ) (**i).initrun();
  for ( ReweightVector::iterator i = preweights.begin();
	i != preweights.end(); ++i ) (**i).initrun();
}

ClassDescription<SubProcessHandler> SubProcessHandler::initSubProcessHandler;

void SubProcessHandler::Init() {

  static ClassDocumentation<SubProcessHandler> documentation
    ("This object contains information about a set of possible sub-processes "
     "to be generated from inside ThePEG. It must contain a "
     "<interface>PartonExtractor</interface> do describe how the partons "
     "entering into the hard sub-process are extracted from the beam "
     "particles. It must also include at least one matrix element object "
     "in <interface>MatrixElements</interface> and a "
     "<interface>Cuts</interface> object describing the kinematical cuts "
     "imposed on the sub-process generation.");

  static Reference<SubProcessHandler,PartonExtractor> interfacePartonExtractor
    ("PartonExtractor",
     "The PartonExtractor object to describe the way partons are extracted "
     "from the incoming particles.",
     &SubProcessHandler::thePartonExtractor, false, false, true, false);

  static RefVector<SubProcessHandler,MEBase> interfaceMEs
    ("MatrixElements",
     "A list of MEBase objects describing the \\f$2\\rightarrow n\\f$ hard "
     "matrix elements.",
     &SubProcessHandler::theMEs, 0, false, false, true, false);

  static Reference<SubProcessHandler,Cuts> interfaceCuts
    ("Cuts",
     "Common kinematical cuts for this SubProcessHandler. These cuts "
     "overides those in a EventHandler.",
     &SubProcessHandler::theCuts, false, false, true, true);

  static RefVector<SubProcessHandler,ReweightBase> interfaceReweights
    ("Reweights",
     "A list of ThePEG::ReweightBase objects to modify all matrix elements "
     "in this SubProcessHandler.",
     &SubProcessHandler::reweights, 0, false, false, true, false);

  static RefVector<SubProcessHandler,ReweightBase> interfacePreweights
    ("Preweights",
     "A list of ThePEG::ReweightBase objects to bias the phase space for all "
     "matrix elements without in this SubProcessHandler influencing the "
     "actual cross section.",
     &SubProcessHandler::preweights, 0, false, false, true, false);

  ThePEG_DECLARE_PREPOST_OBJECTS(SubProcessHandler, SubProcessHandler,
				  Post, after);

  ThePEG_DECLARE_GROUPINTERFACE_OBJECTS(SubProcessHandler, CascadeHandler);
  ThePEG_DECLARE_GROUPINTERFACE_OBJECTS(SubProcessHandler,
					 MultipleInteractionHandler);
  ThePEG_DECLARE_GROUPINTERFACE_OBJECTS(SubProcessHandler,
					 HadronizationHandler);
  ThePEG_DECLARE_GROUPINTERFACE_OBJECTS(SubProcessHandler, DecayHandler);

  interfacePartonExtractor.rank(10);
  interfaceMEs.rank(9);
  interfaceCuts.rank(8);

}


ThePEG_IMPLEMENT_PREPOST_GROUP(SubProcessHandler,SubProcessHandler,
				theSubprocessGroup,Post)\

ThePEG_IMPLEMENT_GROUPINTERFACE(SubProcessHandler,CascadeHandler,
				 theCascadeGroup,CascHdlPtr) \

ThePEG_IMPLEMENT_GROUPINTERFACE(SubProcessHandler,MultipleInteractionHandler,
				 theMultiGroup,MIHdlPtr)  \

ThePEG_IMPLEMENT_GROUPINTERFACE(SubProcessHandler,HadronizationHandler,
				 theHadronizationGroup,HadrHdlPtr) \

ThePEG_IMPLEMENT_GROUPINTERFACE(SubProcessHandler,DecayHandler,
				 theDecayGroup,DecayHdlPtr) \

