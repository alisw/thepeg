// -*- C++ -*-
//
// Containers.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

// This file contains the implementations of the container
// declarations in Containers.h.

#include "ThePEG.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/MatcherBase.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/InterfacedBase.h"
#include "ThePEG/Repository/EventGenerator.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Containers.tcc"
#endif

using namespace ThePEG;

ThePEG_IMPLEMENT_SET(PDPTr,ParticleDataSet)
ThePEG_IMPLEMENT_SET(PMPtr,MatcherSet)
ThePEG_IMPLEMENT_SET(DMPtr,DecayModeSet)
ThePEG_IMPLEMENT_SET(tDMPtr,DecaySet)
ThePEG_IMPLEMENT_SET(IBPtr,ObjectSet)
ThePEG_IMPLEMENT_SET(IBPtr,DependencySet)
ThePEG_IMPLEMENT_MAP(long,PDPtr,ParticleMap)
ThePEG_IMPLEMENT_MAP(string,IBPtr,ObjectMap)
ThePEG_IMPLEMENT_MAP(IBPtr,DependencySet,DependencyMap)
ThePEG_IMPLEMENT_MAP(string,const InterfaceBase *,InterfaceMap)
ThePEG_IMPLEMENT_MAP(string,EGPtr,GeneratorMap)
ThePEG_IMPLEMENT_SET(string,StringSet)
// typedef set<const InterfaceBase *> InterfaceSet;
