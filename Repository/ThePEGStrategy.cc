// -*- C++ -*-
//
// ThePEGStrategy.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThePEGStrategy class.
//

#include "ThePEGStrategy.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/ParticleData.h"
#include "Repository.h"

using namespace ThePEG;

const std::string ThePEGStrategy::versionstring() const {
	return Repository::version();
}

IBPtr ThePEGStrategy::clone() const {
  return new_ptr(*this);
}

IBPtr ThePEGStrategy::fullclone() const {
  return new_ptr(*this);
}

NoPIOClassDescription<ThePEGStrategy> ThePEGStrategy::initThePEGStrategy;

void ThePEGStrategy::Init() {
  static ClassDocumentation<ThePEGStrategy> interfaceDescription
    ("This class represents the default ThePEG strategy", "", "");


}

