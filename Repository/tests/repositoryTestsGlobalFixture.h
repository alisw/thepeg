// -*- C++ -*-
//
// repositoryTestGlobalFixture.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad, 2015 Marco A. Harrendorf
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Repository_Tests_GlobalFixture_H
#define ThePEG_Repository_Tests_GlobalFixture_H

#include <boost/test/unit_test.hpp>

#include "ThePEG/Repository/StandardRandom.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Config/Unitsystem.h"

struct FixGlobal1 {
  ThePEG::StandardRandom srng;
  ThePEG::UseRandom urng;

  FixGlobal1() : srng(), urng(&srng) {
    BOOST_TEST_MESSAGE( "setup global fixture for repositoryTest" ); 
  }
  
  ~FixGlobal1() { 
  	BOOST_TEST_MESSAGE( "teardown global fixture for repositoryTest" ); 
  }
};

BOOST_GLOBAL_FIXTURE(FixGlobal1);

#endif
