// -*- C++ -*-
//
// repositoryTestMain.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad, 2015 Marco A. Harrendorf
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

/**
 * The following part should be included only once. 
*/
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#define BOOST_TEST_MODULE repositoryTest

/**
 * Include global fixture
 * 
 * Global fixture initializes the randomNumber generator
 */
#include "ThePEG/Repository/tests/repositoryTestsGlobalFixture.h"

/**
 * Include here the sub tests
 */
#include "ThePEG/Repository/tests/repositoryTestRandomGenerator.h"


/**
 * Debug and development part
 */
BOOST_AUTO_TEST_CASE(fail)
{
  //BOOST_FAIL("Ende");
}
