// -*- C++ -*-
//
// TestLHAPDF.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_TestLHAPDF_H
#define ThePEG_TestLHAPDF_H
//
// This is the declaration of the <!id>TestLHAPDF<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:.html">.h</a>,
// <a href="http:.html">.h</a>.
// 

#include "ThePEG/Repository/Main.h"

namespace ThePEG {

/**
 * Class derived from Main to do simple tests of the LHAPDF class.
 */
class TestLHAPDF: public Main {

public:

  /**
   * Standard Init function used to initialize the interfaces. This is
   * where all the tests are done.
   */
  static void Init();

private:

  /**
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<TestLHAPDF> initTestLHAPDF;

  /**
   * Private and non-existent assignment operator.
   */
  TestLHAPDF & operator=(const TestLHAPDF &) = delete;

};

}

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

// The following template specialization informs ThePEG about the
// base class of TestLHAPDF.
template <>
struct BaseClassTrait<ThePEG::TestLHAPDF,1>: public ClassTraitsType {
  typedef Main NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<ThePEG::TestLHAPDF>
  : public ClassTraitsBase<ThePEG::TestLHAPDF> {
  static string className() { return "ThePEG::TestLHAPDF"; }
  // Return the class name.
  static string library() { return "TestLHAPDF.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

/** @endcond */

}

#endif /* ThePEG_TestLHAPDF_H */
