// -*- C++ -*-
#ifndef THEPEG_TestDecays_H
#define THEPEG_TestDecays_H
//
// This is the declaration of the <!id>TestDecays<!!id> class.
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
// #include "TestDecays.fh"
// #include "TestDecays.xh"

namespace ThePEG {

class TestDecays: public Main {

public:

  inline TestDecays();
  inline TestDecays(const TestDecays &);
  virtual ~TestDecays();
  // Standard ctors and dtor.

public:

  static void testDecay(tEGPtr eg, long id, int N);
  static void testGG(tEGPtr eg, Energy W, int N);

  static void Init();
  // Standard Init function used to initialize the interfaces.

private:

  static NoPIOClassDescription<TestDecays> initTestDecays;
  // Describe a concrete class without persistent data.

  TestDecays & operator=(const TestDecays &);
  // Private and non-existent assignment operator.

};

}

// CLASSDOC OFF

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

// The following template specialization informs ThePEG about the
// base class of TestDecays.
template <>
struct BaseClassTrait<TestDecays,1>: public ClassTraitsType {
  typedef Main NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<TestDecays>
  : public ClassTraitsBase<TestDecays> {
  static string className() { return "ThePEG::TestDecays"; }
  // Return the class name.
  static string library() { return "TestDecays.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

/** @endcond */

}

#include "TestDecays.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TestDecays.tcc"
#endif

#endif /* THEPEG_TestDecays_H */
