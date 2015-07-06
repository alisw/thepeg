// -*- C++ -*-
#ifndef THEPEG_TestJunctions_H
#define THEPEG_TestJunctions_H
//
// This is the declaration of the <!id>TestJunctions<!!id> class.
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
// #include "TestJunctions.fh"
// #include "TestJunctions.xh"

namespace ThePEG {

class TestJunctions: public Main {

public:

  inline TestJunctions();
  inline TestJunctions(const TestJunctions &);
  virtual ~TestJunctions();
  // Standard ctors and dtor.

public:

  static void testSimple(tEGPtr eg, int N);
  static void testSimpleG(tEGPtr eg, int N);
  static void testConnected(tEGPtr eg, int N);
  static void testConnectedG(tEGPtr eg, int N);

  static void Init();
  // Standard Init function used to initialize the interfaces.

private:

  static NoPIOClassDescription<TestJunctions> initTestJunctions;
  // Describe a concrete class without persistent data.

  TestJunctions & operator=(const TestJunctions &);
  // Private and non-existent assignment operator.

};

}

// CLASSDOC OFF

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

// The following template specialization informs ThePEG about the
// base class of TestJunctions.
template <>
struct BaseClassTrait<TestJunctions,1>: public ClassTraitsType {
  typedef Main NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<TestJunctions>
  : public ClassTraitsBase<TestJunctions> {
  static string className() { return "ThePEG::TestJunctions"; }
  // Return the class name.
  static string library() { return "TestJunctions.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

/** @endcond */

}

#include "TestJunctions.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TestJunctions.tcc"
#endif

#endif /* THEPEG_TestJunctions_H */
