// -*- C++ -*-
#ifndef THEPEG_TestMisc_H
#define THEPEG_TestMisc_H
//
// This is the declaration of the <!id>TestMisc<!!id> class.
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
// #include "TestMisc.fh"
// #include "TestMisc.xh"

namespace ThePEG {

class TestMisc: public Main {

public:

  inline TestMisc();
  inline TestMisc(const TestMisc &);
  virtual ~TestMisc();
  // Standard ctors and dtor.

public:

  static void Init();
  // Standard Init function used to initialize the interfaces.

private:

  static NoPIOClassDescription<TestMisc> initTestMisc;
  // Describe a concrete class without persistent data.

  TestMisc & operator=(const TestMisc &);
  // Private and non-existent assignment operator.

};

}

// CLASSDOC OFF

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

// The following template specialization informs ThePEG about the
// base class of TestMisc.
template <>
struct BaseClassTrait<TestMisc,1>: public ClassTraitsType {
  typedef Main NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<TestMisc>
  : public ClassTraitsBase<TestMisc> {
  static string className() { return "ThePEG::TestMisc"; }
  // Return the class name.
  static string library() { return "TestMisc.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

/** @endcond */

}

#include "TestMisc.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TestMisc.tcc"
#endif

#endif /* THEPEG_TestMisc_H */
