// -*- C++ -*-
//
// ThePEG.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_H
#define ThePEG_H

/** \file ThePEG.h
 * This is the main config header file for ThePEG. Do not make
 * changes in this file. If you need to modify anything, edit a copy
 * of the file which can be included instead of this file using the
 * macro <code>ThePEG_ALTERNATE_CONFIG</code>.
 */

#ifndef ThePEG_ALTERNATE_CONFIG

#include "ThePEG/Pointer/Ptr.h"
#include "ThePEG/Pointer/PtrTraits.h"
#include "ThePEG/Pointer/RCPtr.h"
#include "ThePEG/Utilities/Rebinder.fh"
#include "ThePEG/Utilities/Interval.fh"
#include "ThePEG/Utilities/ClassDescription.fh"
#include "ThePEG/Interface/InterfaceBase.fh"
#include "ThePEG/Persistency/PersistentOStream.fh"
#include "ThePEG/Persistency/PersistentIStream.fh"
#include "TemplateTools.h"
#include "Complex.h"
#include "Unitsystem.h"
#include "Constants.h"
#include "std.h"


/**
 * This is the main namespace within which all identifiers in ThePEG
 * are declared. External packages based on ThePEG should not
 * introduce identifiers in the ThePEG namespace, but in a separate
 * namespace which need not be nested within the ThePEG namespace.
 */
namespace ThePEG {

// Introduce some identifiers in the ThePEG namespace/
using namespace ThePEG::Pointer;
using ThePEG::Pointer::Ptr;
using namespace ThePEG::Units;

/**
 * Define the base class from which all (polymorphic) classes in
 * ThePEG are derived.
 */
struct Base: public ReferenceCounted {

  /** The virtual destructor */
  virtual ~Base() {}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init() {}

  /**
   * Print out debugging information for this object on std::cerr. To
   * be called from within a debugger. Simply calls the virtual
   * debugme() function.
   */
  void debug() const;

  /**
   * Print out debugging information for this object on std::cerr. To
   * be called from within a debugger via the debug() function.
   */
  virtual void debugme() const;

};

/**
 * Define the base class from which all persistent classes in
 * ThePEG are derived.
 */
typedef Base PersistentBase;

/**
 * TraitsType is an empty, non-polymorphic, base class. It is used as
 * a base class of all traits classes in ThePEG in order to group them
 * together in the documentation. It currently serves no other
 * purpose.
 */
struct TraitsType {};

/**
 * A standard exception class to be used for vetoing a whole event.
 */
struct Veto {

  /** the default constructor. */
  Veto();
};

/**
 * A standard exception class to be used to temporarily stop the
 * generation of an event.
 */
struct Stop {};

/**
 * The square function should really have been included in the
 * standard C++ library.
 */
template <typename T>
inline typename BinaryOpTraits<T,T>::MulT sqr(const T& x) {
  return x*x;
}

/**
 * The templated abs function should really have been included in the
 * standard C++ library.
 */
template <typename T>
inline T abs(const T & t) {
  return t < T()? -t: t;
}

// Debugging in ThePEG may be swithced off completely by this
// compilation swithc, eliminating possible overhead in error
// checking.
#ifndef ThePEG_NO_DEBUG
/** Macro for accessing debug functions to enable compile-time disabling. */
#define ThePEG_DEBUG_LEVEL Debug::level
/** Macro for accessing debug functions to enable compile-time disabling. */
#define ThePEG_DEBUG_ITEM(item) Debug::debugItem(item)
#else
/** Macro for accessing debug functions to enable compile-time disabling. */
#define ThePEG_DEBUG_LEVEL 0
/** Macro for accessing debug functions to enable compile-time disabling. */
#define ThePEG_DEBUG_ITEM(item) false
#endif

}

#include "Pointers.h"
#include "Containers.h"

#else

#include ThePEG_ALTERNATIVE_CONFIG

#endif

#endif /* ThePEG_H */

