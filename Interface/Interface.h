// -*- C++ -*-
//
// Interface.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_Interface_H
#define THEPEG_Interface_H
//
// This is the declaration of the Interface namespace.
//

namespace ThePEG {

/**
 * The Interface namespace declares a number of enums to set options
 * for subclasses of InteracedBase.
 */
namespace Interface {

/**
 * Determine whether an interface is dependency safe or
 * not. Dependency safe means that other objects do not depend on the
 * variable being interfaced.
 */
enum DepSafe {
  unsafe = false, /**< The interface is not dependency safe. */
  safe = true     /**< The interface is dependency safe. */
};

/**
 * Determine whether an interface is read-only or not.
 */
enum ReadOnly {
  readwrite = false, /**< The interface is mutable. */
  readonly = true    /**< The interface is read-only. */
};

/**
 * Determine whether a Parameter or ParVector is limited, either
 * upper, lower or both.
 */
enum Limits {
  nolimits = 0, /**< The parameter is not limited. */
  limited = 1,  /**< The parameter is limited (both up- and downwards. */
  upperlim = 2, /**< The parameter has only an upper limit. */
  lowerlim = 3  /**< The parameter has only an lower limit. */
};

/**
 * Determine whether a the objects referred to by a Reference or a
 * RefVector should be automaticlly rebound (i.e. do not need to be
 * explicitly rebound in the rebind() function).
 */
enum Rebind {
  norebind = true, /**< The reference is not automatically rebound. */
  rebind = false     /**< The reference is automatically rebound. */
};

/**
 * Determine whether a Reference or RefVector object may be null.
 */
enum Nullable {
  nonull = false, /**< The reference may not be null. */
  nullok = true   /**< The reference may be null. */
};

/**
 * Determine whether a null reference should be given a default value
 * if suitable object is registered as default in the Strategy object
 * of a run.
 */
enum NullDefault {
  nodefnull = false, /**< The reference will not be set to default if null. */
  defnull = true     /**< The reference will be set to default if null. */
};


}

}

#endif /* THEPEG_Interface_H */
