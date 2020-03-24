// -*- C++ -*-
//
// TemplateTools.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2006-2019 David Grellscheid, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Template_Tools_H
#define Template_Tools_H

#include <type_traits>

/**
 * @file TemplateTools.h 
 * Useful template machinery. Based on Alexandrescu, "Modern C++ Design".
 */

namespace ThePEG {

/// Conversion between integers and types.
template <int v>
struct Int2Type
{
  enum { value = v };
};

/// Dummy type for ambiguous function signatures.
struct DummyType {};

/** @endcond */

/// Selection mechanism for type-dependent implementations.
enum ImplSelector { Dimensioned, Standard, Enumerated };

/// Typedef for dimensioned types.
typedef Int2Type<Dimensioned> DimensionT;

/// Typedef for non-dimensioned types.
typedef Int2Type<Standard> StandardT;

/// Typedef for non-dimensioned types.
typedef Int2Type<Enumerated> EnumT;


/// Type traits for built-in types
template <typename T>
struct TypeTraits
{
  /// Boolean flag. Is true for physical quantities.
  enum { hasDimension = false };
  /// Implementation selector
  typedef StandardT DimType;
  /// Base unit for arithmetic types
  // construction with extra U type is necessary to make
  // enable_if work before concepts are supported properly
  template <typename U = T>
  static constexpr typename
  std::enable_if< (std::is_arithmetic<U>::value && 
                   std::is_same<U, T>::value), U>::type 
  baseunit() 
  { return static_cast<U>(1); }
};

}

#endif
