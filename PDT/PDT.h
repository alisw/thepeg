// -*- C++ -*-
//
// PDT.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_PDT_H
#define ThePEG_PDT_H
// This is the declaration of the PDT class.

#include "ThePEG/Config/ThePEG.h"

namespace ThePEG {

/**
 * PDT is a helper class implementing enumerations for charge, colour
 * and spin to be used by the ParticleData class. In addition, some
 * static utility functions are provided.
 *
 * @see ParticleData
 */
class PDT {

public:

  /**
   * Definition of enumerated values used for spin information. THe
   * integer values are given according to 2s+1.
   */
  enum Spin {
    SpinNA = -1,       /**< Spin is not applicable. */
    SpinUnknown = 0,   /**< Unknown spin */
    SpinUndefined = 0, /**< Undefined spin */
    Spin0 = 1,         /**< Spin zero. */
    Spin1Half = 2,     /**< Spin 1/2. */
    Spin1 = 3,         /**< Spin 1. */
    Spin3Half = 4,     /**< Spin 3/2. */
    Spin2 = 5,         /**< Spin 2. */
    Spin5Half = 6,     /**< Spin 5/2. */
    Spin3 = 7,         /**< Spin 4. */
    Spin7Half = 8,     /**< Spin 7/2. */
    Spin4 = 9          /**< Spin 5. */
  };

  /**
   * Definition of enumerated values used for charge information. The
   * integer values are given in units of e/3.
   */
  enum Charge {
    ChargeUnknown = -999999,   /**< Unknown charge. */
    ChargeUndefined = -999999, /**< Undefined charge. */
    Charged = 999990,          /**< Is charged. */ 
    Positive = 900000,         /**< Is positively charged. */ 
    Negative = -900000,        /**< Is negatively charged. */ 
    ChargeNeutral = 0,         /**< Uncharged. */ 
    Charge0 = 0,               /**< Uncharged. */ 
    Plus1Third = 1,            /**< e/3. */ 
    Plus2Third = 2,            /**< 2e/3. */ 
    Plus1 = 3,                 /**< e. */ 
    Minus1Third = -1,          /**< -e/3. */ 
    Minus2Third = -2,          /**< -2e/3. */ 
    Minus1 = -3,               /**< -e. */ 
    Plus4Third = 4,            /**< 4e/3. */ 
    Plus5Third = 5,            /**< 5e/3. */ 
    Plus2 = 6,                 /**< 2e. */ 
    Minus4Third = -4,          /**< -4e/3. */ 
    Minus5Third = -5,          /**< -5e/3. */ 
    Minus2 = -6,               /**< -3e. */ 
    Plus7Third = 7,            /**< 7e/3. */ 
    Plus8Third = 8,            /**< 8e/3. */ 
    Plus3 = 9,                 /**< 3e. */ 
    Minus7Third = -7,          /**< -7e/3. */ 
    Minus8Third = -8,          /**< -8e/3. */ 
    Minus3 = -9,               /**< -3e. */ 
    Plus4 = 12,                /**< 4e. */ 
    Plus5 = 15,                /**< 5e. */ 
    Plus6 = 18,                /**< 6e. */ 
    Plus7 = 21,                /**< 7e. */ 
    Plus8 = 24,                /**< 8e. */ 
    Minus4 = -12,              /**< -4e. */ 
    Minus5 = -15,              /**< -5e. */ 
    Minus6 = -18,              /**< -6e. */ 
    Minus7 = -21,              /**< -7e. */ 
    Minus8 = -24               /**< -8e. */ 
  };

  /**
   *Definition of enumerated values used for colour information.
   */
  enum Colour {
    ColourUnknown = -1,   /**< Unknown colour */
    ColourUndefined = -1, /**< Undefined colour */
    ColourNeutral = 0,    /**< Colour-singlet */
    Colour0 = 0,          /**< Colour-singlet */
    Coloured = 1,         /**< Coloured */
    Colour3 = 3,          /**< Colour-triplet */
    Colour3bar = -3,      /**< Colour-anti-triplet */
    Colour6 = 6,          /**< Colour-sextet */
    Colour6bar = -6,      /**< Colour-anti-sextet */
    Colour8 = 8           /**< Colour-octet */
  };

  /**
   * True if the argument corresponds to a non-zero charge.
   */
  static bool charged(Charge c) {
    return c != ChargeNeutral && c != ChargeUndefined;
  }

  /**
   * True if the argument corresponds to a positive charge.
   */
  static bool positive(Charge c) {
    return c > ChargeNeutral && c != Charged;
  }

  /**
   * True if the argument corresponds to a negative charge.
   */
  static bool negative(Charge c) {
    return c < ChargeNeutral && c != ChargeUndefined;
  }

  /**
   * True if the argument corresponds to a non-zero colour charge.
   */
  static bool coloured(Colour c) {
    return c != ColourNeutral && c != ColourUnknown;
  }

  /**
   * Return the anti-colour of the specified colour.
   */
  static Colour antiColour(Colour c) {
    if ( c == Colour3 || c == Colour3bar ) return Colour(-c);
    if ( c == Colour6 || c == Colour6bar ) return Colour(-c);
    return c;
  }

  /**
   * Return the flavour content of the given particle. The flavours
   * will be given in decreasing mass with flavour before
   * anti-flavour.
   */
  static vector<long> flavourContent(long id);

  /**
   * Return the flavour content of the given particle. The flavours
   * will be given in decreasing mass with flavour before
   * anti-flavour.
   */
  static vector<long> flavourContent(tcPDPtr);

  /**
   * Return the flavour content of the given particle. The flavours
   * will be given in decreasing mass with flavour before
   * anti-flavour.
   */
  static vector<long> flavourContent(tcPPtr);

  /**
   * Return the flavour content of the given particle. The flavours
   * will be given in decreasing mass with flavour before
   * anti-flavour.
   */
  static vector<long> flavourContent(const ParticleData &);

  /**
   * Return the flavour content of the given particle. The flavours
   * will be given in decreasing mass with flavour before
   * anti-flavour.
   */
  static vector<long> flavourContent(const Particle &);

};

/** Input a colour from a stream. */
template <typename IStream>
IStream & operator>>(IStream & is, PDT::Colour & c) {
  int ci;
  is >> ci;
  c = PDT::Colour(ci);
  return is;
}

/** Input a charge from a stream. */
template <typename IStream>
IStream & operator>>(IStream & is, PDT::Charge & c) {
  int ci;
  is >> ci;
  c = PDT::Charge(ci);
  return is;
}

/** Input a spin from a stream. */
template <typename IStream>
IStream & operator>>(IStream & is, PDT::Spin & s) {
  int si;
  is >> si;
  s = PDT::Spin(si);
  return is;
}

/// Type traits for built-in types
template <> 
struct TypeTraits<PDT::Spin>
{
  /** Enum for dimensions*/
  enum { hasDimension = false };
  /// Type switch set to standard type.
  typedef EnumT DimType;
  /// Base unit
  static constexpr PDT::Spin baseunit() { return PDT::Spin(1); }
};

/// Type traits for built-in types
template <> 
struct TypeTraits<PDT::Charge>
{
  /** Enum for dimensions*/
  enum { hasDimension = false };
  /// Type switch set to standard type.
  typedef EnumT DimType;
  /// Base unit
  static constexpr PDT::Charge baseunit() { return PDT::Charge(1); }
};

/// Type traits for built-in types
template <> 
struct TypeTraits<PDT::Colour>
{
  /** Enum for dimensions*/
  enum { hasDimension = false };
  /// Type switch set to standard type.
  typedef EnumT DimType;
  /// Base unit
  static constexpr PDT::Colour baseunit() { return PDT::Colour(3); }
};

}

#endif /* ThePEG_PDT_H */
