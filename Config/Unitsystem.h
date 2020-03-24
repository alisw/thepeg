// -*- C++ -*-
//
// Unitsystem.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad, David Grellscheid
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Units_H
#define ThePEG_Units_H

#include "ThePEG/Vectors/Lorentz5Vector.fh"
#include "ThePEG/Vectors/LorentzVector.fh"
#include "ThePEG/Vectors/ThreeVector.fh"
#include "ThePEG/Vectors/Transverse.fh"

#include "PhysicalQty.h"
#include "PhysicalQtyOps.h"
#include "PhysicalQtyComplex.h"


namespace ThePEG {

/**
 * The Units namespace contains the declaration of a number of classes
 * for variables with dimension. Currently they are all typedefs of
 * double, but in the future the SIUnits package will be used.
 *
 * The file Utilities/UnitIO.h defines helper-classes and helper
 * functions to read and write variables with dimensions. As an
 * example, to read and write an energy variable <code>e</code> in
 * units of GeV, use: <code>os << ounit(e, GeV)</code> and <code>is >>
 * iunit(e, GeV)</code>
*/
namespace Units {

/// adapter for the old style of naming quantities
template<long int L, long int E, long int Q, long int DL=1, long int DE=1, long int DQ=1>
using Qty = ThePEG::Qty<std::ratio<L,DL>, std::ratio<E,DE>, std::ratio<Q,DQ>>;

/** Energy. */
typedef Qty<0,1,0> Energy;

/** Mass has the same unit as Energy <=> c == 1. */
typedef Energy Mass;

/** Length. */
typedef Qty<1,0,0> Length;

/** Time has the same unit as Length. <=> c == 1. */
typedef Length Time;

/** Inverse Length. */
typedef Qty<-1,0,0> InvLength;

/** Velocities are dimensionless fractions of c. */
typedef double Velocity;

/** Charge. */
typedef Qty<0,0,1> Charge;

/** Angular momentum. */
typedef Qty<1,1,0> AngularMomentum;

/** Tension. */
typedef Qty<-1,1,0> Tension;

/** Area will be assumed to be Length\f$^2\f$. */
typedef Qty<2,0,0> Area;

/** Inverse Area. */
typedef Qty<-2,0,0> InvArea;

/** Cross section is an area. */
typedef Area CrossSection;

/**
 * @name Higher powers of energy.
 * Even higher powers can be created with similar typedefs.
 */
//@{
typedef Qty<0, 2, 0> Energy2;
typedef Qty<0, 3, 0> Energy3;
typedef Qty<0, 4, 0> Energy4;
typedef Qty<0, 5, 0> Energy5;
typedef Qty<0, 6, 0> Energy6;
typedef Qty<0, 7, 0> Energy7;
typedef Qty<0, 8, 0> Energy8;
typedef Qty<0, 9, 0> Energy9;
typedef Qty<0,10, 0> Energy10;
typedef Qty<0,11, 0> Energy11;
typedef Qty<0,12, 0> Energy12;

typedef Qty<0, 1,0, 1,2,1> SqrtEnergy;
typedef Qty<0,-1,0, 1,2,1> InvSqrtEnergy;

typedef Qty<0, -1, 0> InvEnergy;
typedef Qty<0, -2, 0> InvEnergy2;
typedef Qty<0, -3, 0> InvEnergy3;
typedef Qty<0, -4, 0> InvEnergy4;
typedef Qty<0, -5, 0> InvEnergy5;
typedef Qty<0, -6, 0> InvEnergy6;
typedef Qty<0, -7, 0> InvEnergy7;
typedef Qty<0, -8, 0> InvEnergy8;
typedef Qty<0, -9, 0> InvEnergy9;
typedef Qty<0,-10, 0> InvEnergy10;
typedef Qty<0,-11, 0> InvEnergy11;
typedef Qty<0,-12, 0> InvEnergy12;
//@}

/** CrossSection*Energy2. */
typedef Qty<2,2,0> Energy2XSec;

/** CrossSection/Energy2. */
typedef Qty<2,-2,0> DiffXSec;

/** CrossSection/Energy4. */
typedef Qty<2,-4,0> Diff2XSec;

/** CrossSection/Energy6 */
typedef Qty<2,-6,0> Diff3XSec;

/** Scale is the same as a squared energy. */
typedef Energy2 Scale;

/** A point in three-dimensional euclidean space. */
typedef ThreeVector<Length> Point;

/** A distance in three-dimensional euclidean space. */
typedef ThreeVector<Length> Distance;

/** A direction in three-dimensional euclidean space. */
typedef ThreeVector<double> Axis;

/** A momentum in three-dimensional euclidean space. */
typedef ThreeVector<Energy> Momentum3;

/** A three-dimensional boost vector. */
typedef ThreeVector<double> Boost;

/** A distance in four-dimensional space-time. */
typedef LorentzVector<Length> LorentzDistance;

/** A distance in four-dimensional space-time with an explicit
 *  invariant time component. */
typedef Lorentz5Vector<Length> Lorentz5Distance;

/** A point in four-dimensional space-time. */
typedef LorentzVector<Length> LorentzPoint;

/** A momentum in four-dimensional space-time. */
typedef LorentzVector<Energy> LorentzMomentum;

/** A momentum in four-dimensional space-time with an explicit
 *  invariant mass component. */
typedef Lorentz5Vector<Energy> Lorentz5Momentum;

/** Transverse components of a momentum. */
typedef Transverse<Energy> TransverseMomentum;

/// @name Pre-defined basic units.
//@{
constexpr Length operator "" _mm( long double x ) {
  return Length{Length::baseunit(), static_cast<double>(x)};
}
constexpr Length operator "" _mm( unsigned long long x ) {
  return Length{Length::baseunit(), static_cast<double>(x)};
}

constexpr Length meter      = 1.0e+3_mm;
constexpr Length millimeter = 1_mm;
constexpr Length mm         = 1_mm;
constexpr Length centimeter = 10_mm;
constexpr Length micrometer = 1.0e-3_mm;
constexpr Length nanometer  = 1.0e-6_mm;
constexpr Length picometer  = 1.0e-9_mm;
constexpr Length femtometer = 1.0e-12_mm;
 
constexpr Energy operator "" _MeV( long double x ) {
  return Energy{Energy::baseunit(), static_cast<double>(x)};
}
constexpr Energy operator "" _MeV( unsigned long long x ) {
  return Energy{Energy::baseunit(), static_cast<double>(x)};
}

constexpr Energy operator "" _GeV( long double x ) {
  return Energy{1000_MeV, static_cast<double>(x)};
}
constexpr Energy operator "" _GeV( unsigned long long x ) {
  return Energy{1000_MeV, static_cast<double>(x)};
}

constexpr Energy operator "" _TeV( long double x ) {
  return Energy{1000_GeV, static_cast<double>(x)};
}
constexpr Energy operator "" _TeV( unsigned long long x ) {
  return Energy{1000_GeV, static_cast<double>(x)};
}

constexpr Energy keV = 1.0e-3_MeV;
constexpr Energy MeV = 1_MeV;
constexpr Energy GeV = 1_GeV;
constexpr Energy TeV = 1_TeV;



constexpr Energy2 operator "" _MeV2( long double x ) {
  return Energy2{Energy2::baseunit(),	static_cast<double>(x)};
}
constexpr Energy2 operator "" _MeV2( unsigned long long x ) {
  return Energy2{Energy2::baseunit(),	static_cast<double>(x)};
}

constexpr Energy2 operator "" _GeV2( long double x ) {
  return Energy2{1.0e+6_MeV2, static_cast<double>(x)};
}
constexpr Energy2 operator "" _GeV2( unsigned long long x ) {
  return Energy2{1.0e+6_MeV2, static_cast<double>(x)};
}

constexpr Energy2 MeV2 = 1_MeV2;
constexpr Energy2 GeV2 = 1_GeV2;

constexpr InvEnergy InvGeV = 1/GeV;


constexpr Area operator "" _pb( long double x ) {
  return Area{1.0e-34 * Area::baseunit(), static_cast<double>(x)};
}
constexpr Area operator "" _pb( unsigned long long x ) {
  return Area{1.0e-34 * Area::baseunit(), static_cast<double>(x)};
}

constexpr Area femtobarn = 1.0e-03_pb;
constexpr Area picobarn  = 1_pb;
constexpr Area nanobarn  = 1.0e+03_pb;
constexpr Area microbarn = 1.0e+06_pb;
constexpr Area millibarn = 1.0e+09_pb;
constexpr Area barn      = 1.0e+12_pb; 

constexpr Charge eplus = Charge::baseunit();
//@}

/// Planck's constant times c (PDG 2006 value 197.326968(17) MeV fm)
constexpr Qty<1,1,0> hbarc = 197.326968e-15 * MeV * meter;
/// Planck's constant (PDG 2006 value 197.326968(17) MeV fm)
constexpr Qty<1,1,0> hbar_Planck = hbarc / 1.0; // c is one
}

/** 
 * Use symbols from this namespace to make forced breaks of unit
 * consistency explicit.
 */
namespace UnitRemoval {
  /// @name Helper units to make breaks of unit consistency explicit.
  //@{
  constexpr Units::Energy E = Units::Energy::baseunit();

  constexpr Units::Energy2 E2 = E*E;
  constexpr Units::Energy3 E3 = E*E2;
  constexpr Units::Energy4 E4 = E2*E2;

  constexpr Units::InvEnergy InvE = 1.0/E;
  constexpr Units::InvEnergy2 InvE2 = 1.0/E2;
  constexpr Units::InvEnergy3 InvE3 = 1.0/E3;
  constexpr Units::InvEnergy4 InvE4 = 1.0/E4;

  constexpr Units::SqrtEnergy SqrtE = Units::SqrtEnergy::baseunit();
  constexpr Units::InvSqrtEnergy InvSqrtE = Units::InvSqrtEnergy::baseunit();
  //@}
}

}

#endif /* ThePEG_Units_H */
