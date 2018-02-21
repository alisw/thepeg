// -*- C++ -*-
//
// Lorentz5Vector.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Lorentz5Vector_H
#define ThePEG_Lorentz5Vector_H

// This is the declaration of the Lorentz5vector class.

#include "LorentzVector.h"
#include "Lorentz5Vector.fh"
#include "ThePEG/Utilities/Maths.h"
#include "ThePEG/Utilities/Direction.h"
#include "ThePEG/Utilities/UnitIO.h"
#include "LorentzRotation.h"

namespace ThePEG {

template <typename Value>
/**
 * The Lorentz5Vector inherits from the
 * <code>LorentzVector</code> class. It is templated on the
 * type of the member variables. The <code>Lorentz5Vector</code> class is a
 * <code>LorentzVector</code> with an extra member for the invariant
 * length/mass of the vector. Note that an object of the
 * <code>Lorentz5Vector</code> class may be internally inconsistent in
 * that the invariant length/mass of the <code>LorentzVector</code>
 * class need not be the same as the member variable representing the
 * invariant length/mass. The degree of inconsistency can be accessed
 * with the <code>massError()</code>, <code>energyError()</code> and
 * <code>rhoError()</code> methods and an object can be made consistent
 * using the <code>rescaleMass()</code>, <code>rescaleEnergy()</code> or
 * <code>rescaleRho()</code> methods.
 *
 * @see Math
 * 
 */
class Lorentz5Vector: public LorentzVector<Value> {

public:

  /** Template argument typedef. */
  typedef typename BinaryOpTraits<Value,Value>::MulT Value2;

public:
  /// Component access.
  //@{
  Value x() const { return LorentzVector<Value>::x(); }
  Value y() const { return LorentzVector<Value>::y(); }
  Value z() const { return LorentzVector<Value>::z(); }
  Value t() const { return LorentzVector<Value>::t(); }
  //@}

public:

  /** @name Constructors and destructor. */
  //@{
  /**
   * Constructor giving the null vector.
   */
  Lorentz5Vector() : mm() {}

  /**
   * Constructor giving the invariant length.
   */
  Lorentz5Vector(Value m) 
    : LorentzVector<Value>(Value(), Value(), Value(), m), mm(m) {}

  /**
   * Constructor giving the components x, y, z, t. The invariant
   * length is set to LorentzVector::mag().
   */
  Lorentz5Vector(Value x, Value y, Value z, Value t = Value())
    : LorentzVector<Value>(x, y, z, t) { rescaleMass(); }

  /**
   * Constructor giving the components x, y, z, t and invariant length.
   * May result in an inconsistent Lorentz5Vector.
   */
  Lorentz5Vector(Value x, Value y, Value z, Value t, Value tau)
    : LorentzVector<Value>(x, y, z, t), mm(tau) {}

  /**
   * Constructor giving a 3-Vector and a time component. The invariant
   * length is set to LorentzVector::mag().
   */
  Lorentz5Vector(const ThreeVector<Value> & p, Value e)
    : LorentzVector<Value>(p, e) { rescaleMass(); }

  /**
   * Constructor giving an invariant length and a 3-Vector
   * component. The time component is set to the corresponding value.
   */
  Lorentz5Vector(Value m, const ThreeVector<Value> & p) 
    : LorentzVector<Value>(p, sqrt(p.mag2() + m*m)), mm(m) {}

  /**
   * Constructor giving a 3-Vector, a time component and an invariant
   * length. May result in an inconsistent Lorentz5Vector.
   */
  Lorentz5Vector(const ThreeVector<Value> & p, Value t, Value tau)
    : LorentzVector<Value>(p, t), mm(tau) {}

  /**
   * Constructor giving a LorentzVector and an invariant length.
   * May result in an inconsistent Lorentz5Vector.
   */
  Lorentz5Vector(const LorentzVector<Value> & p, Value m) 
    : LorentzVector<Value>(p), mm(m) {}

  /**
   * Copy from HepLorentzVector constructor. The invariant
   * length is set to LorentzVector::mag().
   */
  Lorentz5Vector(const LorentzVector<Value> & p)
    : LorentzVector<Value>(p) { rescaleMass(); } 

  /**
   * Construct from value type U convertible to Value.
   */
  template<class U>
  Lorentz5Vector(const Lorentz5Vector<U> & p)
    : LorentzVector<Value>(p), mm(p.m) {}
  //@}

  /** @name Assignment and set functions. */
  //@{
  /**
   * Set invariant length/mass.
   */
  void setTau(Value a) { mm = a; }

  /**
   * Set invariant length/mass.
   */
  void setMass(Value a) { mm = a; }

  /**
   * Assignment. The invariant length is kept fixed. May result in an
   * inconsistent Lorentz5Vector.
   */
  Lorentz5Vector & operator=(const LorentzVector<Value> & q) {
    LorentzVector<Value>::operator=(q);
    return *this;
  }
  //@}

  /** @name Rescale functions to make consistent. */
  //@{
  /**
   * Rescale energy, so that the invariant length/mass of the
   * LorentzVector agrees with the current one.
   */
  void rescaleEnergy() {
    LorentzVector<Value>::setT(sqrt(LorentzVector<Value>::vect().mag2() + mass2()));
  }

  /**
   * Rescale spatial component, so that the invariant length/mass of
   * the LorentzVector agrees with the current one.
   */
  void rescaleRho() {
    LorentzVector<Value>::setRho(sqrt(t()*t() - mass2()));
  }

  /**
   * Set the invariant length/mass member, so that it agrees with the
   * invariant length/mass of the LorentzVector.
   */
  void rescaleMass() {
    mm = LorentzVector<Value>::m();
  }
  //@}

  /** @name Check consistency. */
  //@{
  /**
   * Return the relative inconsistency in the mass component.
   */
  double massError() const {
    return sqrt(abs(Math::relativeError(mass2(), 
					LorentzVector<Value>::m2())));
  }

  /**
   * Return the relative inconsistency in the energy component.
   */
  double energyError() const {
    return sqrt(abs(Math::relativeError(t()*t(), mass2() 
					+ LorentzVector<Value>::vect().mag2())));
  }

  /**
   * Return the relative inconsistency in the spatial components.
   */
  double rhoError() const {
    return sqrt(abs(Math::relativeError(LorentzVector<Value>::vect().mag2(), 
					t()*t() - mass2())));
  }
  //@}

  /** @name Access components. */
  //@{
  /**
   * Mass/invariant length component squared. m2() gives
   * the same calculated from the LorentzVector
   */
  Value2 mass2() const { return mm > Value() ? mm*mm: -mm*mm; }

  /**
   * Mass/invariant length component squared. m2() gives
   * the same calculated from the LorentzVector
   */
  Value2 tau2() const { return mass2(); }

  /**
   * Mass/invariant length component. m() gives the same
   * calculated from the LorentzVector
   */
  Value mass() const { return mm; }


  /**
   * Mass/invariant length component. m() gives the same
   * calculated from the LorentzVector
   */
  Value tau() const { return mass(); }

  /**
   * Return the positive negative light-cone components (depending on
   * the value of Direction<0>.
   */
  Value dirPlus() const {
    return Direction<0>::pos() ? 
      LorentzVector<Value>::plus() 
      : 
      LorentzVector<Value>::minus();
  }

  /**
   * Return the positive negative light-cone components (depending on
   * the value of Direction<0>.
   */
  Value dirMinus() const {
    return Direction<0>::neg() ? 
      LorentzVector<Value>::plus() 
      : 
      LorentzVector<Value>::minus();
  }
  //@}

  /**
   *  Perform a Lorentz transformation
   */
  Lorentz5Vector & transform(const LorentzRotation & r) 
  {
    LorentzVector<Value>::transform(r.one());
    return *this;
  }

private:

  /** The invariant mass/length member. */
  Value mm;

};

/** Output a Lorentz5Vector to a stream. */
template <typename OStream, typename T, typename UT>
void ounitstream(OStream & os, const Lorentz5Vector<T> & p, UT & u) {
  os << ounit(p.x(), u) << ounit(p.y(), u) << ounit(p.z(), u)
     << ounit(p.e(), u) << ounit(p.mass(), u);
}

/** Input a Lorentz5Vector from a stream. */
template <typename IStream, typename T, typename UT>
void iunitstream(IStream & is, Lorentz5Vector<T> & p, UT & u) {
  T x, y, z, e, mass;
  is >> iunit(x, u) >> iunit(y, u) >> iunit(z, u) >> iunit(e, u)
     >> iunit(mass, u);
  p = Lorentz5Vector<T>(x, y, z, e, mass);
}

template <typename T, typename U>
struct BinaryOpTraits;

/**
 *  Template for multiplication by scalar
 */
template <typename T>
struct BinaryOpTraits<Lorentz5Vector<T>, double> {
  /** The type resulting from multiplication of the template type with
      itself. */
  typedef Lorentz5Vector<T> MulT;
  /** The type resulting from division of one template type with
      another. */
  typedef Lorentz5Vector<T> DivT;
};

/**
 *  Template for multiplication by scalar
 */
template <typename U>
struct BinaryOpTraits<double, Lorentz5Vector<U> > {
  /** The type resulting from multiplication of the template type with
      itself. */
  typedef Lorentz5Vector<U> MulT;
};

/**
 * Template for multiplication for complex and Lorentz5Vector
 */
template <typename T, typename U>
struct BinaryOpTraits<Lorentz5Vector<T>, std::complex<U> > {
  /** The type resulting from multiplication of the template type with
      itself. */
  typedef Lorentz5Vector<std::complex<typename BinaryOpTraits<T,U>::MulT> > MulT;
  /** The type resulting from division of one template type with
      another. */
  typedef Lorentz5Vector<std::complex<typename BinaryOpTraits<T,U>::DivT> > DivT;
};

/**
 *  Template for multiplication by scalar
 */
template <typename T, typename U>
struct BinaryOpTraits<std::complex<T>, Lorentz5Vector<U> > {
  /** The type resulting from multiplication of the template type with
      itself. */
  typedef Lorentz5Vector<std::complex<typename BinaryOpTraits<T,U>::MulT> > MulT;
};

/**
 *  Template for multiplication of two Lorentz5Vectors
 */
template <typename T, typename U>
struct BinaryOpTraits<Lorentz5Vector<T>, Lorentz5Vector<U> > {
  /** The type resulting from multiplication of the template type with
      itself. */
  typedef typename BinaryOpTraits<T,U>::MulT MulT;
};

/**
 * Template for multiplication for LorentzVector and Lorentz5Vector
 */
template <typename T, typename U>
struct BinaryOpTraits<LorentzVector<T>, Lorentz5Vector<U> > {
  /** The type resulting from multiplication of the template type with
      itself. */
  typedef typename BinaryOpTraits<T,U>::MulT MulT;
};

/**
 * Template for multiplication for LorentzVector and Lorentz5Vector
 */
template <typename T, typename U>
struct BinaryOpTraits<Lorentz5Vector<T>, LorentzVector<U> > {
  /** The type resulting from multiplication of the template type with
      itself. */
  typedef typename BinaryOpTraits<T,U>::MulT MulT;
};

/// @name Dot product overloads.
//@{
template <typename ValueA, typename ValueB>
inline typename BinaryOpTraits<ValueA,ValueB>::MulT 
operator*(const Lorentz5Vector<ValueA> & a, const Lorentz5Vector<ValueB> & b) {
  return a.dot(b);
}

template <typename ValueA, typename ValueB>
inline typename BinaryOpTraits<ValueA,ValueB>::MulT 
operator*(const LorentzVector<ValueA> & a, const Lorentz5Vector<ValueB> & b) {
  return a.dot(b);
}

template <typename ValueA, typename ValueB>
inline typename BinaryOpTraits<ValueA,ValueB>::MulT 
operator*(const Lorentz5Vector<ValueA> & a, const LorentzVector<ValueB> & b) {
  return a.dot(b);
}

template <typename Value>
inline typename BinaryOpTraits<Value,Value>::MulT 
operator*(const Lorentz5Vector<Value> & a, const Lorentz5Vector<Value> & b) {
  return a.dot(b);
}
//@}
}

#endif /* ThePEG_Particle_H */
