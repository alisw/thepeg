// -*- C++ -*-
//
// ParticleTraits.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_ParticleTraits_H
#define ThePEG_ParticleTraits_H
// This is the declaration of the ParticleTraits class.

#include "ThePEG/Config/ThePEG.h"
// #include "ParticleTraits.fh"
// #include "ParticleTraits.xh"

namespace ThePEG {

template <typename PType>
/**
 * ParticleTraits is a templated class defining a general interface to
 * any particle class. To make another particle type
 * <code>PType</code> available to some general ThePEG routines, the
 * ParticleTraits should be specialized to that class implementing
 * relevant methods of the general ParticleTraits class
 * below. Typically one needs specialisation both for the class itself
 * and of pointers to the class.
 * 
 * @see Particle
 * @see Lorentz5Vector
 * 
 */
struct ParticleTraits: public TraitsType {

  /**
   * Return a reference to the particle.
   */
  static PType & ref(PType & p) {
    return p;
  }

  /**
   * Return the momentum of particle \a p.
   */
  static LorentzMomentum momentum(const PType & p) {
    return p.momentum();
  }

  /**
   * Return the mass of particle \a p.
   */
  static Energy mass(const PType & p) {
    return p.mass();
  }

  /**
   * Perform a Lorentz transformation on particle \a p.
   */
  static void transform(PType & p, const LorentzRotation & r) {
    p.transform(r);
  }

  /**
   * Set the momentum and mass of a particle.
   */
  static void set5Momentum(PType & p, const Lorentz5Momentum & q) {
    p.set5Momentum(q);
  }

  /**
   * Set the 3-momentum of a particle. The energy is rescaled to
   * preserve invariant mass.
   */
  static void set3Momentum(PType & p, const Momentum3 & q) {
    p.set3Momentum(q);
  }

  /**
   * Return charge of particle \a p in units of e/3.
   */
  static int iCharge(const PType & p) {
    return p.data().iCharge();
  }

};

/** @cond TRAITSPECIALIZATIONS */

/** Specialization of ParticleTraits for pointer to Particle. */
template <>
struct ParticleTraits<PPtr>: public TraitsType {

  /**
   * Return a reference to the particle.
   */
  static Particle & ref(tPPtr p) {
    return *p;
  }

  /**
   * Return the momentum of particle \a p.
   */
  static const LorentzMomentum & momentum(tPPtr p) {
    return p->momentum();
  }

  /**
   * Return the mass of particle \a p.
   */
  static Energy mass(tPPtr p) {
    return p->mass();
  }

  /**
   * Perform a Lorentz transformation on particle \a p.
   */
  static void transform(tPPtr p, const LorentzRotation & r) {
    p->transform(r);
  }

  /**
   * Set the momentum and mass of a particle.
   */
  static void set5Momentum(tPPtr p, const Lorentz5Momentum & q) {
    p->set5Momentum(q);
  }

  /**
   * Set the 3-momentum of a particle. The energy is rescaled to
   * preserve invariant mass.
   */
  static void set3Momentum(tPPtr p, const Momentum3 & q) {
    p->set3Momentum(q);
  }

  /**
   * Return charge of particle \a p in units of e/3.
   */
  static int iCharge(tPPtr p) {
    return p->data().iCharge();
  }
};  

/** Specialization of ParticleTraits for pointer to const Particle. */
template <>
struct ParticleTraits<cPPtr>: public TraitsType {

  /**
   * Return a const reference to the particle.
   */
  static const Particle & ref(tcPPtr p) {
    return *p;
  }

  /**
   * Return the momentum of particle \a p.
   */
  static const LorentzMomentum & momentum(tcPPtr & p) {
    return p->momentum();
  }

  /**
   * Return the mass of particle \a p.
   */
  static Energy mass(tcPPtr p) {
    return p->mass();
  }

  /**
   * Return charge of particle \a p in units of e/3.
   */
  static int iCharge(tcPPtr & p) {
    return p->data().iCharge();
  }
};  

/** Specialization of ParticleTraits for transient pointer to Particle. */
template <>
struct ParticleTraits<tPPtr>: public TraitsType {

  /**
   * Return a reference to the particle.
   */
  static Particle & ref(tPPtr p) {
    return *p;
  }

  /**
   * Return the momentum of particle \a p.
   */
  static const LorentzMomentum & momentum(tPPtr p) {
    return p->momentum();
  }

  /**
   * Return the mass of particle \a p.
   */
  static Energy mass(tPPtr p) {
    return p->mass();
  }

  /**
   * Perform a Lorentz transformation on particle \a p.
   */
  static void transform(tPPtr p, const LorentzRotation & r) {
    p->transform(r);
  }

  /**
   * Set the momentum and mass of a particle.
   */
  static void set5Momentum(tPPtr p, const Lorentz5Momentum & q) {
    p->set5Momentum(q);
  }

  /**
   * Set the 3-momentum of a particle. The energy is rescaled to
   * preserve invariant mass.
   */
  static void set3Momentum(tPPtr p, const Momentum3 & q) {
    p->set3Momentum(q);
  }

  /**
   * Return charge of particle \a p in units of e/3.
   */
  static int iCharge(tPPtr p) {
    return p->data().iCharge();
  }
};  

/** Specialization of ParticleTraits for transient pointer to const Particle. */
template <>
struct ParticleTraits<tcPPtr>: public TraitsType {

  /**
   * Return a const reference to the particle.
   */
  static const Particle & ref(tcPPtr p) {
    return *p;
  }

  /**
   * Return the momentum of particle \a p.
   */
  static const LorentzMomentum & momentum(tcPPtr p) {
    return p->momentum();
  }

  /**
   * Return the mass of particle \a p.
   */
  static Energy mass(tcPPtr p) {
    return p->mass();
  }

  /**
   * Return charge of particle \a p in units of e/3.
   */
  static int iCharge(tcPPtr p) {
    return p->data().iCharge();
  }
};

/** Partial specialization of ParticleTraits for bare pointer to anything. */
template <typename T>
struct ParticleTraits<T*>: public TraitsType {

  /**
   * Return a reference to the particle.
   */
  static Particle & ref(T * p) {
    return *p;
  }

  /**
   * Return the momentum of particle \a p.
   */
  static const LorentzMomentum & momentum(T * p) {
    return ParticleTraits<T>::momentum(*p);
  }

  /**
   * Return the mass of particle \a p.
   */
  static Energy mass(T * p) {
    return ParticleTraits<T>::mass(*p);
  }

  /**
   * Perform a Lorentz transformation on particle \a p.
   */
  static void transform(T * p, const LorentzRotation & r) {
    ParticleTraits<T>::transform(*p, r);
  }

  /**
   * Set the momentum and mass of a particle.
   */
  static void set5Momentum(T * p, const Lorentz5Momentum & q) {
    ParticleTraits<T>::set5Momentum(*p, q);
  }

  /**
   * Set the 3-momentum of a particle. The energy is rescaled to
   * preserve invariant mass.
   */
  static void set3Momentum(T * p, const Momentum3 & q) {
    ParticleTraits<T>::set3Momentum(*p, q);
  }

  /**
   * Return charge of particle \a p in units of e/3.
   */
  static int iCharge(T * p) {
    return ParticleTraits<T>::iCharge(*p);
  }
};  

/** Partial specialization of ParticleTraits for bare pointer to const
    anything. */
template <typename T>
struct ParticleTraits<const T *>: public TraitsType {

  /**
   * Return a const reference to the particle.
   */
  static const Particle & ref(const T * p) {
    return *p;
  }

  /**
   * Return the momentum of particle \a p.
   */
  static const LorentzMomentum & momentum(const T * p) {
    return ParticleTraits<T>::momentum(*p);
  }

  /**
   * Return the mass of particle \a p.
   */
  static Energy mass(const T * p) {
    return ParticleTraits<T>::mass(*p);
  }

  /**
   * Return charge of particle \a p in units of e/3.
   */
  static int iCharge(const T * p) {
    return ParticleTraits<T>::iCharge(*p);
  }
};

/** Specialization of ParticleTraits for LorentzMomentum. In this way
 *  a LorentzMomentum can be used where only the momentum parts of a
 *  Particle is required. */
template <>
struct ParticleTraits<LorentzMomentum>: public TraitsType {

  /**
   * Return a reference to the LorentzMomentum.
   */
  static LorentzMomentum & ref(LorentzMomentum & p) {
    return p;
  }

  /**
   * Return the momentum.
   */
  static const LorentzMomentum & momentum(const LorentzMomentum & p) {
    return p;
  }

  /**
   * Return the mass.
   */
  static Energy mass(const LorentzMomentum & p) {
    return p.m();
  }

  /**
   * Perform a Lorentz transformation.
   */
  static void transform(LorentzMomentum & p, const LorentzRotation & r) {
    p.transform(r);
  }

  /**
   * Set the momentum and mass.
   */
  static void set5Momentum(LorentzMomentum & p, const Lorentz5Momentum & q) {
    p = q;
  }

  /**
   * Set the 3-momentum. The energy is rescaled to
   * preserve invariant mass.
   */
  static void set3Momentum(LorentzMomentum & p, const Momentum3 & q) {
    p = LorentzMomentum(q, sqrt(q.mag2() + p.m2()));
  }
};  

/** Specialization of ParticleTraits for Lorentz5Momentum. In this way
 *  a Lorentz5Momentum can be used where only the momentum parts of a
 *  Particle is required. */
template <>
struct ParticleTraits<Lorentz5Momentum>: public TraitsType {

  /**
   * Return a reference to the Lorentz5Momentum.
   */
  static Lorentz5Momentum & ref(Lorentz5Momentum & p) {
    return p;
  }

  /**
   * Return the momentum.
   */
  static const LorentzMomentum & momentum(const Lorentz5Momentum & p) {
    return p;
  }

  /**
   * Return the mass.
   */
  static Energy mass(const Lorentz5Momentum & p) {
    return p.mass();
  }

  /**
   * Perform a Lorentz transformation.
   */
  static void transform(Lorentz5Momentum & p, const LorentzRotation & r) {
    p.transform(r);
  }

  /**
   * Set the momentum and mass.
   */
  static void set5Momentum(Lorentz5Momentum & p, const Lorentz5Momentum & q) {
    p = q;
  }

  /**
   * Set the 3-momentum. The energy is rescaled to
   * preserve invariant mass.
   */
  static void set3Momentum(Lorentz5Momentum & p, const Momentum3 & q) {
    p = Lorentz5Momentum(p.mass(), q);
  }
};  

/** @endcond */

/** A helper class to be used in <code>std::</code> algorithms to
 *  transform a range of particles. */
struct Transformer {
  /** Constructor taking a reference to the Lorentz rotation to be
   *  performed. */
  Transformer(const LorentzRotation & rin) : r(rin) {}
  /** Copy constructor. */
  Transformer(const Transformer & t) : r(t.r) {}
  /** Perform the rotation on a given particle. */
  template <typename PType>
  void operator()(const PType & p) {
    ParticleTraits<PType>::transform(p, r);
  }
  /** A reference to the Lorentz rotation to be performed. */
  const LorentzRotation & r;
};

}

// #include "ParticleTraits.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ParticleTraits.tcc"
#endif

#endif /* ThePEG_ParticleTraits_H */
