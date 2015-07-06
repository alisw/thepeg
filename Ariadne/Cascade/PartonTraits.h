// -*- C++ -*-
#ifndef Ariadne5_PartonTraits_H
#define Ariadne5_PartonTraits_H
// Here is some specializations of the ParticleTraits class.

#include "ThePEG/EventRecord/ParticleTraits.h"
#include "Parton.h"

namespace ThePEG {


/** @cond TRAITSPECIALIZATIONS */

/** Specialization of ParticleTraits for pointer to Parton. */
template <>
struct ParticleTraits<Ariadne5::ParPtr>: public TraitsType {

  /**
   * Return a reference to the particle.
   */
  static Ariadne5::Parton & ref(Ariadne5::tParPtr p) {
    return *p;
  }

  /**
   * Return the momentum of particle \a p.
   */
  static const LorentzMomentum & momentum(Ariadne5::tParPtr p) {
    return p->momentum();
  }

  /**
   * Return the mass of particle \a p.
   */
  static Energy mass(Ariadne5::tParPtr p) {
    return p->momentum().mass();
  }

  /**
   * Perform a Lorentz transformation on particle \a p.
   */
  static void transform(Ariadne5::tParPtr p, const LorentzRotation & r) {
    p->momentum().transform(r);
  }

  /**
   * Set the momentum and mass of a particle.
   */
  static void set5Momentum(Ariadne5::tParPtr p, const Lorentz5Momentum & q) {
    p->momentum() = q;
  }

  /**
   * Set the 3-momentum of a particle. The energy is rescaled to
   * preserve invariant mass.
   */
  static void set3Momentum(Ariadne5::tParPtr p, const Momentum3 & q) {
    p->momentum().setVect(q);
    p->momentum().rescaleEnergy();
  }

  /**
   * Return charge of particle \a p in units of e/3.
   */
  static int iCharge(Ariadne5::tParPtr p) {
    return p->data().iCharge();
  }
};  

/** Specialization of ParticleTraits for pointer to const Parton. */
template <>
struct ParticleTraits<Ariadne5::cParPtr>: public TraitsType {

  /**
   * Return a const reference to the particle.
   */
  static const Ariadne5::Parton & ref(Ariadne5::tcParPtr p) {
    return *p;
  }

  /**
   * Return the momentum of particle \a p.
   */
  static const LorentzMomentum & momentum(Ariadne5::tcParPtr & p) {
    return p->momentum();
  }

  /**
   * Return the mass of particle \a p.
   */
  static Energy mass(Ariadne5::tcParPtr p) {
    return p->momentum().mass();
  }

  /**
   * Return charge of particle \a p in units of e/3.
   */
  static int iCharge(Ariadne5::tcParPtr & p) {
    return p->data().iCharge();
  }
};  

/** Specialization of ParticleTraits for transient pointer to Parton. */
template <>
struct ParticleTraits<Ariadne5::tParPtr>: public TraitsType {

  /**
   * Return a reference to the particle.
   */
  static Ariadne5::Parton & ref(Ariadne5::tParPtr p) {
    return *p;
  }

  /**
   * Return the momentum of particle \a p.
   */
  static const LorentzMomentum & momentum(Ariadne5::tParPtr p) {
    return p->momentum();
  }

  /**
   * Return the mass of particle \a p.
   */
  static Energy mass(Ariadne5::tParPtr p) {
    return p->momentum().mass();
  }

  /**
   * Perform a Lorentz transformation on particle \a p.
   */
  static void transform(Ariadne5::tParPtr p, const LorentzRotation & r) {
    p->momentum().transform(r);
  }

  /**
   * Set the momentum and mass of a particle.
   */
  static void set5Momentum(Ariadne5::tParPtr p, const Lorentz5Momentum & q) {
    p->momentum() = q;
  }

  /**
   * Set the 3-momentum of a particle. The energy is rescaled to
   * preserve invariant mass.
   */
  static void set3Momentum(Ariadne5::tParPtr p, const Momentum3 & q) {
    p->momentum().setVect(q);
    p->momentum().rescaleEnergy();
  }

  /**
   * Return charge of particle \a p in units of e/3.
   */
  static int iCharge(Ariadne5::tParPtr p) {
    return p->data().iCharge();
  }
};  

/** Specialization of ParticleTraits for transient pointer to const Parton. */
template <>
struct ParticleTraits<Ariadne5::tcParPtr>: public TraitsType {

  /**
   * Return a const reference to the particle.
   */
  static const Ariadne5::Parton & ref(Ariadne5::tcParPtr p) {
    return *p;
  }

  /**
   * Return the momentum of particle \a p.
   */
  static const LorentzMomentum & momentum(Ariadne5::tcParPtr p) {
    return p->momentum();
  }

  /**
   * Return the mass of particle \a p.
   */
  static Energy mass(Ariadne5::tcParPtr p) {
    return p->momentum().mass();
  }

  /**
   * Return charge of particle \a p in units of e/3.
   */
  static int iCharge(Ariadne5::tcParPtr p) {
    return p->data().iCharge();
  }
};

/** @endcond */

}

#endif /* Ariadne5_PartonTraits_H */
