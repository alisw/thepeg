// -*- C++ -*-
#ifndef Ariadne_PartonTraits_H
#define Ariadne_PartonTraits_H
// Here is some specializations of the ParticleTraits class.

#include "ThePEG/EventRecord/ParticleTraits.h"
#include "Ariadne/DipoleCascade/Parton.h"

namespace ThePEG {


/** @cond TRAITSPECIALIZATIONS */

/** Specialization of ParticleTraits for pointer to Particle. */
template <>
struct ParticleTraits<Ariadne::ParPtr>: public TraitsType {

  /**
   * Return a reference to the particle.
   */
  static Ariadne::Parton & ref(Ariadne::tParPtr p) {
    return *p;
  }

  /**
   * Return the momentum of particle \a p.
   */
  static const LorentzMomentum & momentum(Ariadne::tParPtr p) {
    return p->momentum();
  }

  /**
   * Return the mass of particle \a p.
   */
  static Energy mass(Ariadne::tParPtr p) {
    return p->momentum().mass();
  }

  /**
   * Perform a Lorentz transformation on particle \a p.
   */
  static void transform(Ariadne::tParPtr p, const LorentzRotation & r) {
    p->momentum().transform(r);
  }

  /**
   * Set the momentum and mass of a particle.
   */
  static void set5Momentum(Ariadne::tParPtr p, const Lorentz5Momentum & q) {
    p->momentum() = q;
  }

  /**
   * Set the 3-momentum of a particle. The energy is rescaled to
   * preserve invariant mass.
   */
  static void set3Momentum(Ariadne::tParPtr p, const Momentum3 & q) {
    p->momentum().setVect(q);
    p->momentum().rescaleEnergy();
  }

  /**
   * Return charge of particle \a p in units of e/3.
   */
  static int iCharge(Ariadne::tParPtr p) {
    return p->data().iCharge();
  }
};  

/** Specialization of ParticleTraits for pointer to const Particle. */
template <>
struct ParticleTraits<Ariadne::cParPtr>: public TraitsType {

  /**
   * Return a const reference to the particle.
   */
  static const Ariadne::Parton & ref(Ariadne::tcParPtr p) {
    return *p;
  }

  /**
   * Return the momentum of particle \a p.
   */
  static const LorentzMomentum & momentum(Ariadne::tcParPtr & p) {
    return p->momentum();
  }

  /**
   * Return the mass of particle \a p.
   */
  static Energy mass(Ariadne::tcParPtr p) {
    return p->momentum().mass();
  }

  /**
   * Return charge of particle \a p in units of e/3.
   */
  static int iCharge(Ariadne::tcParPtr & p) {
    return p->data().iCharge();
  }
};  

/** Specialization of ParticleTraits for transient pointer to Particle. */
template <>
struct ParticleTraits<Ariadne::tParPtr>: public TraitsType {

  /**
   * Return a reference to the particle.
   */
  static Ariadne::Parton & ref(Ariadne::tParPtr p) {
    return *p;
  }

  /**
   * Return the momentum of particle \a p.
   */
  static const LorentzMomentum & momentum(Ariadne::tParPtr p) {
    return p->momentum();
  }

  /**
   * Return the mass of particle \a p.
   */
  static Energy mass(Ariadne::tParPtr p) {
    return p->momentum().mass();
  }

  /**
   * Perform a Lorentz transformation on particle \a p.
   */
  static void transform(Ariadne::tParPtr p, const LorentzRotation & r) {
    p->momentum().transform(r);
  }

  /**
   * Set the momentum and mass of a particle.
   */
  static void set5Momentum(Ariadne::tParPtr p, const Lorentz5Momentum & q) {
    p->momentum() = q;
  }

  /**
   * Set the 3-momentum of a particle. The energy is rescaled to
   * preserve invariant mass.
   */
  static void set3Momentum(Ariadne::tParPtr p, const Momentum3 & q) {
    p->momentum().setVect(q);
    p->momentum().rescaleEnergy();
  }

  /**
   * Return charge of particle \a p in units of e/3.
   */
  static int iCharge(Ariadne::tParPtr p) {
    return p->data().iCharge();
  }
};  

/** Specialization of ParticleTraits for transient pointer to const Particle. */
template <>
struct ParticleTraits<Ariadne::tcParPtr>: public TraitsType {

  /**
   * Return a const reference to the particle.
   */
  static const Ariadne::Parton & ref(Ariadne::tcParPtr p) {
    return *p;
  }

  /**
   * Return the momentum of particle \a p.
   */
  static const LorentzMomentum & momentum(Ariadne::tcParPtr p) {
    return p->momentum();
  }

  /**
   * Return the mass of particle \a p.
   */
  static Energy mass(Ariadne::tcParPtr p) {
    return p->momentum().mass();
  }

  /**
   * Return charge of particle \a p in units of e/3.
   */
  static int iCharge(Ariadne::tcParPtr p) {
    return p->data().iCharge();
  }
};

/** @endcond */

}

#endif /* Ariadne_PartonTraits_H */
