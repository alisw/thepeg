// -*- C++ -*-
//
// SimplePhaseSpace.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_SimplePhaseSpace_H
#define ThePEG_SimplePhaseSpace_H

#include "ThePEG/Config/ThePEG.h"

#include "ThePEG/Vectors/LorentzRotation.h"
#include "ThePEG/Vectors/LorentzRotation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/ParticleTraits.h"
#include "ThePEG/Repository/UseRandom.h"
#include "SimplePhaseSpace.xh"
#include <numeric>

namespace ThePEG {

/**
 * SimplePhaseSpace defines a set of static functions to be used for
 * distributing momenta evenly in phase space. In most cases pointers
 * and references to both particle and momentum objects can be used as
 * arguments as long as the ParticleTraits class is specialized
 * properly. When needed, random numbers are generated with the
 * generator given by the static UseRandom class.
 */
struct SimplePhaseSpace {

  /**
   * Set two momenta in their center of mass system. Their total
   * invariant mass squared is given by s, and their direction is
   * distributed isotropically.
   * @param s the total invariant mass squared.
   * @param p1 pointer or reference to the first momentum. Its
   * invariant mass will be preserved.
   * @param p2 pointer or reference to the second momentum. Its
   * invariant mass will be preserved.
   * @throw ImpossibleKinematics if the sum of the invariant masses was
   * larger than the given invariant mass (\f$\sqrt{s}\f$).
   */
  template <typename PType>
  static void CMS(Energy2 s, PType & p1, PType & p2);

  /**
   * Set two momenta in their center of mass system. Their total
   * invariant mass squared is given by s, and their direction is
   * given in terms of the polar and azimuth angle of the first
   * momenta.
   * @param s the total invariant mass squared.
   * @param p1 pointer or reference to the first momentum. Its
   * invariant mass will be preserved.
   * @param p2 pointer or reference to the second momentum. Its
   * invariant mass will be preserved.
   * @param cosTheta cosine of the azimuth angle of the first momentum.
   * @param phi azimuth angle of the first momentum.
   * @throw ImpossibleKinematics if the sum of the invariant masses was
   * larger than the given invariant mass (\f$\sqrt{s}\f$).
   */
  template <typename PType>
  static void CMS(PType & p1, PType & p2, Energy2 s,
		  double cosTheta, double phi);

  /**
   * Set two momenta in their center of mass system. Their total
   * invariant mass squared is given by s. The helper momentum p0 is
   * used so that afterwards \f$t=(p0-p1)^2\f$ and p1 has the azimuth
   * angle phi around p0.
   * @param p1 pointer or reference to the first momentum. Its
   * invariant mass will be preserved.
   * @param p2 pointer or reference to the second momentum. Its
   * invariant mass will be preserved.
   * @param s the total invariant mass squared.
   * @param t \f$=(p0-p1)^2\f$.
   * @param phi azimuth angle of the first momentum around p0.
   * @param p0 pointer or reference to an auxiliary momentum.
   * @throw ImpossibleKinematics if the sum of the invariant masses was
   * larger than the given invariant mass (\f$\sqrt{s}\f$).
   */
  template <typename PType>
  static void CMS(PType & p1, PType & p2, Energy2 s, Energy2 t, double phi,
		  const PType & p0);

  /**
   * Set two momenta in their center of mass system. Their total
   * invariant mass squared is given by s. p1 will be along the z-axis.
   * @param p1 pointer or reference to the first momentum. Its
   * invariant mass will be preserved.
   * @param p2 pointer or reference to the second momentum. Its
   * invariant mass will be preserved.
   * @param s the total invariant mass squared.
   * @throw ImpossibleKinematics if the sum of the invariant masses was
   * larger than the given invariant mass (\f$\sqrt{s}\f$).
   */
  template <typename PType>
  static void CMS(PType & p1, PType & p2, Energy2 s);

  /**
   * Set two momenta in their center of mass system. Their total
   * invariant mass squared is given by s. The first will be along the
   * z-axis.
   * @param p a pair of pointers or references to the two momenta. Their
   * invariant masses will be preserved.
   * @param s the total invariant mass squared.
   * @throw ImpossibleKinematics if the sum of the invariant masses was
   * larger than the given invariant mass (\f$\sqrt{s}\f$).
   */
  template <typename PPairType>
  static void CMS(const PPairType & p, Energy2 s)
  {
    CMS(*p.first, *p.second, s);
  }

  /**
   * Set three momenta in their center of mass system. Their total
   * invariant mass squared is given by s. The energy fraction of
   * particle p1(3) is x1(3) of the total energy and the angles of the
   * system is distributed isotropically.
   * @param p1 pointer or reference to the first momentum. Its
   * invariant mass will be preserved.
   * @param p2 pointer or reference to the second momentum. Its
   * invariant mass will be preserved.
   * @param p3 pointer or reference to the second momentum. Its
   * invariant mass will be preserved.
   * @param s the total invariant mass squared.
   * @param x1 the energy fraction \f$2e_1/\sqrt{s}\f$.
   * @param x3 the energy fraction \f$2e_3/\sqrt{s}\f$.
   * @throw ImpossibleKinematics if the sum of the invariant masses was
   * larger than the given invariant mass (\f$\sqrt{s}\f$).
   */
  template <typename PType>
  static void CMS(PType & p1, PType & p2, PType & p3, Energy2 s,
		  double x1, double x3);

  /**
   * Set three momenta in their center of mass system. Their total
   * invariant mass squared is given by s. The energy fraction of
   * particle p1(3) is x1(3) of the total energy. Particle p1 is
   * initially placed along the z-axis and particle p2 is given
   * azimuth angle phii. Then the system is then rotated with
   * theta and phi respectively.
   * @param p1 pointer or reference to the first momentum. Its
   * invariant mass will be preserved.
   * @param p2 pointer or reference to the second momentum. Its
   * invariant mass will be preserved.
   * @param p3 pointer or reference to the second momentum. Its
   * invariant mass will be preserved.
   * @param s the total invariant mass squared.
   * @param x1 the energy fraction \f$2e_1/\sqrt{s}\f$.
   * @param x3 the energy fraction \f$2e_3/\sqrt{s}\f$.
   * @param phii the azimuth angle of p2 around p1.
   * @param theta the polar angle of p1.
   * @param phi the azimuth angle of p1.
   * @throw ImpossibleKinematics if the sum of the invariant masses was
   * larger than the given invariant mass (\f$\sqrt{s}\f$).
   */
  template <typename PType>
  static void CMS(PType & p1, PType & p2, PType & p3, Energy2 s,
		  double x1, double x3, double phii = 0.0,
		  double theta = 0.0, double phi = 0.0);

  /**
   * Calculate the absolute magnitude of the momenta of two particles
   * with masses m1 and m2 when put in their CMS of total invariant
   * mass squared s.
   * @param s the total invariant mass squared.
   * @param m1 the mass of particle 1.
   * @param m2 the mass of particle 2.
   * @throw ImpossibleKinematics if the sum of the masses was
   * larger than the given invariant mass (\f$\sqrt{s}\f$).
   */
  static Energy getMagnitude(Energy2 s, Energy m1, Energy m2);

  /**
   * Return a three-vector given the absolute momentum, cos(theta) and
   * phi.
   * @param p the magnitude of the momentum.
   * @param costheta the cosine of the polar angle.
   * @param phi the azimuth angle.
   */
  static Momentum3 polar3Vector(Energy p, double costheta, double phi)
  {
    return Momentum3(p*sqrt(1.0 - sqr(costheta))*sin(phi),
		     p*sqrt(1.0 - sqr(costheta))*cos(phi),
		     p*costheta);
  }

  /**
   * Get a number of randomly distributed momenta.
   * Given a number specified invariant masses and a
   * total invariant mass m0, return corresponding four-momenta
   * randomly distributed according to phase space. 
   * @param m0 the
   * total invariant mass of the resulting momenta.
   * @param m a vector
   * of invariant masses of the resulting momenta.
   * @return a vector
   * of momenta with the given masses randomly distributed.
   * @throw ImpossibleKinematics if the sum of the masses was
   * larger than the given invariant mass (\f$\sqrt{s}\f$).
   */
  static vector<LorentzMomentum>
  CMSn(Energy m0, const vector<Energy> & m);

  /**
   * Set the momentum of a number of particles. Given a number of
   * particles and a total invariant mass m0, distribute their
   * four-momenta randomly according to phase space.
   * @param particles a container of particles or pointers to
   * particles. The invariant mass of these particles will not be
   * chaned.
   * @param m0 the
   * total invariant mass of the resulting momenta.
   * @throw ImpossibleKinematics if the sum of the masses was
   * larger than the given invariant mass (\f$\sqrt{s}\f$).
   */
  template <typename Container>
  static void CMSn(Container & particles, Energy m0);

};

}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "SimplePhaseSpace.tcc"
#endif

#endif /* ThePEG_SimplePhaseSpace_H */
