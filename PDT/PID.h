// -*- C++ -*-
//
// PID.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_PID_H
#define ThePEG_PID_H

#include "EnumParticles.h"

namespace ThePEG {

/**
 * PID is a helper class implementing the type of PDG particle ids. At
 * the moment it just converts to a long, and other conversions are
 * intercepted.
 *
 * @see ParticleData
 */
class PID {
public:
  /**
   * Generic construction is disallowed.
   */
  template <typename T>
  PID(T t) { 
    t.ERROR_only_use_long_as_Particle_ID_type(); 
  }

  /**
   * Casting to generic types is disallowed.
   */
  template <typename T>
  operator T() const {
    T t; 
    t.ERROR_only_use_long_as_Particle_ID_type(); 
    return t;
  }

  /**
   * The negation operator
   */
  PID operator-() const;

private:
  /**
   * The particle id.
   */
  long id;
};

/// Specialized constructor for 'long'
template <> inline PID::PID(long t) : id(t) {} 

/// Specialized constructor for 'int'
template <> inline PID::PID(int t) : id(t) {} 

/// Specialized constructor for Particle Code enum
template <> inline PID::PID(ParticleID::ParticleCodes t) : id(t) {} 

/// Specialized cast for 'long'
template <> inline PID::operator long() const { return id; }

/**
 * The negation operator
 */
inline PID PID::operator-() const {
  return -id;
}

}

#endif /* ThePEG_PID_H */
