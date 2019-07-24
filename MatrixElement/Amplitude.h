// -*- C++ -*-
//
// Amplitude.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Amplitude_H
#define ThePEG_Amplitude_H
// This is the declaration of the Amplitude class.


#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/EventRecord/Particle.h"

namespace ThePEG {

/**
 * The Amplitude class is the abstract base class for all the classes
 * representing complex amplitudes associated with either a hard
 * 2\f$\rightarrow\f$ N subprocess or a decay 1\f$\rightarrow\f$ N
 * process.  The returned value should be dimensionless suitable
 * scaled by the total invariant mass squared (shat), which is always
 * computable from the specified momenta of the particles in the
 * vertex.  Notice that the amplitude for splitting
 * 1\f$\rightarrow\f$ N processes is instead represented in other
 * classes (derived from the SplitFun class).
 *
 * @see \ref AmplitudeInterfaces "The interfaces"
 * defined for Amplitude.
 */
class Amplitude: public HandlerBase {

  /** @name Main virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return the amplitude. Given the ParticleData objects in \a
   * particles, their \a momenta and \a helicities of all the
   * particles in the vertex, return the complex amplitude.  The
   * convention is the order of the vectors is that first there is the
   * incoming particle(s) and then the outgoing ones.  For the
   * helicities, the convention is to number them starting from 0 (no
   * negative values, because they are used as vector indeces), for
   * example, for a massive particle of spin S, 0 <= helicity <= 2*S.
   * The returned value should be dimensionless suitable scaled by the
   * total invariant mass squared (\f$\hat{s}\f$), which is always
   * computable from the specified \a momenta of the particles in the
   * vertex.
   */
  virtual Complex value(const tcPDVector & particles,
			const vector<Lorentz5Momentum> & momenta, 
			const vector<int> & helicities) = 0;

  /**
   * Return an overestimated amplitude. Same as value(const tcPDVector
   * &, const vector<Lorentz5Momentum> &, const vector<int> &), but it
   * provides an overestimate of the complex amplitude, that is:
   * <code>abs( overestimaValue() ) >= abs(value()) </code> The
   * default definition just returns value(), but it can be overriden
   * by derived classes.
   */
  virtual Complex overestimateValue(const tcPDVector & particles,
				    const vector<Lorentz5Momentum> & momenta, 
				    const vector<int> & helicities);
  //@}

  /** @name Alternative interface to main virtual functions. */
  //@{
  /**
   * Return the amplitude. Calls value(const tcPDVector &, const
   * vector<Lorentz5Momentum> &, const vector<int> &) and should not
   * be overridden.
   */

  Complex value(const PVector & particles, const vector<int> & helicities);

  /**
   * Return an overestimated amplitude. Calls overestimateValue(const
   * tcPDVector &, const vector<Lorentz5Momentum> &, const vector<int>
   * &)
   */
  Complex overestimateValue(const PVector & particles,
			    const vector<int> & helicities);
  //@}

public:

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

  /**
   * Describe an abstract base class with persistent data.
   */
  static AbstractNoPIOClassDescription<Amplitude> initAmplitude;

  /**
   *  Private and non-existent assignment operator.
   */
  Amplitude & operator=(const Amplitude &) = delete;

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the
 * base class of Amplitude.
 */
template <>
struct BaseClassTrait<Amplitude,1>: public ClassTraitsType {
  /** Typedef of the base class of Amplitude. */
  typedef HandlerBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * Amplitude class.
 */
template <>
struct ClassTraits<Amplitude>: public ClassTraitsBase<Amplitude> {
  /** Return the class name. */
  static string className() { return "ThePEG::Amplitude"; }
};

/** @endcond */

}

#endif /* ThePEG_Amplitude_H */
