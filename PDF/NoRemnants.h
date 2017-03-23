// -*- C++ -*-
//
// NoRemnants.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_NoRemnants_H
#define ThePEG_NoRemnants_H
// This is the declaration of the NoRemnants class.

#include "ThePEG/PDF/RemnantHandler.h"
// #include "NoRemnants.fh"
// #include "NoRemnants.xh"

namespace ThePEG {

/**
 * NoRemnants inherits from RemnantHandler but can only handle
 * particles without sub-structure with the parton density given by a
 * NoPDF object and will never give any remnants.
 *
 *
 * @see RemnantHandler
 * @see NoPDF
 * 
 */
class NoRemnants: public RemnantHandler {

public:

  /** @name Virtual functions mandated by the RemnantHandler base class. */
  //@{
  /**
   * Return true if this remnant handler can handle extracting all
   * specified partons. The NoRemnants will return false if any
   * partons are given.
   */
  virtual bool canHandle(tcPDPtr, const cPDVector & partons) const {
    return partons.empty();
  }

  /**
   * Generate Remnants. Will not generate remnants and will throw a
   * RemnantHandlerException if the extracted parton is not the
   * incomining particle with x=1.
   */
  virtual Lorentz5Momentum generate(PartonBinInstance & pb, const double * r,
				    Energy2 scale,
				    const LorentzMomentum & p,
				    bool fixedPartonMomentum = false) const;

  /**
   * Generate the momentum of the extracted parton with the \a parent
   * momentum given by the last argument. If the \a scale is negative,
   * it means that the doScale in the previous call to nDim() was
   * true, otherwise the given \a scale should be the virtuality of
   * the extracted parton. \a shat is the total invariant mass squared
   * of the hard sub-system produced by the extracted parton and the
   * primary parton entering from the other side. Generated quantities
   * which are not returned in the momentum may be saved in the
   * PartonBinInstance, \a pb, for later use. In particular, if the
   * nDim() random numbers, \a r, are not enough to generate with
   * weight one, the resulting weight should be stored with the
   * remnantWeight() method of the parton bin.
   */
  virtual Lorentz5Momentum generate(PartonBinInstance & pb, const double * r,
				    Energy2 scale, Energy2 shat,
				    const LorentzMomentum & parent,
				    bool fixedPartonMomentum = false) const;
  //@}

public:

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<NoRemnants> initNoRemnants;

  /**
   *  Private and non-existent assignment operator.
   */
  NoRemnants & operator=(const NoRemnants &);

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of NoRemnants. */
template <>
struct BaseClassTrait<NoRemnants,1>: public ClassTraitsType {
  /** Typedef of the first base class of NoRemnants. */
  typedef RemnantHandler NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  NoRemnants class. */
template <>
struct ClassTraits<NoRemnants>: public ClassTraitsBase<NoRemnants> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::NoRemnants"; }
};

/** @endcond */

}

#endif /* ThePEG_NoRemnants_H */
