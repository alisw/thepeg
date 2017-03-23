// -*- C++ -*-
//
// ScalarSpinInfo.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2011 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_ScalarSpinInfo_H
#define ThePEG_ScalarSpinInfo_H
// This is the declaration of the ScalarSpinInfo class.

#include "ThePEG/EventRecord/SpinInfo.h"
#include "ScalarSpinInfo.fh"

namespace ThePEG {
namespace Helicity {

/**
 * The ScalarSpinInfo class is designed to be the implementation of
 * the spin information for a scalar particle. Obviously it is pretty
 * trival in this case.
 *
 * @author Peter Richardson
 *
 */
class ScalarSpinInfo: public SpinInfo {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  ScalarSpinInfo() : SpinInfo(PDT::Spin0) {}

  /**
   * Standard Constructor.
   * @param p the production momentum.
   * @param time true if the particle is time-like.
   */
  ScalarSpinInfo(const Lorentz5Momentum & p, bool time) 
    : SpinInfo(PDT::Spin0, p, time) {}
  //@}

public:

  /**
   * Standard Init function.
   */
  static void Init();

  /**
   * Standard clone methods.
   */
  virtual EIPtr clone() const
  {
    tcSpinPtr temp = this;
    return const_ptr_cast<SpinPtr>(temp);
  }
  
  /**
   * Perform a lorentz rotation of the spin information
   */
  virtual void transform(const LorentzMomentum &,const LorentzRotation &);

private:

  /**
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<ScalarSpinInfo> initScalarSpinInfo;

  /**
   * Private and non-existent assignment operator.
   */
  ScalarSpinInfo & operator=(const ScalarSpinInfo &);

};

}
}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * ScalarSpinInfo.
 */
template <>
struct BaseClassTrait<ThePEG::Helicity::ScalarSpinInfo,1>
  : public ClassTraitsType {
  /** Typedef of the base class of ScalarSpinInfo. */
  typedef ThePEG::SpinInfo NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * ScalarSpinInfo class and the shared object where it is defined.
 */
template <>
struct ClassTraits<ThePEG::Helicity::ScalarSpinInfo>
  : public ClassTraitsBase<ThePEG::Helicity::ScalarSpinInfo> {
  /**
   * Return the class name.
   */
  static string className() { return "ThePEG::Helicity::ScalarSpinInfo"; }
};

/** @endcond */

}

#endif /* ThePEG_ScalarSpinInfo_H */
