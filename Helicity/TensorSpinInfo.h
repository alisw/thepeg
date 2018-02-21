// -*- C++ -*-
//
// TensorSpinInfo.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_TensorSpinInfo_H
#define THEPEG_TensorSpinInfo_H
// This is the declaration of the TensorSpinInfo class.

#include "ThePEG/EventRecord/SpinInfo.h"
#include "ThePEG/Helicity/LorentzTensor.h"
#include "TensorSpinInfo.fh"
// #include "TensorSpinInfo.xh"
#include <array>

namespace ThePEG {
namespace Helicity {

/**
 *  The TensorSpinInfo class is the implementation of the spin
 *  information for tensor particles.  It inherits from the SpinInfo
 *  class and implements the storage of the basis tensors.
 *
 *  These basis states should be set by either matrix elements or
 *  decayers which are capable of generating spin correlation
 *  information.
 *
 *  The basis states in the rest frame of the particles can then be
 *  accessed by decayers to produce the correct correlation.
 *
 *  N.B. in our convention 0 is the \f$-2\f$ helicity state,
 *  1 is the \f$-1\f$ helicity state,
 *  2 is the \f$0\f$ helicity state,
 *  3 is the \f$+1\f$ helicity state and
 *  4 is the \f$+2\f$ helicity state.
 *
 * @author Peter Richardson
 *
 */
class TensorSpinInfo: public SpinInfo {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  TensorSpinInfo() : SpinInfo(PDT::Spin2), _decaycalc(false) {}
  
  /**
   * Standard Constructor.
   * @param p the production momentum.
   * @param time true if the particle is time-like.
   */
  TensorSpinInfo(const Lorentz5Momentum & p, bool time)
    : SpinInfo(PDT::Spin2, p, time), _decaycalc(false) {}
  //@}

public:

  /** @name Access the basis states. */
  //@{
  /**
   * Set the basis state, this is production state.
   * @param hel the helicity (0,1,2,3,4 as described above.)
   * @param in the LorentzTensor for the given helicity.
   */
  void setBasisState(unsigned int hel, LorentzTensor<double> in) const {
    assert(hel<5);
    _productionstates[hel]=in;
    _currentstates   [hel]=in;
  }

  /**
   * Set the basis state for the decay.
   * @param hel the helicity (0,1,2,3,4 as described above.)
   * @param in the LorentzTensor for the given helicity.
   */
  void setDecayState(unsigned int hel, LorentzTensor<double> in) const {
    assert(hel<5);
    _decaycalc = true;
    _decaystates[hel] = in;
  }

  /**
   * Get the basis state for the production for the given helicity, \a
   * hel  (0,1,2,3,4 as described above.)
   */
  const LorentzTensor<double> & getProductionBasisState(unsigned int hel) const {
    assert(hel<5);
    return _productionstates[hel];
  }

  /**
   * Get the basis state for the decay for the given helicity, \a hel
   * (0,1,2,3,4 as described above.)
   */
  const LorentzTensor<double> & getDecayBasisState(unsigned int hel) const {
    assert(hel<5);
    if(!_decaycalc) {
      for(unsigned int ix=0;ix<5;++ix)
	_decaystates[ix]=_currentstates[ix].conjugate();
      _decaycalc=true;
    }
    return _decaystates[hel];
  }
  //@}

  /**
   * Perform a lorentz rotation of the spin information
   */
  virtual void transform(const LorentzMomentum &,const LorentzRotation &);

public:

  /**
   * Standard Init function.
   */
  static void Init();

  /**
   * Standard clone method.
   */
  virtual EIPtr clone() const;

private:

  /**
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<TensorSpinInfo> initTensorSpinInfo;

  /**
   * Private and non-existent assignment operator.
   */
  TensorSpinInfo & operator=(const TensorSpinInfo &);

private:

  /**
   * Basis states in the frame in which the particle was produced.
   */
  mutable std::array<LorentzTensor<double>,5> _productionstates;

  /**
   * Basis states in the frame in which the particle decays.
   */
  mutable std::array<LorentzTensor<double>,5> _decaystates;

  /**
   * Basis states in the current frame of the particle
   */
  mutable std::array<LorentzTensor<double>,5> _currentstates;

  /**
   * True if the decay state has been set.
   */
  mutable bool _decaycalc;

};

}
}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * TensorSpinInfo.
 */
template <>
struct BaseClassTrait<ThePEG::Helicity::TensorSpinInfo,1>
  : public ClassTraitsType {
  /** Typedef of the base class of ScalarSpinInfo. */
  typedef ThePEG::SpinInfo NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * TensorSpinInfo class and the shared object where it is defined.
 */
template <>
struct ClassTraits<ThePEG::Helicity::TensorSpinInfo>
  : public ClassTraitsBase<ThePEG::Helicity::TensorSpinInfo> {
  /**
   * Return the class name.
   */
  static string className() { return "ThePEG::Helicity::TensorSpinInfo"; }
};

/** @endcond */

}

#endif /* THEPEG_TensorSpinInfo_H */
