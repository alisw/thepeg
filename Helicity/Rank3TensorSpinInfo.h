// -*- C++ -*-
//
// Rank3TensorSpinInfo.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_Rank3TensorSpinInfo_H
#define THEPEG_Rank3TensorSpinInfo_H
// This is the declaration of the Rank3TensorSpinInfo class.

#include "ThePEG/EventRecord/SpinInfo.h"
#include "ThePEG/Helicity/LorentzRank3Tensor.h"
#include "Rank3TensorSpinInfo.fh"
// #include "Rank3TensorSpinInfo.xh"
#include <array>

namespace ThePEG {
namespace Helicity {

/**
 *  The Rank3TensorSpinInfo class is the implementation of the spin
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
 *  N.B. in our convention 
 *  0 is the \f$-3\f$ helicity state,
 *  1 is the \f$-2\f$ helicity state,
 *  2 is the \f$-1\f$ helicity state,
 *  3 is the \f$0\f$ helicity state,
 *  4 is the \f$+1\f$ helicity state,
 *  5 is the \f$+2\f$ helicity state and
 *  6 is the \f$+3\f$ helicity state.
 *
 * @author Peter Richardson
 *
 */
class Rank3TensorSpinInfo: public SpinInfo {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  Rank3TensorSpinInfo() : SpinInfo(PDT::Spin3), _decaycalc(false) {}
  
  /**
   * Standard Constructor.
   * @param p the production momentum.
   * @param time true if the particle is time-like.
   */
  Rank3TensorSpinInfo(const Lorentz5Momentum & p, bool time)
    : SpinInfo(PDT::Spin3, p, time), _decaycalc(false) {}
  //@}

public:

  /** @name Access the basis states. */
  //@{
  /**
   * Set the basis state, this is production state.
   * @param hel the helicity (0,1,2,3,4 as described above.)
   * @param in the LorentzRank3Tensor for the given helicity.
   */
  void setBasisState(unsigned int hel, LorentzRank3Tensor<double> in) const {
    assert(hel<7);
    _productionstates[hel]=in;
    _currentstates   [hel]=in;
  }

  /**
   * Set the basis state for the decay.
   * @param hel the helicity (0,1,2,3,4 as described above.)
   * @param in the LorentzRank3Tensor for the given helicity.
   */
  void setDecayState(unsigned int hel, LorentzRank3Tensor<double> in) const {
    assert(hel<7);
    _decaycalc = true;
    _decaystates[hel] = in;
  }

  /**
   * Get the basis state for the production for the given helicity, \a
   * hel  (0,1,2,3,4 as described above.)
   */
  const LorentzRank3Tensor<double> & getProductionBasisState(unsigned int hel) const {
    assert(hel<7);
    return _productionstates[hel];
  }
  
  /**
   * Get the basis state for the current for the given helicity, \a
   * hel  (0,1,2,3,4 as described above.)
   */
  const LorentzRank3Tensor<double> & getCurrentBasisState(unsigned int hel) const {
    assert(hel<7);
    return _currentstates[hel];
  }

  /**
   * Get the basis state for the decay for the given helicity, \a hel
   * (0,1,2,3,4,5,6 as described above.)
   */
  const LorentzRank3Tensor<double> & getDecayBasisState(unsigned int hel) const {
    assert(hel<7);
    if(!_decaycalc) {
      for(unsigned int ix=0;ix<7;++ix)
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

  /**
   *  Undecay
   */
  virtual void undecay() const {
    _decaycalc=false;
    SpinInfo::undecay();
  }

  /**
   *  Reset
   */
  virtual void reset() {
    undecay();
    _currentstates = _productionstates;
    SpinInfo::reset();
  }

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
   * Private and non-existent assignment operator.
   */
  Rank3TensorSpinInfo & operator=(const Rank3TensorSpinInfo &) = delete;

private:

  /**
   * Basis states in the frame in which the particle was produced.
   */
  mutable std::array<LorentzRank3Tensor<double>,7> _productionstates;

  /**
   * Basis states in the frame in which the particle decays.
   */
  mutable std::array<LorentzRank3Tensor<double>,7> _decaystates;

  /**
   * Basis states in the current frame of the particle
   */
  mutable std::array<LorentzRank3Tensor<double>,7> _currentstates;

  /**
   * True if the decay state has been set.
   */
  mutable bool _decaycalc;

};

}
}


namespace ThePEG {

}
#endif /* THEPEG_Rank3TensorSpinInfo_H */
