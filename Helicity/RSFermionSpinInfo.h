// -*- C++ -*-
//
// RSFermionSpinInfo.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_RSFermionSpinInfo_H
#define THEPEG_RSFermionSpinInfo_H
// This is the declaration of the RSFermionSpinInfo class.

#include "ThePEG/EventRecord/SpinInfo.h"
#include "ThePEG/Helicity/LorentzRSSpinor.h"
#include "RSFermionSpinInfo.fh"
#include <array>

namespace ThePEG {
namespace Helicity {

/**
 *  The RSFermionSpinInfo class inherits from the SpinInfo class and
 *  implements the storage of the basis vector for a spin-3/2 particle.
 *  The basis states are the vector u spinors for a particle and the vector
 *  v-spinors for an antiparticle. The barred spinors can be obtained from these.
 *
 *  These basis states should be set by either the matrixelements or decayers
 *  which are capable of generating spin correlation information.
 *
 *  The basis states in the rest frame of the particles can then be accessed by
 *  the decayers to produce the correct correlations.
 *
 *  N.B. in our convention 0 is the \f$-\frac32\f$ helicity state,
 *  1 is the \f$-\frac12\f$ helicity state,
 *  2 is the \f$+\frac12\f$ helicity state,
 *  3 is the \f$+\frac32\f$ helicity state.
 *
 * @see SpinInfo
 *
 * \author Peter Richardson
 *
 */
class RSFermionSpinInfo: public SpinInfo {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  RSFermionSpinInfo()  : SpinInfo(PDT::Spin3Half), _decaycalc(false) {}

  /**
   * Standard Constructor.
   * @param p the production momentum.
   * @param time true if the particle is time-like.
   */
  RSFermionSpinInfo(const Lorentz5Momentum & p,bool time)
    : SpinInfo(PDT::Spin3Half, p, time), _decaycalc(false) {}
  //@}

public:

  /** @name Set and get methods for the basis state. */
  //@{
  /**
   * Set the basis state, this is production state.
   * @param hel the helicity (0,1,2,3 as described above.)
   * @param in the LorentzRSSpinor for the given helicity.
   */
  void setBasisState(unsigned int hel,
		     const LorentzRSSpinor<SqrtEnergy> & in) const {
    assert(hel<4);
    _productionstates[hel] = in;
    _currentstates   [hel] = in;
  }

  /**
   * Set the basis state for the decay.
   * @param hel the helicity (0,1,2,3 as described above.)
   * @param in the LorentzRSSpinor for the given helicity.
   */
  void setDecayState(unsigned int hel,
		     const LorentzRSSpinor<SqrtEnergy> & in) const {
    assert(hel<4);
    _decaycalc = true;
    _decaystates[hel] = in;
  }

  /**
   * Get the basis state for the production for the given helicity, \a
   * hel (0,1,2,3 as described above.)
   */
  const LorentzRSSpinor<SqrtEnergy> & getProductionBasisState(unsigned int hel) const {
    assert(hel<4);
    return _productionstates[hel];
  }
  
  /**
   * Get the basis state for the current for the given helicity, \a
   * hel (0,1,2,3 as described above.)
   */
  const LorentzRSSpinor<SqrtEnergy> & getCurrentBasisState(unsigned int hel) const {
    assert(hel<4);
    return _currentstates[hel];
  }

  /**
   * Get the basis state for the decay for the given helicity, \a hel
   * (0,1,2,3 as described above.)
   */
  const LorentzRSSpinor<SqrtEnergy> & getDecayBasisState(unsigned int hel) const {
    assert(hel<4);
    if(!_decaycalc) {
      for(unsigned int ix=0;ix<4;++ix) _decaystates[ix]=_currentstates[ix];
      _decaycalc=true;
    }
    return _decaystates[hel];
  }

  /**
   * Perform a lorentz rotation of the spin information
   */
  virtual void transform(const LorentzMomentum &,const LorentzRotation &);
  //@}

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
   * Standard Init function used to initialize the interfaces.
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
  RSFermionSpinInfo & operator=(const RSFermionSpinInfo &) = delete;

private:

  /**
   * Basis states in the frame in which the particle was produced.
   */
  mutable std::array<LorentzRSSpinor<SqrtEnergy>,4> _productionstates;

  /**
   * Basis states in the frame in which the particle decays.
   */
  mutable std::array<LorentzRSSpinor<SqrtEnergy>,4> _decaystates;

  /**
   * Basis states in the current frame of the particle
   */
  mutable std::array<LorentzRSSpinor<SqrtEnergy>,4> _currentstates;

  /**
   * True if the decay state has been set.
   */
  mutable bool _decaycalc;

};

}
}



namespace ThePEG {

}
#endif /* THEPEG_RSFermionSpinInfo_H */
