// -*- C++ -*-
//
// VectorSpinInfo.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_VectorSpinInfo_H
#define THEPEG_VectorSpinInfo_H
// This is the declaration of the VectorSpinInfo class.

#include "ThePEG/EventRecord/SpinInfo.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "VectorSpinInfo.fh"

namespace ThePEG {
namespace Helicity {

/**
 *  The VectorSpinInfo class is the implementation of the spin
 *  information for vector particles.  It inherits from the SpinInfo
 *  class and implements the storage of the basis vectors.
 *
 *  These basis states should be set by either matrixelements or
 *  decayers which are capable of generating spin correlation
 *  information.
 *
 *  The basis states in the rest frame of the particles can then be
 *  accessed by decayers to produce the correct correlation.
 *
 *  N.B. in our convention 0 is the \f$-1\f$ helicity state,
 *  1 is the \f$0\f$ helicity state and
 *  2 is the \f$+1\f$ helicity state.
 *
 * @author Peter Richardson
 *
 */
class VectorSpinInfo: public SpinInfo {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  VectorSpinInfo() : SpinInfo(PDT::Spin1), _decaycalc(false) {}
  
  /**
   * Standard Constructor.
   * @param p the production momentum.
   * @param time true if the particle is time-like.
   */
  VectorSpinInfo(const Lorentz5Momentum & p, bool time)
    : SpinInfo(PDT::Spin1, p, time), _decaycalc(false) {}
  //@}
  
public:

  /** @name Set and get methods for the basis state. */
  //@{
  /**
   * Set the basis state, this is production state.
   * @param hel the helicity (0,1,2 as described above.)
   * @param in the LorentzPolarizationVector for the given helicity.
   */
  void setBasisState(unsigned int hel, 
		     const LorentzPolarizationVector & in) const {
    assert(hel<3);
    _productionstates[hel] = in;
    _currentstates   [hel] = in;
  }

  /**
   * Set the basis state for the decay.
   * @param hel the helicity (0,1,2 as described above.)
   * @param in the LorentzPolarizationVector for the given helicity.
   */
  void setDecayState(unsigned int hel, 
		     const LorentzPolarizationVector & in) const {
    assert(hel<3);
    _decaycalc = true;
    _decaystates[hel] = in;;
  }

  /**
   * Get the basis state for the production for the given helicity, \a
   * hel (0,1,2 as described above.)
   */
  const LorentzPolarizationVector & getProductionBasisState(unsigned int hel) const {
    assert(hel<3);
    return _productionstates[hel];
  }

  /**
   * Get the basis state for the decay for the given helicity, \a hel 
   * (0,1,2 as described above.)
   */
  const LorentzPolarizationVector & getDecayBasisState(unsigned int hel) const {
    assert(hel<3);
    if(!_decaycalc) {
      for(unsigned int ix=0;ix<3;++ix)
	_decaystates[ix]=_currentstates[ix].conjugate();
      _decaycalc=true;
    }
    // return the basis function
    return _decaystates[hel];
  }
  //@}

  /**
   * Perform a Lorentz rotation of the spin information
   */
  virtual void transform(const LorentzMomentum &,const LorentzRotation & );

  /**
   *  Undecay
   */
  virtual void undecay() const {
    _decaycalc=false;
    SpinInfo::undecay();
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
  VectorSpinInfo & operator=(const VectorSpinInfo &) = delete;

private:

  /**
   * Basis states in the frame in which the particle was produced.
   */
  mutable std::array<LorentzPolarizationVector,3> _productionstates;

  /**
   * Basis states in the frame in which the particle decays.
   */
  mutable std::array<LorentzPolarizationVector,3> _decaystates;

  /**
   * Basis states in the current frame of the particle
   */
  mutable std::array<LorentzPolarizationVector,3> _currentstates;

  /**
   * True if the decay state has been set.
   */
  mutable bool _decaycalc;

};

}
}


namespace ThePEG {

}
#endif /* THEPEG_VectorSpinInfo_H */
