// -*- C++ -*-
//
// ScalarSpinInfo.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
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
   * Private and non-existent assignment operator.
   */
  ScalarSpinInfo & operator=(const ScalarSpinInfo &) = delete;

};

}
}


namespace ThePEG {

}
#endif /* ThePEG_ScalarSpinInfo_H */
