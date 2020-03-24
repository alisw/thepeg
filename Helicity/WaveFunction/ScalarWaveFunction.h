// -*- C++ -*-
//
// ScalarWaveFunction.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_ScalarWaveFunction_H
#define ThePEG_ScalarWaveFunction_H
//
// This is the declaration of the ScalarWaveFunction class.

#include "WaveFunctionBase.h"
#include <ThePEG/Helicity/ScalarSpinInfo.h>
#include <ThePEG/EventRecord/Particle.h>
#include <ThePEG/EventRecord/RhoDMatrix.h>

namespace ThePEG {
namespace Helicity {

/** \ingroup Helicity
 *  \author Peter Richardson
 * 
 *  This class is the base class for scalar wavefunctions for use in 
 *  helicity amplitude calculations. The general approach 
 *  is to use a similar philosophy to the FORTRAN HELAS code but with 
 *  additional structure.
 *
 *  This class stores the scalar wavefunction as a complex number and inherits
 *  from the WaveFunctionBase class for the storage of the particles
 *  momentum and type.
 * 
 *  @see WaveFunctionBase
 */
class ScalarWaveFunction : public WaveFunctionBase {

public:

  /** @name Standard constructors and destructors. */
  //@{

  /**
   * Constructor, set the momentum, direction and Wavefunction.
   * @param p The momentum.
   * @param part The ParticleData pointer
   * @param wave The wavefunction.
   * @param dir The direction of the particle.
   */
  ScalarWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,
		     Complex wave,Direction dir=intermediate) 
    : WaveFunctionBase(p,part,dir), _wf(wave)
  {
    assert(iSpin()==1);
  }

  /**
   * Constructor,set the 5-momentum and zero the wavefunction.
   * @param p The 5-momentum.
   * @param part The ParticleData pointer.
   * @param dir The direction of the particle.
   */
  ScalarWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,Direction dir) 
    : WaveFunctionBase(p,part,dir), _wf(1.0)
  {
    assert(iSpin()==1);
  }

  static void calculateWaveFunctions(RhoDMatrix & rho,
				     tPPtr, Direction) {
    rho=RhoDMatrix(PDT::Spin0);
  }

  static void constructSpinInfo(tPPtr part,Direction, bool time) {
    tScalarSpinPtr inspin;
    if(part->spinInfo()) inspin=dynamic_ptr_cast<tScalarSpinPtr>(part->spinInfo());
    if(inspin) return;
    assert(!part->spinInfo());
    ScalarSpinPtr temp = new_ptr(ScalarSpinInfo(part->momentum(),time));
    part->spinInfo(temp);
  }

  /**
   * Default constructor.
   */
  ScalarWaveFunction() : WaveFunctionBase(), _wf(1.0) {}

  /**
   *  Special for spin correlations
   */
  ScalarWaveFunction(tPPtr p,Direction dir,bool time) 
    : WaveFunctionBase(p->momentum(), p->dataPtr(), dir), _wf(1.0)
  {
    assert(iSpin()==1);
    constructSpinInfo(p,dir,time);
  }

  /**
   * Return the wavefunction.
   */
  const Complex & wave() const {return _wf;}

public:

  void transform(const LorentzRotation & r) {
    transformMomentum(r);
  }

private:

  /**
   * Complex number to store the wavefunction.
   */
  Complex _wf;

};
}
}

#endif
